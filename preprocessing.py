import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pairwise', type=str, default='/fs/nexus-scratch/rhaworth/plasmid/pairs.txt')
parser.add_argument('--kmers', type=str, default='/fs/nexus-projects/plasmids/plasmid_host_range_2024_summer_interns/expanded_dataset/4mer_table.txt')
parser.add_argument('--srccs', type=str, default='/fs/nexus-projects/plasmids/plasmid_host_range_2024_summer_interns/expanded_dataset/correlation_table_v4.txt')
parser.add_argument('-k', '--kmer_feats', choices=['none', 'rmsd', 'profile'])
parser.add_argument('-s', '--use_srcc', action='store_true')
parser.add_argument('-o', '--outfile', type=str, default='/fs/nexus-scratch/rhaworth/plasmid/processed.csv')
args = parser.parse_args()

data = pd.read_csv(args.pairwise, sep='\t')

# compute pairwise features
data['Len_ratio'] = data['Chr_len'].div(data['Plasmid_len'])
data['GC_abs_diff'] = data['Chr_GC'].add(data['Plasmid_GC'] * -1).abs()

# fetch kmers
if args.kmer_feats is not 'none':
    df_kmers = pd.read_csv(args.kmers, sep='\t')
    # split into chromosomes and plasmids
    c_idx = df_kmers.index[df_kmers['Seq_type'] == 'chromosome']
    c_kmers = df_kmers.iloc[c_idx]
    p_kmers = df_kmers.iloc[~c_idx]

    # remove unfound labels, get c_kmers in order of data
    labels_c = data['Assembly_ID'].str[:-7]
    labels_kept_c = labels_c.isin(c_kmers['SeqID'])
    labels_c = labels_c[labels_kept_c]

    c_kmers = c_kmers.set_index('SeqID').loc[labels_c]

    # repeat for plasmids; drop rows that have been dropped already
    labels_p = data['Plasmid_ID']
    labels_p = labels_p[labels_kept_c]
    labels_kept_p = labels_p.isin(p_kmers['SeqID'])
    labels_p = labels_p[labels_kept_p]

    p_kmers = p_kmers.set_index('SeqID').loc[labels_p]

    # get corresponding subsets of original data
    data_c = data[labels_kept_c]
    data_p = data_c[labels_kept_p]

    # isolate kmers in df then compute kmer frequencies for plasmids then chromosomes
    p_kmers['SeqID'] = p_kmers.index
    p_kmers['Index'] = data_p.index
    p_kmers = p_kmers.set_index('Index').iloc[:, 2:-1]
    p_sums = p_kmers.sum(axis=1)
    p_freqs = p_kmers.div(p_sums, axis=0)

    c_kmers['Index'] = data_c.index
    c_kmers = c_kmers.set_index('Index')[labels_kept_p]
    c_kmers = c_kmers.iloc[:, 2:]
    c_sums = c_kmers.sum(axis=1)
    c_freqs = c_kmers.div(c_sums, axis=0)

    # compute differences
    diffs = p_freqs - c_freqs

    # if using full profile, get abs distance and concat
    if args.kmer_feats == 'profile':
        data = pd.concat([data_p, diffs.abs()])
    # otherwise, compute RMSD
    else:
        rmsd = diffs ** 2
        rmsd = rmsd.sum(axis=1)
        rmsd = rmsd / 136 # * 1/n
        rmsd = rmsd ** 0.5

        # copy, add column
        data = data_p.copy()
        data['Kmer_RMSD'] = rmsd

# fetch SRCCs
if args.use_srcc:
    # add column, get trimmed assembly IDs
    data['SRCC'] = 0.0
    assemblies_trimmed = data['Assembly_ID'].str[:-7].values

    with open(args.srccs, 'r') as f:
        # get list of plasmids
        plasmid_list = f.readline().split()
        # iterate over rest of file
        index = 0
        while True:
            if index < 6800: # found in practice the first 6800 lines aren't present in data, so this makes things faster
                next(f)
                index += 1
                continue
            line = f.readline()
            if line is None or line == '':
                break
            assembly = line.split('\t', 1)[0]
            assembly = assembly[:-15]
            # check for containment; assume all assemblies found
            if assembly in assemblies_trimmed:
                rows_found = data[data['Assembly_ID'].str[:-7] == assembly]
                for i, plasmid in enumerate(rows_found['Plasmid_ID']):
                    if plasmid not in plasmid_list:
                        # drop row if plasmid not found
                        data.drop(rows_found.iloc[i].name, inplace=True)
                    else:
                        # if both assembly and plasmid are found, split line further to extract plasmid
                        plasmid_idx = plasmid_list.index(plasmid)
                        # if idx=1, we need 3 splits: 1 for assembly ID, 2 for ID + idx 0, 3 for ID + idx 0-1
                        # idx+1 indexing works in last index edge case, -2 indexing doesn't
                        srcc = line.split('\t', plasmid_idx+2)[plasmid_idx+1]
                        # update dataframe
                        data.loc[rows_found.iloc[i].name, 'SRCC'] = float(srcc)

# write
data.to_csv(args.outfile)
print('wrote to', args.outfile)