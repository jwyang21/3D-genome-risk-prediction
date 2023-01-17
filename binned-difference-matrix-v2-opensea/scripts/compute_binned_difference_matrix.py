import argparse
import pandas as pd
import numpy as np

from collections import defaultdict

CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    parser.add_argument('-b', '--binsize', type=int, default=1e6)
    parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    parser.add_argument('-o', '--output', help='Output.', required=True)

    return parser.parse_args()

def compute_binned_matrix(beta_df, chrom, chrom2size, binsize, n_min_cpgs):
    # Get chromosome-wise methylation profile.
    beta_df_chrom = beta_df[beta_df.chrom == chrom].reset_index(drop=True).copy()

    # Assign bin labels.
    beta_df_chrom['bin'] = beta_df_chrom['start'] // binsize

    # Bin-CpG list map.
    bin2cpgs = defaultdict(list)
    for r in beta_df_chrom.to_records():
        bin2cpgs[r['bin']].append(r.beta)

    n_bins = int(chrom2size[chrom] // binsize) + 1

    binned_diffmat = np.zeros([n_bins, n_bins])  # Binned difference matrix.
    mask = np.zeros(n_bins).astype(bool)  # If # CpGs are not enough for bin i, True, otherwise False.

    for i in range(n_bins):
        cpgs_i = bin2cpgs[i]

        if len(cpgs_i) < n_min_cpgs:
            # Not enough CpGs in the bin.
            mask[i] = True
            continue

        for j in range(i, n_bins):
            cpgs_j = bin2cpgs[j]
            if len(cpgs_j) < n_min_cpgs:
                # Not enough CpGs in the bin.
                continue

            median_i = np.median(cpgs_i)
            median_j = np.median(cpgs_j)
            binned_diffmat[i, j] = binned_diffmat[j, i] = np.abs(median_i - median_j)

    return binned_diffmat, mask

if __name__ == '__main__':
    args = parse_arguments()

    chrom2size = pd.read_csv(args.chrom_size, sep='\t', names=['chrom', 'size'])
    chrom2size = {r.chrom:r['size'] for r in chrom2size.to_records()}

    result = {}

    beta_df = pd.read_csv(args.input, sep='\t', names=['chrom', 'start', 'end', 'beta'])
    for chrom in CHROMOSOMES:
        binned_diffmat, mask = compute_binned_matrix(beta_df, chrom, chrom2size, args.binsize, args.n_min_cpgs)

        result[f'{chrom}'] = binned_diffmat
        result[f'{chrom}_mask'] = mask

    np.savez(args.output, **result)
        
