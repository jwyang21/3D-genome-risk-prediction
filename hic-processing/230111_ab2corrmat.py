import sys
import numpy as np
import fanc

ab_file = sys.argv[1]
output_prefix = ab_file[:-3]

AB_mat = fanc.load(ab_file)

CHROMS = [f'chr{i}' for i in range(1, 23)]
for chrom in CHROMS:
    ab_chr = AB_mat.matrix((chrom, chrom))
    np.save(f'{output_prefix}.{chrom}.corrmat', ab_chr.data)