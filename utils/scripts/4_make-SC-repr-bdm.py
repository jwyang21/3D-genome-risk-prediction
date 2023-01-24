import pandas as pd
import numpy as np
import glob
import os
import random
random.seed(42)
np.random.seed(42)
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]
cohort = 'PCBC'
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'
cpg_type = 'opensea'

def get_sample_list(cohort):
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = []
    else: 
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: 
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S

T, N, S = get_sample_list(cohort)

cohort_pc1_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/1_compute-score-{cpg_type}/result/{cohort}/pc1'
save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm/{cohort}'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

all_npzs = glob.glob(os.path.join(cohort_pc1_dir, f'*.npz')) 
bdm_pc1_files = []
for i in range(len(all_npzs)):
    if 'inv_exp' not in all_npzs[i]:
        bdm_pc1_files.append(all_npzs[i])

random.shuffle(bdm_pc1_files)

repr_size = 10
num_repr = (len(bdm_pc1_files) // repr_size)
save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm/{cohort}'

for chrom in CHR_LIST:
    print(f"===\n{chrom}")
    k = f'{chrom}_pc1'
    print(f"total {num_repr} repr vectors")
    for i in range(num_repr):
        print(f"---\nMaking repr vector {i}")
        data = []
        for j in range(repr_size):
            current_pointer = i * repr_size + j
            current_pc1_fname = bdm_pc1_files[current_pointer]
            current_pc1 = np.load(current_pc1_fname)[k]
            data.append(current_pc1)
        data = np.vstack(data)
        print(data.shape)
        fname = f'repr_S_vector_{chrom}_{i}'
        np.save(os.path.join(save_dir, fname), data)
        print(os.path.join(save_dir, fname))
        del(data)
        print("---")
