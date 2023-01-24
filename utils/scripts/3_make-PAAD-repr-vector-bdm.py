import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]
cohort = 'TCGA-PAAD'
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
all_npzs = glob.glob(os.path.join(cohort_pc1_dir, f'*.npz')) 
bdm_pc1_files = []

for i in range(len(all_npzs)):
    if 'inv_exp' not in all_npzs[i]:
        bdm_pc1_files.append(all_npzs[i])
        
bdm_N_pc1_files = []

for x in bdm_pc1_files:
    for n in N:
        if n in x:
            bdm_N_pc1_files.append(x)

for chrom in CHR_LIST:
    k = f'{chrom}_pc1'
    globals()[f'{chrom}_array'] = []
    for i in range(len(bdm_N_pc1_files)):
        current_file = bdm_N_pc1_files[i]
        globals()[f'{chrom}_array'].append(np.load(current_file)[k])
    globals()[f'{chrom}_array'] = np.vstack(globals()[f'{chrom}_array'])
    result_fname = f'repr_N_vector_{chrom}_0'
    np.save(os.path.join(save_dir, result_fname), globals()[f'{chrom}_array'])
    print(os.path.join(save_dir, result_fname)+'.npy')
