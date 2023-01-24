import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import os
import glob

from collections import Counter, defaultdict
from tqdm import tqdm
from scipy import stats
from scipy.stats import pearsonr
import random
import seaborn as sns

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

pd.set_option("display.max_columns", None)

def save_figures(f, exts=['png', 'pdf']):
    for ext in exts:
        plt.savefig(f + f'.{ext}', dpi=300, bbox_inches='tight', transparent=True)

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

cpg_type = 'opensea'
assert cpg_type in ['opensea', 'island', 'shelf_shore']

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

np.random.seed(42)
random.seed(42)

save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/tissue-specificity-bdm'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()

repr_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm'

T_colnames = [f'{chrom}_T' for chrom in CHR_LIST]
N_colnames = [f'{chrom}_N' for chrom in CHR_LIST]

T_repr_num_df = pd.DataFrame(index = NORMAL7_COHORT, columns = T_colnames)
N_repr_num_df = pd.DataFrame(index = NORMAL7_COHORT, columns = N_colnames)

for cohort in NORMAL7_COHORT:
    for chrom in CHR_LIST:
        cohort_repr_dir = os.path.join(repr_dir, cohort)
        globals()[f'{cohort}_{chrom}_T'] = glob.glob(os.path.join(cohort_repr_dir, f'repr_T_vector_{chrom}_*.npy'))
        globals()[f'{cohort}_{chrom}_N'] = glob.glob(os.path.join(cohort_repr_dir, f'repr_N_vector_{chrom}_*.npy'))
        
        T_repr_num_df.loc[cohort][f'{chrom}_T'] = len(glob.glob(os.path.join(cohort_repr_dir, f'repr_T_vector_{chrom}_*.npy')))
        N_repr_num_df.loc[cohort][f'{chrom}_N'] = len(glob.glob(os.path.join(cohort_repr_dir, f'repr_N_vector_{chrom}_*.npy')))

T_criterion = 18 
cohort_few_T = []
for cohort in T_repr_num_df.index.values:
    if T_repr_num_df.loc[cohort]['chr1_T'] < T_criterion:
        cohort_few_T.append(cohort)
        
N_criterion = 4
cohort_few_N = []
for cohort in N_repr_num_df.index.values:
    if N_repr_num_df.loc[cohort]['chr1_N'] < N_criterion:
        cohort_few_N.append(cohort)

print(f'cohrot having less than {T_criterion} representative T vectors: {cohort_few_T}')
print(f'cohrot having less than {N_criterion} representative N vectors: {cohort_few_N}')

NORMAL7_T_COHORT = list(set(NORMAL7_COHORT)-set(cohort_few_T))
NORMAL7_N_COHORT = list(set(NORMAL7_COHORT)-set(cohort_few_N))

print(f"cohorts for tumor-tumor comparison: \n {NORMAL7_T_COHORT}")
print(f"cohorts for normal-normal comparison: \n {NORMAL7_N_COHORT}")

for cohort in NORMAL7_T_COHORT:
    for chrom in CHR_LIST:
        globals()[f'{cohort}_{chrom}_T_sampled'] = random.sample(globals()[f'{cohort}_{chrom}_T'], min(T_criterion, T_repr_num_df.loc[cohort][f'{chrom}_T']))

for cohort in NORMAL7_N_COHORT:
    for chrom in CHR_LIST:
        globals()[f'{cohort}_{chrom}_N_sampled'] = random.sample(globals()[f'{cohort}_{chrom}_N'], min(N_criterion, N_repr_num_df.loc[cohort][f'{chrom}_N']))

repr_T_sampled = {}
for cohort in NORMAL7_COHORT:
    if cohort not in cohort_few_T:
        for chrom in CHR_LIST:
            k = f'{cohort}_{chrom}_T_sampled'
            repr_T_sampled[k] = globals()[f'{cohort}_{chrom}_T_sampled']
np.savez(os.path.join(save_dir, 'T_sampled-reprs'), **repr_T_sampled)
print(os.path.join(save_dir, 'T_sampled-reprs.npz'))

repr_N_sampled = {}
for cohort in NORMAL7_COHORT:
    if cohort not in cohort_few_N:
        for chrom in CHR_LIST:
            k = f'{cohort}_{chrom}_N_sampled'
            repr_N_sampled[k] = globals()[f'{cohort}_{chrom}_N_sampled']

np.savez(os.path.join(save_dir, 'N_sampled-reprs'), **repr_N_sampled)
print(os.path.join(save_dir, 'N_sampled-reprs.npz'))

print("Tumor-Tumor comparison")
type_full = 'Tumor'
type_ = 'T'
pcc_key = []
pval_key = []
for i in range(len(NORMAL7_T_COHORT)):
    target_cohort = NORMAL7_T_COHORT[i]
    print(f"===\nTarget: {target_cohort}")

    target_key = target_cohort.split('-')[1]

    for j in range(i, len(NORMAL7_T_COHORT)):

        probe_cohort = NORMAL7_T_COHORT[j]

        probe_key = probe_cohort.split('-')[1]

        for chrom in CHR_LIST:
            target_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{target_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
            probe_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{probe_cohort}_diffmat_bins.npz')[f'{chrom}_bins']

            globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'] = []
            globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'] = []

            if target_key != probe_key:                        
                intersecting_bins = np.intersect1d(target_bins, probe_bins)

                target_mask = [x in intersecting_bins for x in target_bins]
                probe_mask = [x in intersecting_bins for x in probe_bins]

            for k in range(len(globals()[f'{target_cohort}_{chrom}_{type_}_sampled'])):
                current_target_fname = globals()[f'{target_cohort}_{chrom}_{type_}_sampled'][k]
                current_target = np.load(current_target_fname).mean(axis = 0)

                for p in range(len(globals()[f'{probe_cohort}_{chrom}_{type_}_sampled'])):
                    current_probe_fname = globals()[f'{probe_cohort}_{chrom}_{type_}_sampled'][p]
                    current_probe = np.load(current_probe_fname).mean(axis = 0)

                    if target_key != probe_key:
                        masked_target = current_target[target_mask] 
                        masked_probe = current_probe[probe_mask]
                    else:
                        masked_target = current_target.copy()
                        masked_probe = current_probe.copy()

                    assert len(masked_target) == len(masked_probe)

                    pcc, pval = pearsonr(masked_target, masked_probe)
                    globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'].append(pcc)
                    globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'].append(pval)  
                    current_pcc_key = f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'
                    current_pval_key = f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'
                    if current_pcc_key not in pcc_key:
                        pcc_key.append(current_pcc_key)
                    if current_pval_key not in pval_key:
                        pval_key.append(current_pval_key)

print(f'num_cohorts for Tumor-Tumor comparison: {len(NORMAL7_T_COHORT)}')
print(f'num_all_pcc: {len(pcc_key)}')
print(f'num_all_pval: {len(pval_key)}')

type_full = 'Tumor'
globals()[f'{type_full}_pcc_all'] = {}
globals()[f'{type_full}_pval_all'] = {}
for x in pcc_key:
    globals()[f'{type_full}_pcc_all'][x] = globals()[x]
for y in pval_key:
    globals()[f'{type_full}_pval_all'][y] = globals()[y]
    
print(f"num_total_pcc ({type_full}): {len(globals()[f'{type_full}_pcc_all'])}")
print(f"num_total_pval ({type_full}): {len(globals()[f'{type_full}_pval_all'])}")

np.savez(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pcc'), **globals()[f'{type_full}_pcc_all'])
np.savez(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pval'), **globals()[f'{type_full}_pval_all'])

print(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pcc.npz'))
print(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pval.npz'))

type_full = 'Tumor'
for chrom in CHR_LIST:
    globals()[f'{type_full}_{chrom}_pcc_key'] = []
    globals()[f'{type_full}_{chrom}_pval_key'] = []
    for x in pcc_key:
        if chrom+'_' in x:
            globals()[f'{type_full}_{chrom}_pcc_key'].append(x)
    for y in pval_key:
        if chrom+'_' in y:
            globals()[f'{type_full}_{chrom}_pval_key'].append(y)

type_ = 'T'
type_full = 'Tumor'
indices = [x.split('-')[1] for x in NORMAL7_T_COHORT]
for chrom in CHR_LIST:
    globals()[f'{type_}_{chrom}_pcc_df'] = pd.DataFrame(index = indices, columns = indices)
    globals()[f'{type_}_{chrom}_num_sig_pval_df'] = pd.DataFrame(index = indices, columns = indices)
    
    for i in range(len(indices)):
        target = indices[i]
        for j in range(i, len(indices)):
            probe = indices[j]
            
            current_pval_key = f'{type_full}_{target}_{probe}_{chrom}_pval'
            if current_pval_key not in pval_key:
                current_pval_key = f'{type_full}_{probe}_{target}_{chrom}_pval'
            tmp = globals()[current_pval_key]
            pval_mask = [x <= 5e-2 for x in tmp]

            current_pcc_key = current_pval_key.replace('pval', 'pcc')
            tmp2 = globals()[current_pcc_key]
            sig_pcc = np.array(tmp2)[pval_mask].copy()
            sig_pcc_mean = sig_pcc.mean()

            globals()[f'{type_}_{chrom}_pcc_df'].loc[target][probe] = float(sig_pcc_mean)
            globals()[f'{type_}_{chrom}_pcc_df'].loc[probe][target] = float(sig_pcc_mean)

            globals()[f'{type_}_{chrom}_num_sig_pval_df'].loc[target][probe] = int(np.array(pval_mask).sum())
            globals()[f'{type_}_{chrom}_num_sig_pval_df'].loc[probe][target] = int(np.array(pval_mask).sum())
    
    pcc_fname = f'{type_}_{chrom}_pcc_df.csv'
    pval_fname = f'{type_}_{chrom}_num_sig_pval_df'
    
    globals()[f'{type_}_{chrom}_pcc_df'].to_csv(os.path.join(save_dir, pcc_fname))
    print(os.path.join(save_dir, pcc_fname))
    
    globals()[f'{type_}_{chrom}_num_sig_pval_df'].to_csv(os.path.join(save_dir, pval_fname))
    print(os.path.join(save_dir, pval_fname))

type_full = 'Normal'
type_ = 'N'
pcc_key = []
pval_key = []

print(f"{type_full}-{type_full} comparison")

for i in range(len(NORMAL7_N_COHORT)):
    target_cohort = NORMAL7_N_COHORT[i]
    print(f"===\nTarget: {target_cohort}")

    target_key = target_cohort.split('-')[1]

    for j in range(i, len(NORMAL7_N_COHORT)):

        probe_cohort = NORMAL7_N_COHORT[j]

        probe_key = probe_cohort.split('-')[1]

        for chrom in CHR_LIST:
            target_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{target_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
            probe_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{probe_cohort}_diffmat_bins.npz')[f'{chrom}_bins']

            globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'] = []
            globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'] = []

            if target_key != probe_key:                        
                intersecting_bins = np.intersect1d(target_bins, probe_bins)

                target_mask = [x in intersecting_bins for x in target_bins]
                probe_mask = [x in intersecting_bins for x in probe_bins]

            for k in range(len(globals()[f'{target_cohort}_{chrom}_{type_}_sampled'])):
                current_target_fname = globals()[f'{target_cohort}_{chrom}_{type_}_sampled'][k]
                current_target = np.load(current_target_fname).mean(axis = 0)

                for p in range(len(globals()[f'{probe_cohort}_{chrom}_{type_}_sampled'])):
                    current_probe_fname = globals()[f'{probe_cohort}_{chrom}_{type_}_sampled'][p]
                    current_probe = np.load(current_probe_fname).mean(axis = 0)

                    if target_key != probe_key:
                        masked_target = current_target[target_mask] 
                        masked_probe = current_probe[probe_mask]
                    else:
                        masked_target = current_target.copy()
                        masked_probe = current_probe.copy()

                    assert len(masked_target) == len(masked_probe)

                    pcc, pval = pearsonr(masked_target, masked_probe)
                    globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'].append(pcc)
                    globals()[f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'].append(pval)  
                    current_pcc_key = f'{type_full}_{target_key}_{probe_key}_{chrom}_pcc'
                    current_pval_key = f'{type_full}_{target_key}_{probe_key}_{chrom}_pval'
                    if current_pcc_key not in pcc_key:
                        pcc_key.append(current_pcc_key)
                    if current_pval_key not in pval_key:
                        pval_key.append(current_pval_key)

print(f'num_cohorts for {type_full}-{type_full} comparison: {len(NORMAL7_N_COHORT)}')
print(f'num_all_pcc: {len(pcc_key)}')
print(f'num_all_pval: {len(pval_key)}')

assert type_full == 'Normal'
globals()[f'{type_full}_pcc_all'] = {}
globals()[f'{type_full}_pval_all'] = {}
for x in pcc_key:
    globals()[f'{type_full}_pcc_all'][x] = globals()[x]
for y in pval_key:
    globals()[f'{type_full}_pval_all'][y] = globals()[y]
    
print(f"num_total_pcc ({type_full}): {len(globals()[f'{type_full}_pcc_all'])}")
print(f"num_total_pval ({type_full}): {len(globals()[f'{type_full}_pval_all'])}")

assert type_full == 'Normal'
np.savez(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pcc'), **globals()[f'{type_full}_pcc_all'])
np.savez(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pval'), **globals()[f'{type_full}_pval_all'])

print(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pcc.npz'))
print(os.path.join(save_dir, f'all_{cpg_type}_{type_full}_pval.npz'))

assert type_full == 'Normal'

for chrom in CHR_LIST:
    globals()[f'{type_full}_{chrom}_pcc_key'] = []
    globals()[f'{type_full}_{chrom}_pval_key'] = []
    for x in pcc_key:
        if chrom+'_' in x:
            globals()[f'{type_full}_{chrom}_pcc_key'].append(x)
    for y in pval_key:
        if chrom+'_' in y:
            globals()[f'{type_full}_{chrom}_pval_key'].append(y)

assert type_full == 'Normal'
type_ = type_full[0]

indices = [x.split('-')[1] for x in NORMAL7_N_COHORT]
for chrom in CHR_LIST:
    globals()[f'{type_}_{chrom}_pcc_df'] = pd.DataFrame(index = indices, columns = indices)
    globals()[f'{type_}_{chrom}_num_sig_pval_df'] = pd.DataFrame(index = indices, columns = indices)
    
    for i in range(len(indices)):
        target = indices[i]
        for j in range(i, len(indices)):
            probe = indices[j]
            
            current_pval_key = f'{type_full}_{target}_{probe}_{chrom}_pval'
            if current_pval_key not in pval_key:
                current_pval_key = f'{type_full}_{probe}_{target}_{chrom}_pval'
            tmp = globals()[current_pval_key]
            pval_mask = [x <= 5e-2 for x in tmp]

            current_pcc_key = current_pval_key.replace('pval', 'pcc')
            tmp2 = globals()[current_pcc_key]
            sig_pcc = np.array(tmp2)[pval_mask].copy()
            sig_pcc_mean = sig_pcc.mean()

            globals()[f'{type_}_{chrom}_pcc_df'].loc[target][probe] = float(sig_pcc_mean)
            globals()[f'{type_}_{chrom}_pcc_df'].loc[probe][target] = float(sig_pcc_mean)

            globals()[f'{type_}_{chrom}_num_sig_pval_df'].loc[target][probe] = int(np.array(pval_mask).sum())
            globals()[f'{type_}_{chrom}_num_sig_pval_df'].loc[probe][target] = int(np.array(pval_mask).sum())
    
    pcc_fname = f'{type_}_{chrom}_pcc_df.csv'
    pval_fname = f'{type_}_{chrom}_num_sig_pval_df'
    
    globals()[f'{type_}_{chrom}_pcc_df'].to_csv(os.path.join(save_dir, pcc_fname))
    print(os.path.join(save_dir, pcc_fname))
    
    globals()[f'{type_}_{chrom}_num_sig_pval_df'].to_csv(os.path.join(save_dir, pval_fname))
    print(os.path.join(save_dir, pval_fname))

tissue_specificity_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/tissue-specificity-bdm'

T_pcc_files = glob.glob(os.path.join(tissue_specificity_dir, f'T_*_pcc_df.csv'))
N_pcc_files = glob.glob(os.path.join(tissue_specificity_dir, f'N_*_pcc_df.csv'))
print(len(T_pcc_files))
print(len(N_pcc_files))

type_full = 'Tumor'
fig = plt.figure(figsize = (3*6, 3*4))
mins = []
maxs = []

for i in range(len(CHR_LIST)):
    chrom = CHR_LIST[i]
    ax = fig.add_subplot(4, 6, i+1)
    fnames = glob.glob(os.path.join(tissue_specificity_dir, f'{type_full[0]}*{chrom}_pcc_df.csv'))
    assert len(fnames) == 1
    fname = fnames[0]
    df = pd.read_csv(fname, index_col = 0)
    cohorts = np.array([x for x in df.index.values])
      
    ax = sns.heatmap(df, cmap = 'Reds')
    ax.set_title(chrom, fontsize = 13)
    
    mins.append(df.min().min())
    maxs.append(df.max().max())

fig.tight_layout()
fig_name = f'{cpg_type}_{type_full}_heatmaps.png'
plt.savefig(os.path.join(save_dir, fig_name))
print(os.path.join(save_dir, fig_name))

type_full = 'Normal'
fig = plt.figure(figsize = (3*6, 3*4))
mins = []
maxs = []

for i in range(len(CHR_LIST)):
    chrom = CHR_LIST[i]
    ax = fig.add_subplot(4, 6, i+1)
    fnames = glob.glob(os.path.join(tissue_specificity_dir, f'{type_full[0]}*{chrom}_pcc_df.csv'))
    assert len(fnames) == 1
    fname = fnames[0]
    df = pd.read_csv(fname, index_col = 0)
    cohorts = np.array([x for x in df.index.values])
    
    ax = sns.heatmap(df, cmap = 'Reds')
    ax.set_title(chrom, fontsize = 13)
    
    mins.append(df.min().min())
    maxs.append(df.max().max())

fig.tight_layout()
fig_name = f'{cpg_type}_{type_full}_heatmaps.png'
plt.savefig(os.path.join(save_dir, fig_name))
print(os.path.join(save_dir, fig_name))
