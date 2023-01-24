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

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
pd.set_option("display.max_columns", None)

def save_figures(f, exts=['png', 'pdf']):
    for ext in exts:
        plt.savefig(f + f'.{ext}', dpi=300, bbox_inches='tight', transparent=True)

SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

cpg_type = 'opensea'
assert cpg_type in ['opensea', 'island', 'shelf_shore']

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

np.random.seed(42)

def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()

tcga_fire_cohort = pd.read_csv('/data/project/3dith/data/etc/tcga-fire-cohorts.csv', index_col = 0)
tcga_cohort_info = pd.read_csv('/data/project/3dith/data/etc/manifest.csv', index_col = 0)
save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm'
cohort_few_normal = ['TCGA-CHOL', 'TCGA-ESCA', 'TCGA-PAAD']

for chrom in CHR_LIST:
    print(f"==\nchrom: {chrom}")
    for cohort in NORMAL7_COHORT:
        all_npzs = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*.npz') 
        files = []
        for i in range(len(all_npzs)):
            if 'inv_exp' not in all_npzs[i]:
                files.append(all_npzs[i])
        if chrom == CHR_LIST[0]:
            if not os.path.exists(os.path.join(save_dir, cohort)):
                os.makedirs(os.path.join(save_dir, cohort))
        data_T = []
        data_N = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            if int(sample[13:15])<=9:
                data_T.append(standardize(np.load(f)[f'{chrom}_pc1']))
            else:
                data_N.append(standardize(np.load(f)[f'{chrom}_pc1']))
        data_T = np.vstack(data_T)
        data_N = np.vstack(data_N)

        idx_T = np.arange(len(data_T))
        np.random.shuffle(idx_T)

        idx_N = np.arange(len(data_N))
        np.random.shuffle(idx_N)

        data_T = data_T[idx_T]
        data_N = data_N[idx_N]

        num_repr_T_vectors = (len(data_T) // 10)
        num_repr_N_vectors = (len(data_N) // 10)

        for i in range(num_repr_T_vectors):
            current_T_repr = data_T[(i) * 10 : (i+1)*10]
            assert len(current_T_repr) == 10
            np.save(f'{save_dir}/{cohort}/repr_T_vector_{chrom}_{str(i)}.npy', current_T_repr)

        for i in range(num_repr_N_vectors):
            if cohort not in cohort_few_normal:
                current_N_repr = data_N[(i) * 10 : (i+1)*10]
                assert len(current_N_repr) == 10
                np.save(f'{save_dir}/{cohort}/repr_N_vector_{chrom}_{str(i)}.npy', current_T_repr)
