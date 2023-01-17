#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#%run imports.ipynb


# In[1]:


# contents in imports.ipynb
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

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)

plt.rc('figure', figsize=(1.5, 1.5))

pd.set_option("display.max_columns", None)

def save_figures(f, exts=['png', 'pdf']):
    for ext in exts:
        plt.savefig(f + f'.{ext}', dpi=300, bbox_inches='tight', transparent=True)


# ## global variables

# In[2]:


NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[3]:


cpg_type = 'opensea'
assert cpg_type in ['opensea', 'island', 'shelf_shore']


# In[4]:


CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


# In[5]:


# fix random seed
np.random.seed(42)


# ## functions

# In[6]:


def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()


# In[7]:


tcga_fire_cohort = pd.read_csv('/data/project/3dith/data/etc/tcga-fire-cohorts.csv', index_col = 0)
# tcga_fire_cohort.loc['TCGA-LIHC','FIRE'] # -> 'LI'
# print(tcga_fire_cohort.columns) #Index(['FIRE', 'TCGA_disease', 'FIRE_tissue'], dtype='object')



tcga_cohort_info = pd.read_csv('/data/project/3dith/data/etc/manifest.csv', index_col = 0)





cpg_type='opensea'


# In[12]:


save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors'

cohort_few_normal = ['TCGA-CHOL', 'TCGA-ESCA', 'TCGA-PAAD']


# 각 cohort마다, tumor/normal group 각각을 random shuffle해서 index list를 얻은 후, 맨 앞에서부터 10개씩 잘라서 tumor/normal representative vector를 만든다. 

# In[ ]:


for chrom in CHR_LIST:
    print(f"==\nchrom: {chrom}")
    for cohort in NORMAL7_COHORT:
        files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*_inv_exp.npz') 
        if chrom == CHR_LIST[0]:
            #print(f"{len(files)} samples in {cohort}")
            
            if not os.path.exists(os.path.join(save_dir, cohort)):
                os.makedirs(os.path.join(save_dir, cohort))
            
        # make representative vector
        data_T = []
        data_N = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            #print(int(sample[13:15]) <= 9)
            if int(sample[13:15])<=9:
                data_T.append(standardize(np.load(f)[f'{chrom}_pc1']))
            else:
                data_N.append(standardize(np.load(f)[f'{chrom}_pc1']))

            #print(len(data_T), len(data_N))
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
            #print(f'{save_dir}/{cohort}/repr_T_vector_{chrom}_{str(i)}.npy')

        for i in range(num_repr_N_vectors):
            if cohort not in cohort_few_normal:
                current_N_repr = data_N[(i) * 10 : (i+1)*10]
                assert len(current_N_repr) == 10
                np.save(f'{save_dir}/{cohort}/repr_N_vector_{chrom}_{str(i)}.npy', current_T_repr)
                #print(f'{save_dir}/{cohort}/repr_N_vector_{chrom}_{str(i)}.npy')

