#!/usr/bin/env python
# coding: utf-8

# In[12]:


import pandas as pd
import numpy as np
import glob
import os


# In[2]:


import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)


# In[3]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]


# In[4]:


cohort = 'TCGA-PAAD'


# In[5]:
# 230119_make-TCGA-repr-bdm.py만 돌리면 PAAD의 repr N vector들이 만들어지지 않음. 

SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


# In[8]:


cpg_type = 'opensea'


# In[6]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = []
    else: #TCGA cohort
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S
    # usage: T, N, S = get_sample_list(args.cohort)


# In[7]:


T, N, S = get_sample_list(cohort)


# In[9]:


cohort_pc1_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/1_compute-score-{cpg_type}/result/{cohort}/pc1'


# In[35]:


save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm/{cohort}'


# In[13]:


all_npzs = glob.glob(os.path.join(cohort_pc1_dir, f'*.npz')) 
bdm_pc1_files = []
for i in range(len(all_npzs)):
    if 'inv_exp' not in all_npzs[i]:
        bdm_pc1_files.append(all_npzs[i])


# In[15]:


bdm_N_pc1_files = []
for x in bdm_pc1_files:
    for n in N:
        if n in x:
            bdm_N_pc1_files.append(x)


# In[38]:


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


# In[ ]:




