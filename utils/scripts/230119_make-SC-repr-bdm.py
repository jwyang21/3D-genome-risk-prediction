#!/usr/bin/env python
# coding: utf-8

# In[18]:


import pandas as pd
import numpy as np
import glob
import os
import random
random.seed(42)
np.random.seed(42)


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


cohort = 'PCBC'


# In[5]:


SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


# In[6]:


cpg_type = 'opensea'


# In[7]:


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


# In[10]:


T, N, S = get_sample_list(cohort)


# In[11]:


cohort_pc1_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/1_compute-score-{cpg_type}/result/{cohort}/pc1'


# In[12]:


save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm/{cohort}'


# In[13]:


if not os.path.exists(save_dir):
    os.makedirs(save_dir)


# In[14]:


#iebdm_pc1_files = glob.glob(os.path.join(cohort_pc1_dir, f'*_inv_exp.npz'))
all_npzs = glob.glob(os.path.join(cohort_pc1_dir, f'*.npz')) 

bdm_pc1_files = []
for i in range(len(all_npzs)):
    if 'inv_exp' not in all_npzs[i]:
        bdm_pc1_files.append(all_npzs[i])

# In[20]:


random.shuffle(bdm_pc1_files)


# In[22]:


repr_size = 10


# In[25]:


num_repr = (len(bdm_pc1_files) // repr_size)


# In[29]:


save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors-bdm/{cohort}'


# In[ ]:


for chrom in CHR_LIST:
    print(f"===\n{chrom}")
    k = f'{chrom}_pc1'
    print(f"total {num_repr} repr vectors")
    for i in range(num_repr):
        print(f"---\nMaking repr vector {i}")
        data = []
        for j in range(repr_size):
            current_pointer = i * repr_size + j
            #print(current_pointer)
            current_pc1_fname = bdm_pc1_files[current_pointer]
            #print(current_pc1_fname)
            current_pc1 = np.load(current_pc1_fname)[k]
            data.append(current_pc1)
        data = np.vstack(data)
        print(data.shape)
        fname = f'repr_S_vector_{chrom}_{i}'
        np.save(os.path.join(save_dir, fname), data)
        print(os.path.join(save_dir, fname))
        del(data)
        print("---")


# In[ ]:




