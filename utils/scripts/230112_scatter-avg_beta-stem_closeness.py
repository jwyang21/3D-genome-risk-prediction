#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl


# In[2]:


mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[3]:


score_df = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)


# In[4]:


cpg_type = 'opensea'
avg_beta_fname = f'{cpg_type}_all_samples_avg_beta.csv'
fig_save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result'
fig_name = 'scatter-avg_beta-stem_closeness-ALL-cohorts.png'


# In[6]:


fig = plt.figure(figsize = (2*5, 2*3))
for i in range(score_df.shape[0]):#real
#for i in range(1):#debug

    cohort = score_df.index.values[i]
    cohort_score_fname = score_df[f'filename_{cpg_type}'].values[i]
    cohort_score_df = pd.read_csv(cohort_score_fname, index_col = 0)
    avg_beta_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/{cohort}'
    cohort_avg_beta = pd.read_csv(os.path.join(avg_beta_dir, avg_beta_fname), index_col = 0)

    assert cohort_score_df.shape[0] == cohort_avg_beta.shape[0]
    for j in range(cohort_score_df.shape[0]):
        assert cohort_score_df.index.values[j] == cohort_avg_beta.index.values[j]

    ax = fig.add_subplot(3, 5, i+1)
    ax.scatter(cohort_avg_beta.values.flatten(), cohort_score_df.cos_radian.values.flatten())
    ax.set_xlabel('avg beta')
    ax.set_ylabel('stem closeness')
    ax.set_title(f'{cohort}', fontsize = 8)
fig.tight_layout()  
full_fig_name = os.path.join(fig_save_dir, fig_name)
print(full_fig_name)
plt.savefig(full_fig_name)


# In[ ]:




