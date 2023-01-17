#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr


# In[ ]:


import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
#plt.rc('xtick', labelsize=7)
#plt.rc('ytick', labelsize=7)


# In[ ]:


scores = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)


# ## confirm that all normal and stem distance values are nonnegative

# In[ ]:


for cpg_type in ['island', 'opensea', 'shelf_shore']:
    print(f"===\n{cpg_type}")
    for i in range(scores.shape[0]):
        print("---")
        print(scores.index.values[i])
        fname = scores[f'filename_{cpg_type}'].values[i]
        f_ = pd.read_csv(fname, index_col = 0)
        #display(f_)
        print('negative value in normal distance: ', (f_.normal_distance.values < 0).sum())
        print('negative value in stem distance: ', (f_.normal_distance.values < 0).sum())


# In[ ]:


cpg_types = ['opensea', 'island', 'shelf_shore']
lowest_subdir = 'normal-stem-distance-scatter' # name of the lowest subdirectory where the results will be saved in. 


# In[ ]:


result_fname_s = 'scatter-normal-stem-distances'
result_fname_h = 'histogram-pcc-normal-stem-distances'


# ## scatter plot of normal and stem distancaes

# In[ ]:


for cpg_type in cpg_types:
    print(f"===\n{cpg_type}")
    scatter_values = {}
    result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result'
    save_dir = f'{result_dir}/{lowest_subdir}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    #print(save_dir)
    
    fig = plt.figure(figsize = (2*5, 2*3))
    distance_corr_df = pd.DataFrame(np.zeros((len(scores.index.values), 1), dtype = float), index = scores.index.values, columns = ['corr_N_S_distance'])
    for i in range(scores.shape[0]):
        ax = fig.add_subplot(3, 5, i+1)
        cohort = scores.index.values[i]
        #print(f'---\n{cohort}')
        #print(cohort)
        score = pd.read_csv(scores[f'filename_{cpg_type}'].values[i], index_col = 0)
        ax.scatter(score.normal_distance.values, score.stem_distance.values)
        ax.set_xlabel('normal distance')
        ax.set_ylabel('stem distance')
        ax.set_title(cohort)
        distance_corr_df.loc[cohort] = pearsonr(score.normal_distance.values, score.stem_distance.values)[0]
        scatter_values[cohort] = score[['normal_distance', 'stem_distance']]
    fig.tight_layout()
    plt.savefig(os.path.join(save_dir, result_fname_s+'.png'))
    plt.clf()
    np.savez(os.path.join(save_dir, result_fname_s), **scatter_values)
    del(scatter_values)
    print(os.path.join(save_dir, result_fname_s+'.png'))
    print(os.path.join(save_dir, result_fname_s+'.npz'))
    
    fig = plt.figure(figsize = (2.5, 1.5))
    ax = fig.add_subplot(111)
    ax.hist(distance_corr_df.values.flatten())
    ax.set_xlabel('pcc (normal distance, stem distance)')
    ax.set_ylabel('frequency')
    ax.set_title(f'Histogram of correlation values ({cpg_type})', fontsize = 8)
    distance_corr_df.to_csv(os.path.join(save_dir, result_fname_h+'.csv'))
    fig.tight_layout()
    plt.savefig(os.path.join(save_dir, result_fname_h+'.png'))
    print(os.path.join(save_dir, result_fname_h+'.csv'))
    print(os.path.join(save_dir, result_fname_h+'.png'))
    display(distance_corr_df)


# ## histogram of pcc (normal distance, stem distance)

# In[ ]:


fig = plt.figure(figsize = (2.5, 1.5))
plt.rc('xtick', labelsize=5)
plt.rc('ytick', labelsize=5)
ax = fig.add_subplot(111)
ax.hist(distance_corr_df.values.flatten())
ax.set_xlabel('PCC (normal distance, stem distance)', fontsize = 5)
ax.set_ylabel('Frequency', fontsize = 5)
ax.set_title('Histogram of correlation between normal and stem distances', fontsize = 5)

print(f'/data/project/3dith/figure/{cpg_type}_N_S_distance_corr_histogram.png')


fig.tight_layout()
plt.savefig(f'/data/project/3dith/figure/{cpg_type}_N_S_distance_corr_histogram.png')


# In[ ]:




