#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pickle
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse


# stacked bar chart: https://matplotlib.org/stable/gallery/lines_bars_and_markers/bar_stacked.html 참고

# matplotlib colors: https://matplotlib.org/stable/gallery/color/named_colors.html

# In[2]:


mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)


# ## set target variables

# In[3]:


dmr_type = 'TN'
cpg_type = 'shelf_shore'
#cpg_type='island'
#cpg_type='opensea'


# ## global variables

# In[4]:


working_dir = '/data/project/3dith/pipelines/utils/scripts'
result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/3_dmr-{cpg_type}/result/'#figure 그릴 데이터들의 위치
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
THRESHOLD = 'mean_std'
SMALL_CATEGORY_FNAME = '/data/project/3dith/data/etc/dmr-feature-small-category.npz'
#save_dir = '/data/project/3dith/figure/'#figure saving directory
#if not os.path.exists(save_dir):
#    os.makedirs(save_dir)


# In[5]:


os.chdir(working_dir)


# In[6]:


print("result_dir: {}".format(result_dir))


# In[7]:


#print("save_dir: {}".format(save_dir))


# In[8]:


# category names
if os.path.exists(SMALL_CATEGORY_FNAME):
    print("Load small category names.")
    SMALL_CATEGORY = np.load(SMALL_CATEGORY_FNAME, allow_pickle = True)
else:
    SMALL_CATEGORY = {}
    SMALL_CATEGORY['GENE'] = ['gene','transcript']
    SMALL_CATEGORY['REG'] = ['open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', 'enhancer', 'promoter',                             'promoter_flanking_region']
    SMALL_CATEGORY['EPI'] = ['18_Quies', '13_Het', '17_ReprPCWk', '16_ReprPC', '14_TssBiv', '2_TssFlnk',     '12_ZNF/Rpts', '11_EnhWk', '1_TssA', '6_TxWk', '5_Tx', '9_EnhA1', '7_EnhG1',     '4_TssFlnkD' ,'15_EnhBiv', '10_EnhA2', '3_TssFlnkU', '8_EnhG2']
    np.savez(SMALL_CATEGORY_FNAME,**SMALL_CATEGORY)


# In[9]:


# mapping TCGA cohort to EID (epigenome ID) based on tissue type.
cohort2eid = pd.read_csv('/data/project/3dith/data/etc/cohort2eid.txt', sep = '\t', header = None)
cohort2eid.columns = ['cohort', 'eid']
eid_cohorts = cohort2eid.cohort.values


# In[10]:


# legend colors
stacked_bar_color_fname = '/data/project/3dith/data/etc/stacked_bar_color.pickle'
if os.path.exists(stacked_bar_color_fname):
    print("Load legend colors")
    with open(stacked_bar_color_fname, 'rb') as f:
        stacked_bar_color = pickle.load(f)
    f.close()
else:
    stacked_bar_color = {}
    stacked_bar_color['1_TssA'] = 'red'
    stacked_bar_color['2_TssFlnk'] = 'orangered'
    stacked_bar_color['3_TssFlnkU'] = 'orangered'
    stacked_bar_color['4_TssFlnkD'] = 'orangered'
    stacked_bar_color['5_Tx'] = 'green'
    stacked_bar_color['6_TxWk'] = 'darkgreen'
    stacked_bar_color['7_EnhG1'] = 'limegreen'
    stacked_bar_color['8_EnhG2'] = 'limegreen'
    stacked_bar_color['9_EnhA1'] = 'orange'
    stacked_bar_color['10_EnhA2'] = 'orange'
    stacked_bar_color['11_EnhWk'] = 'yellow'
    stacked_bar_color['12_ZNF-Rpts'] = 'mediumseagreen'
    stacked_bar_color['13_Het'] = 'blueviolet'
    stacked_bar_color['14_TssBiv'] = 'brown'
    stacked_bar_color['15_EnhBiv'] = 'darkkhaki'
    stacked_bar_color['16_ReprPC'] = 'cornflowerblue'
    stacked_bar_color['17_ReprPCWk'] = 'royalblue'
    stacked_bar_color['18_Quies'] = 'blue'
    with open(stacked_bar_color_fname, 'wb') as f:
        pickle.dump(stacked_bar_color, f)
    f.close()


# In[11]:


chromatin_states_fname = '/data/project/3dith/data/chromatin_states.npy'
rownames = np.load(chromatin_states_fname, allow_pickle = True)
print("chromatin states: ", rownames)
for k in range(len(rownames)):
    globals()[rownames[k]] = []


# ## import data for figure

# In[12]:


data = pd.read_csv(f'{result_dir}/EPI-category-proportion-stacked-bar-chart-{cpg_type}-{dmr_type}.csv', index_col = 0)


# In[13]:


result_dir


# In[14]:


stacked_bar_color


# # draw stacked bar chart

# In[15]:


cohorts = data.index.values
cohorts = [x.split('-')[-1] for x in cohorts]
labels = cohorts

state2description = pd.read_csv('/data/project/3dith/data/etc/chromatin-states-description.csv', index_col = 0)
state2description_dict = state2description.to_dict()

for i in range(data.shape[1]):
    globals()[data.columns[i]] = data.iloc[:,i].values.flatten()

fig = plt.figure(figsize=(4,2))
ax = fig.add_subplot(111)
sum_ = [0]*data.shape[0]
for i in range(data.shape[1]):
    # Make sure that only the first letter of first word is capitalized.
    current_label = state2description_dict['description'][data.columns[i]]
    current_label_split = current_label.split(' ')
    if len(current_label_split) >= 2:
        for j in np.arange(1, len(current_label_split)):
            if current_label_split[j].lower() != 'tss' and current_label_split[j].lower() != 'polycomb':
                current_label_split[j] = current_label_split[j].lower()
    current_label = ' '.join(current_label_split)
    ax.bar(labels, data.iloc[:,i].values.flatten(), width = 0.7, label = current_label, bottom = sum_, color = stacked_bar_color[data.columns[i]])
    sum_ += data.iloc[:,i].values.flatten()
ax.set_ylim([0,1])

ax.set_ylabel('Proportion of each chromatin state', fontsize = 7)

# legend 라벨 순서 변경: https://zephyrus1111.tistory.com/24
handles, labels = ax.get_legend_handles_labels()
dict_labels_handles = dict(zip(labels, handles))
labels.reverse()
handles = [dict_labels_handles[l] for l in labels]
#print(labels)

plt.legend(handles, labels, bbox_to_anchor = [1.03, 1.05], ncols = 1, fontsize = 4.5)#, loc = 'center left')
fig.tight_layout()
fig_fname = os.path.join(result_dir, f'EPI-category-proportion-stacked-bar-chart-{cpg_type}-{dmr_type}.png')

plt.savefig(fig_fname)
print(fig_fname)


# ===================================

# In[ ]:




