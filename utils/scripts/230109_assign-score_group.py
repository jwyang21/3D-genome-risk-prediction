#!/usr/bin/env python
# coding: utf-8

# 최종 score (refer to '/data/project/3dith/data/cohort-1-best-score-km.csv')를 사용해서, 각 cohort 내의 tumor sample들을 High group과 Low group으로 나누기.

# In[1]:


import pandas as pd
import numpy as np
import os


# In[2]:


score_files = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)


# In[3]:


num_groups = 2


# In[4]:


group_labels = ['Low', 'High']


# In[5]:


all_cpg_types = ['opensea', 'island', 'shelf_shore']


# In[6]:


score_files.index.values


# In[ ]:


#for i in range(1):
for cpg_type in all_cpg_types:
    for i in range(score_files.shape[0]):
        cohort = score_files.index.values[i]

        print("===\ncohort: {}".format(cohort))

        score_fname = score_files[f'filename_{cpg_type}'].values[i]

        score_df = pd.read_csv(score_fname, index_col = 0)
        print("score_df: ", score_df.shape)
        # tumor만 남겨야 함!!!
        tumor_mask = np.array([int(x[13:15]) <= 9 for x in score_df.index.values]) #tumor barcode ends with 01 ~ 09
        score_df = score_df.iloc[tumor_mask, :]
        print("score_df after extracting tumors only: ", score_df.shape)

        sc_df = pd.DataFrame(score_df.cos_radian.values.flatten(), index = score_df.index.values, columns = ['stem_closeness'])
        sc_df['group'] = pd.qcut(sc_df.stem_closeness, q = num_groups, labels = list(range(num_groups))) #0, 1
        group_name = [group_labels[x] for x in sc_df.group.values]
        sc_df['group_name'] = group_name

        result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/{cohort}'
        sc_df_fname = 'sc-group.csv'
        print(os.path.join(result_dir, sc_df_fname))
        sc_df.to_csv(os.path.join(result_dir, sc_df_fname))
        print("sc_df shape: ", sc_df.shape)
        #display(sc_df.head(3))


# In[ ]:




