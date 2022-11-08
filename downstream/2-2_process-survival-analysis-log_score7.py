#!/usr/bin/env python
# coding: utf-8

# In[1]:

# description: parse log file of survival analysis and output resulting .csv file (row: cohort, column: each clinical variable) 


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl


# In[2]:


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# ### README
# - 생존분석 결과 정리
# - score5 (score2,3,4의 조화평균)
# - log file name:
#     - 2_survival-analysis_v3.log
# - survival_analysis_v3.py로 계산한 결과

# In[3]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/')


# In[4]:


#FNAME = os.path.join(os.getcwd(), 'log', '2_score2-fire-cosine.log')
FNAME = os.path.join(os.getcwd(), 'log', '2_score7.log')#score5

CLINICAL_VARS = ['treatment_outcome_first_course', 'new_tumor_event_site', 'new_tumor_event_type', 
                 'tumor_status', 'vital_status', 'histological_grade', 'ajcc_pathologic_tumor_stage', 'race', 'gender']
P_THRESHOLD = 5e-2


# In[5]:


FNAME.split('_')[-1].split('.')[0]+ '_clinical_var_result.csv'


# ## Relationship to clinical variables

# In[6]:


#Step1. Parse log file and categorize each item.
cohort_list = []
items_list = []
f = open(FNAME, 'r')
while True:
    line = f.readline()
    if not line: break
    if 'cohort: TCGA' in line:
        current_cohort = line.split('-')[1].strip()
        cohort_list.append(current_cohort)
        #cnt += 1
        globals()[current_cohort] = [] #values #initialize value list of current cohort
        globals()[current_cohort+'_items'] = [] #items #initialzie item list of current cohort

    if 'P_val' in line:
        pvalue_str = line.split('P_val=')[1].split(' ')[0]
        if float(pvalue_str) < P_THRESHOLD:
            globals()[current_cohort].append(pvalue_str)
            globals()[current_cohort+'_items'].append(line.split(':')[0].strip())
            if line.split(':')[0].strip() not in items_list:
                items_list.append(line.split(':')[0].strip())
f.close()


# In[7]:


# Step2. 결과들을 모아서 하나의 df로 만들기.
TCGA_cohort_list = ['TCGA-'+x for x in cohort_list]

df = pd.DataFrame(index = TCGA_cohort_list, columns = items_list)

for i in range(len(cohort_list)):
    values = globals()[cohort_list[i]]
    items = globals()[cohort_list[i]+'_items']
    if len(values) != len(items):
        raise ValueError
    for j in range(len(items)):
        item_index = items_list.index(items[j]) # index of this item in df column list.
        df.iloc[i, item_index] = values[j]


# In[ ]:


'''
#step 3-1. 2-level columnn name 사용해서 결과를 깔끔하게 정리하고 싶을때.
NA_cohorts = df[df.isna().sum(axis = 1) == df.shape[1]].index.values
df2 = df.drop(NA_cohorts, axis = 0)

NA_vars = df2.columns[df2.isna().sum() == df2.shape[0]].values
df3 = df2.drop(NA_vars, axis = 1)

#df3.to_csv(os.path.join(os.getcwd(), 'result', FNAME.split('_')[-1].split('.')[0]+ '_clinical_var_result.csv'), index = True)



level1_col = ['histological_grade', 'race','vital_status', 'ajcc_pathologic_tumor_stage', 'gender', 'tumor_status'] 

level2_col = df3.columns.values

if len(level1_col) != len(level2_col):
    raise ValueError

df4 = pd.DataFrame(df3.values, index = df3.index.values, columns = [level1_col, level2_col])

# save file
fname = os.path.join(os.getcwd(), 'result', FNAME.split('_')[-1].split('.')[0]+ '_clinical_var_result.csv')
print("fname: {}".format(fname))
df4.to_csv(fname, index = True)

# check saved file
print("[FIRE, cosine]")
pd.read_csv(fname, index_col = 0)
'''


# In[8]:


#step 3-2. #그냥 저장하고 싶을때
# save file
NA_cohorts = df[df.isna().sum(axis = 1) == df.shape[1]].index.values
df2 = df.drop(NA_cohorts, axis = 0)

NA_vars = df2.columns[df2.isna().sum() == df2.shape[0]].values
df3 = df2.drop(NA_vars, axis = 1)

fname = os.path.join(os.getcwd(), 'result', FNAME.split('_')[-1].split('.')[0]+ '_clinical_var_result.csv')
print("fname: {}".format(fname))
df3.to_csv(fname, index = True)

