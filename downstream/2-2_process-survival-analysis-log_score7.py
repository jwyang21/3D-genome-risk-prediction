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


# In[9]:


#to check
df3


# -----------------
# below this line are codes for score2 only!!
## Tumor_normal ttest results

# In[18]:


ref = ['normal', 'fire']
metric = ['euclidean', 'cosine-similarity']
ALL_COHORT = 'TCGA-CESC TCGA-ESCA TCGA-KIRC TCGA-LIHC TCGA-OV TCGA-READ TCGA-TGCT TCGA-UCS TCGA-ACC TCGA-CHOL TCGA-GBM TCGA-KIRP TCGA-LUAD TCGA-PAAD TCGA-SARC TCGA-THCA TCGA-UVM TCGA-BLCA TCGA-COAD TCGA-HNSC TCGA-LAML TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-THYM TCGA-BRCA TCGA-DLBC TCGA-KICH TCGA-LGG TCGA-MESO TCGA-PRAD TCGA-STAD TCGA-UCEC'.split(' ')
FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
TUMOR_NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[16]:


'''
score2_normal_cosine-similarity_survival_analysis_pvalues.npz
score2_normal_cosine-similarity_tumor_normal_ttest.npz

score2_fire_cosine-similarity_survival_analysis_pvalues.npz
score2_fire_cosine-similarity_tumor_normal_ttest.npz

'''
tumor_normal_ttest_score5.npz
survival_analysis_pvalues_score5.npz
#clinical_score5_merged_TumorOnly.pickle
#subtype_dictionary_score5.pickle 


# In[22]:


# NORMAL cohort들의 각 metric에서의 tumor v.s. normal ttest 결과
tumor_normal_ttest_items = ['tumor_mean', 'tumor_std', 'normal_mean', 'normal_std', 'tumor_normal_ttest_p']
tumor_normal_ttest_df = pd.DataFrame(np.zeros((len(TUMOR_NORMAL_COHORT), len(tumor_normal_ttest_items)), dtype = float), index = TUMOR_NORMAL_COHORT, columns = tumor_normal_ttest_items)
for cohort in TUMOR_NORMAL_COHORT:
    ttest_npz_fname = 'tumor_normal_ttest_score5.npz'
    ttest_npz_fullname = os.path.join(os.getcwd(), 'result', cohort, ttest_npz_fname)
    ttest_values = np.load(ttest_npz_fullname)['values']
    tumor_normal_ttest_df.loc[cohort] = ttest_values
df_fname = os.path.join(os.getcwd(), 'result', 'tumor_normal_ttest_score5.csv')
print("df_fname: {}".format(df_fname))
tumor_normal_ttest_df.to_csv(df_fname, index = True)        


# In[24]:


tumor_normal_ttest_df.head(3)


# In[41]:


# Plot results
fig = plt.figure(figsize = (4*2, 4))
ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(tumor_normal_ttest_df.tumor_mean.values.flatten(), label = 'Tumor', linewidth = 3)
ax1.plot(tumor_normal_ttest_df.normal_mean.values.flatten(), label = 'Normal', linewidth = 3)
#ax1.legend()
ax1.set_xlabel('TCGA cohorts')
#ax1.set_ylabel('Values')
ax1.set_ylabel('Mean')
#ax1.set_xticks(np.arange(normal_euclidean2.shape[0]), normal_euclidean2.index.values)
ax1.set_title('score5')

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(tumor_normal_ttest_df.tumor_mean.values.flatten(), label = 'Tumor', linewidth = 3)
ax2.plot(tumor_normal_ttest_df.normal_mean.values.flatten(), label = 'Normal', linewidth = 3)
#ax2.legend()
ax2.set_xlabel('TCGA cohorts')
#ax2.set_ylabel('Values')
ax2.set_ylabel('Std')
#ax2.set_xticks(np.arange(normal_cosine2.shape[0]), normal_cosine2.index.values)
ax2.set_title('score5')

fig.suptitle('Scores of tumor and normal samples (Score5)', fontsize = 15)
ax2.legend(frameon = False, loc = 'center left', bbox_to_anchor = (1.05, 0.5))

fig.tight_layout()
plt.show()
plt.savefig('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/score5_tumor_normal.png')
print('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/score5_tumor_normal.png')


# In[42]:


plt.clf()


# In[ ]:





# ## Survival analysis

# In[44]:


for cohort in TUMOR_NORMAL_COHORT:
    print("{}: ".format(cohort), end = '')
    survival_npz_fname = 'survival_analysis_pvalues_score5.npz'
    survival_npz_fullname = os.path.join(os.getcwd(), 'result', cohort, survival_npz_fname)
    sig_t = np.load(survival_npz_fullname)['significant_target']
    sig_t_pvals = np.load(survival_npz_fullname)['significant_target_pvals']
    if len(sig_t)> 0:
        for i in range(len(sig_t)):
            print(sig_t[i]+' (p = ', end = '')
            print(str(sig_t_pvals[i])+')', end = '')
            if i < len(sig_t)-1:
                print(', ', end = '')


    print("\n")


# In[ ]:




