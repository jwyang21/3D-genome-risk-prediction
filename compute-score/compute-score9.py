#!/usr/bin/env python
# coding: utf-8

# - use fire reference and score9
# - /data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/fire-cohort-reference.csv
# - description: run for FIRE_COHORTS, compute score9 by calculating euclidean distance between score3 of each sample and fire-reference

# In[22]:


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from collections import defaultdict
import matplotlib as mpl
import argparse
#import pickle
from scipy.spatial.distance import euclidean
from scipy.stats import ttest_ind


# In[3]:


# set current working directory
os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[5]:


FIRE_REFERENCE_DF = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/fire-cohort-reference.csv', index_col = 0)
#PCBC_REFERENCE = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/pcbc/pcbc-reference.csv', index_col = 0).values.flatten() #순서대로 chr1~22
SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/{cohort}_score3.pickle
FIRE_COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
ALL_COHORTS = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL_COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[6]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #for score2
    #parser.add_argument('-s2r', '--score2_ref', help = 'Reference of score2', type = str, required = True)
    #parser.add_argument('-s2d', '--score2_dist', help = 'Distance metric for score2', type = str, required = True)
    #parser.add_argument('-s4r', '--score4_ref', help = 'Reference of score4', type = str, required = True)
    #parser.add_argument('-s4d', '--score4_dist', help = 'Distance metric for score4', type = str, required = True)
    
    return parser.parse_args()


# In[11]:


# test in jupyter first
#args = parse_arguments() #args.cohort
cohort = 'TCGA-LIHC'
RESULT_DIR = os.path.join(os.getcwd(), 'result')
SAVEDIR = os.path.join(RESULT_DIR, cohort)

print("cohort: {}".format(cohort))
if not os.path.exists(RESULT_DIR):
    os.makedirs(RESULT_DIR)
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIRS)
print("RESULT_DIR: {}".format(RESULT_DIR))
print("SAVEDIR: {}".format(SAVEDIR))

cohort_score3 = pd.read_pickle(os.path.join(SCORE3_DIR, cohort, cohort+'_score3.pickle'))
samples = cohort_score3.index.values
FIRE_REFERENCE = FIRE_REFERENCE_DF.loc[cohort].values.flatten()

# initialize score8 dataframe
score9_df = pd.DataFrame(np.zeros(len(samples), dtype = float), index = samples, columns = ['score9'])

for s in samples:
    current_sample_score3 = cohort_score3.loc[s].values.flatten()
    score9_df.loc[s] = euclidean(current_sample_score3, FIRE_REFERENCE)

#fname = os.path.join(SAVEDIR, 'score9.csv')
#score9_df.to_csv(fname, index = True) #index is sample name. 
#print(fname)


# In[15]:


score9_df


# In[25]:


'''
# test in jupyter(2)
# check if the hypothes that score9 is higher in tumors than normals. 
tumor_mask = np.array([int(i[13:15]) <= 9 for i in score9_df.index.values])
score9_df_tumor = score9_df.iloc[tumor_mask,:]
score9_df_normal = score9_df.iloc[~tumor_mask,:]
print('tumor score9 mean, std:', np.mean(score9_df_tumor.values.flatten()), np.std(score9_df_tumor.values.flatten()))
print('normal score9 mean, std:', np.mean(score9_df_normal.values.flatten()), np.std(score9_df_normal.values.flatten()))
'''


# In[28]:


# test in jupyter(3)
# check if score9 of tumors and normals differ significantly
#print(ttest_ind(score9_df_tumor.values.flatten(), score9_df_normal.values.flatten()))


# In[20]:


'''
# test in jupyter(4)
# plot histogram of tumor score9 and normal score9
fig = plt.figure(figsize = (6,3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.hist(score9_df_tumor.values.flatten())
ax1.set_title('tumor')
ax2.hist(score9_df_normal.values.flatten())
ax2.set_title('normal')
'''


# In[ ]:


if __name__ == '__main__':
    args = parse_arguments() #args.cohort
    
    RESULT_DIR = os.path.join(os.getcwd(), 'result')
    SAVEDIR = os.path.join(RESULT_DIR, args.cohort)
    
    print("cohort: {}".format(args.cohort))
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    print("RESULT_DIR: {}".format(RESULT_DIR))
    print("SAVEDIR: {}".format(SAVEDIR))
    
    cohort_score3 = pd.read_pickle(os.path.join(SCORE3_DIR, args.cohort, args.cohort+'_score3.pickle'))
    samples = cohort_score3.index.values
    FIRE_REFERENCE = FIRE_REFERENCE_DF.loc[args.cohort].values.flatten()
    
    # initialize score8 dataframe
    score9_df = pd.DataFrame(np.zeros(len(samples), dtype = float), index = samples, columns = ['score9'])
    
    for s in samples:
        current_sample_score3 = cohort_score3.loc[s].values.flatten()
        score9_df.loc[s] = euclidean(current_sample_score3, FIRE_REFERENCE)
    
    fname = os.path.join(SAVEDIR, 'score9.csv')
    score9_df.to_csv(fname, index = True) #index is sample name. 
    print(fname)

