#!/usr/bin/env python
# coding: utf-8

# - use pcbc reference and score3
# - SAVEDIR: /data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/pcbc
# - abs_integration_df filename: /data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/pcbc/integrate-pcbc-abs-pc1.csv
# - pcbc_reference filename: /data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/pcbc/pcbc-reference.csv

# In[6]:


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


# In[ ]:


# set current working directory
os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[ ]:


# description: compute pan-cancer score8 (ALL_COHROTS)


# In[27]:


PCBC_REFERENCE = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/scripts/result/pcbc/pcbc-reference.csv', index_col = 0).values.flatten() #순서대로 chr1~22
SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/{cohort}_score3.pickle
FIRE_COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
ALL_COHORTS = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL_COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[18]:





# In[3]:


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


# In[15]:


'''
# test in jupyter
cohort = 'TCGA-LIHC'
cohort_score3 = pd.read_pickle(os.path.join(SCORE3_DIR, cohort, cohort+'_score3.pickle'))
samples = cohort_score3.index.values

# initialize score8 dataframe
score8_df = pd.DataFrame(np.zeros(len(samples), dtype = float), index = samples, columns = ['score8'])

for s in samples:
    current_sample_score3 = cohort_score3.loc[s].values.flatten()
    score8_df.loc[s] = euclidean(current_sample_score3, PCBC_REFERENCE)
    
fname = os.path.join(SAVEDIR, 'score8.csv')
score8_df.to_csv(fname, index = True)
print(fname)
'''


# In[22]:





# In[ ]:


if __name__ == '__main__':
    args = parse_arguments() #args.cohort
    
    RESULT_DIR = os.path.join(os.getcwd(), 'result')
    SAVEDIR = os.path.join(RESULT_DIR, args.cohort)
    
    print("cohort: {}".format(args.cohort))
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIRS)
    print("RESULT_DIR: {}".format(RESULT_DIR))
    print("SAVEDIR: {}".format(SAVEDIR))
    
    cohort_score3 = pd.read_pickle(os.path.join(SCORE3_DIR, args.cohort, args.cohort+'_score3.pickle'))
    samples = cohort_score3.index.values
    
    # initialize score8 dataframe
    score8_df = pd.DataFrame(np.zeros(len(samples), dtype = float), index = samples, columns = ['score8'])
    
    for s in samples:
        current_sample_score3 = cohort_score3.loc[s].values.flatten()
        score8_df.loc[s] = euclidean(current_sample_score3, PCBC_REFERENCE)
    
    fname = os.path.join(SAVEDIR, 'score8.csv')
    score8_df.to_csv(fname, index = True) #index is sample name. 
    print(fname)


# In[ ]:




