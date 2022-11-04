#!/usr/bin/env python
# coding: utf-8

# score6: score2 * score3 / score4
# score2: normal과의 거리, score3: perturbation의 정도, score4: SC와의 거리
# normal과 거리가 멀 수록, perturbation이 심할 수록, SC에 가까울 수록 cancer cell-like라고 생각.

# In[1]:


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from collections import defaultdict
import lifelines
import matplotlib as mpl
from statannot import add_stat_annotation
import argparse
import pickle


# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[3]:


# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))


# In[4]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] + ['chrX','chrY'] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/' 

CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2

SCORE2_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/score2_{reference}_{metric}.pickle    # normal, euclidean
SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/{cohort}_score3.pickle
SCORE4_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/score4_{ref_type}_{distance}.pickle   # all, euclidean


# In[5]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #for score2
    parser.add_argument('-s2r', '--score2_ref', help = 'Reference of score2', type = str, required = True)
    parser.add_argument('-s2d', '--score2_dist', help = 'Distance metric for score2', type = str, required = True)
    parser.add_argument('-s4r', '--score4_ref', help = 'Reference of score4', type = str, required = True)
    parser.add_argument('-s4d', '--score4_dist', help = 'Distance metric for score4', type = str, required = True)
    
    return parser.parse_args()


# In[6]:


def get_sample_list(directory):
    print("Get_sample_list")
    t = []
    n = []
    for f in os.listdir(directory):
        if 'TCGA' in f and '.npz' in f:
            if f[13:15] in TUMOR_BARCODE:
                t.append(f[:15])
            elif f[13:15] in NORMAL_BARCODE:
                n.append(f[:15])
            else:
                print("Neither tumor nor normal")
    T = list(set(t))
    N = list(set(n))
    S = T + N
    return T, N, S


# In[15]:


def import_score(S, score2_ref, score2_metric, score4_ref, score4_metric, cohort):
    s2_df = pd.read_pickle(os.path.join(SCORE2_DIR, cohort, 'score2_'+score2_ref+'_'+score2_metric+'.pickle')).loc[S]
    s3_df = pd.read_pickle(os.path.join(SCORE3_DIR, cohort, cohort+'_score3.pickle')).loc[S]
    s4_df = pd.read_pickle(os.path.join(SCORE4_DIR, cohort, 'score4_'+score4_ref+'_'+score4_metric+'.pickle')).loc[S]
    df = pd.DataFrame(np.zeros((len(S), 3),dtype=float), index = S, columns = ['score2', 'score3', 'score4'])
    for i, s in enumerate([s2_df, s3_df, s4_df]):
        if 'simple_avg' in s:
            df['score{}'.format(i+2)] = s.simple_avg.values.flatten()
        else:
            avg = s.mean(axis = 1).values.flatten()
            df['score{}'.format(i+2)] = avg
    return df


# In[8]:


def compute_score(df):
    #(1) 각 score를 1D numpy array로 표현.     
    s2 = df.score2.values.flatten()
    s3 = df.score3.values.flatten()
    s4 = df.score4.values.flatten()
    
    # (2) score6 = score2 * score3 / score4
    s6 = s2 * s3 / s4
    
    # (3) df에 score6 column을 추가
    df['score6'] = s6
    
    return df 


# In[11]:

if __name__ == '__main__':
    args = parse_arguments()
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("cohort: {}".format(args.cohort))
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))

    #score_df = import_score(S, 'normal', 'euclidean', 'all', 'euclidean')#row가 sample, column이 score2, score3, score4
    score_df = import_score(S, args.score2_ref, args.score2_dist, args.score4_ref, args.score4_dist, args.cohort)
    
    score_df2 = compute_score(score_df) #score_df에 score5 column을 추가하여 반환
    npz_fname = os.path.join(SAVEDIR, 'score2346')
    np.savez(npz_fname, rownames = score_df2.index.values, colnames = score_df2.columns.values, score2 = score_df2.score2.values.flatten(),
             score3 = score_df2.score3.values.flatten(), score4 = score_df2.score4.values.flatten(), score6 = score_df2.score6.values.flatten())
    print(npz_fname+'.npz')
    print('----\n')

