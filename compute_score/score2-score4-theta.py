#!/usr/bin/env python
# coding: utf-8

# In[47]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from collections import defaultdict
#import lifelines
import matplotlib as mpl
#from statannot import add_stat_annotation
import argparse
import pickle
from sklearn.preprocessing import MinMaxScaler
import math


# ### reference
# https://www.adamsmith.haus/python/answers/how-to-calculate-the-angle-between-a-line-and-the-x-axis-in-python

# In[3]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta')


# In[4]:


# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))


# In[11]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] + ['chrX','chrY'] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]
SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'
#/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-BLCA/score2345.npz
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2



## score2 parameters
#REFERENCE = 'normal'
#METRIC = 'euclidean'


# In[6]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    parser.add_argument('-r2', '--reference_score2', help = 'Reference for score2', type = str, default = 'normal', required = False) #for score2
    parser.add_argument('-d2', '--distance_score2', help = 'Distance metric for score2', type = str, default = 'euclidean', required = False) #for score2
    parser.add_argument('-r4', '--reference_score4', help = 'Reference for score4', type = str, default = 'all', required = False) #for score2
    parser.add_argument('-d4', '--distance_score4', help = 'Distance metric for score4', type = str, default = 'euclidean', required = False) #for score2
    #s2r, s2d, s4r, s4d: normal, euclidean, all, euclidean
    return parser.parse_args()


# In[8]:


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


# In[27]:


def import_score(cohort, score, reference, distance):
    if score=='score2':
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/'+score+'_'+reference+'_'+ distance+'.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        mean_score = raw_score_df.mean(axis=1).values
        score_df = pd.DataFrame(mean_score, index = raw_score_df.index.values, columns = ['score2'])
         
    elif score=='score4':
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/'+score+'_'+reference+'_'+distance+'.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        score_df = pd.DataFrame(raw_score_df.simple_avg.values.flatten(), index = raw_score_df.index.values, columns = ['score4'])
    else: #add here if other types of scores are needed. (score1, score3, ...)
        pass
    return score_df


# In[81]:


def compute_angle(x, y):
    x = score2_norm.copy().flatten()
    y = score4_norm.copy().flatten()
    if len(x) != len(y):
        raise ValueError
    slope = np.array([float(y[i]/x[i]) for i in range(len(x))])
    radian_theta = np.array([math.atan(s) for s in slope])
    degree_theta = np.array([math.degrees(r) for r in radian_theta])
    cosine_radian = np.array([math.cos(r) for r in radian_theta])

    '''
    slope = y / x #slope of line = y / x
    radian_theta = math.atan(slope)
    degree_theta = math.degrees(radian_theta)
    cosine_radian = math.cos(radian_theta) 
    '''
    return radian_theta, degree_theta, cosine_radian



if __name__=='__main__':

    args = parse_arguments()

    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("cohort: {}".format(args.cohort))
    print("savedir: {}".format(SAVEDIR))
    
    print("1. import score 2 and score 4")
    score2_df = import_score(args.cohort, 'score2', args.reference_score2, args.distance_score2) #x: sample, y: score #index: samplename #column : ['score2']
    score4_df = import_score(args.cohort, 'score4', args.reference_score4, args.distance_score4) #x: sample, y:score #index: samplename #column : ['score4']
    
    print("2. normalize score2 and score4 to [0, 1]")
    scaler = MinMaxScaler()
    score2_norm = scaler.fit_transform(score2_df)
    scaler2 = MinMaxScaler()
    score4_norm = scaler2.fit_transform(score4_df)
    
    print("3. calculate theta (x: score2, y: score4)")
    radian_theta, degree_theta, cosine_radian = compute_angle(score2_norm, score4_norm) #input should be in order of x and y values
    
    final_df = pd.DataFrame(zip(score2_norm.flatten(), score4_norm.flatten(), radian_theta, degree_theta, cosine_radian), 
                            index = score2_df.index.values, 
                            columns = ['score2', 'score4', 'angle_radian', 'angle_degree','cos_radian'])#angle: (x, y) = (score2, score4)
    
    print("4. save result")
    result_fname = os.path.join(SAVEDIR, 'score2-score4-cosine.csv')
    final_df.to_csv(result_fname, index=True, header=True) #import 시 index_col = 0
    print('result file: {}'.format(result_fname))
    
    print("----")

