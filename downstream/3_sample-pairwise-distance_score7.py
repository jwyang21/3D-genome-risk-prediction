#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
import argparse
#from scipy.spatial import distance
from sklearn.manifold import MDS
#import matplotlib as mpl
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import jensenshannon
from numpy import dot
from numpy.linalg import norm


# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses')


# In[3]:


'''
# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))
'''


# In[4]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
#SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/' # score3
SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta/result/'#/TCGA-BLCA/score2-score4-cosine.csv

#PSEUDOCOUNT = 1e-6
#score3 pickle file: /data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/{cohort}/{cohort}_score3.pickle

#FIRE_DIR = '/data/project/3dith/data/'
#FIRE_AB_FNAME = 'fire_a_b_labels.csv'
#FIRE_PC1_FNAME = 'fire_pc1.csv'
#COHORT_TCGA_FIRE_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA_FIRE_intersection.csv'
#COHORT_TCGA_FIRE = pd.read_csv(COHORT_TCGA_FIRE_FNAME)
#CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
#P_THRESHOLD = 5e-2
#SUBTYPE_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA_subtypes.csv'
#SUBTYPE = pd.read_csv(SUBTYPE_FNAME)


# In[8]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True) #어떤 score로 구한 PC1으로 sample 간 거리를 계산하는지.
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    parser.add_argument('-d', '--distance', help = 'Distance metric to be used to compute how similar/different two samples are.', type = str, default = 'euclidean', required = False) 
    ###which distance metric you want to use
    return parser.parse_args()


# In[7]:


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


# In[9]:


def cosine_similarity(a, b):
    c = dot(a, b) / (norm(a) * norm(b))
    return c


# In[85]:


'''
# for score1, score2
# 지금은 함수 안에서 pc1을 계산하는데, pre-computed pc1을 import하도록 코드 수정 필요
def sample_pairwise_distance(directory, S, chrom):
    s_s_df = pd.DataFrame(np.zeros((len(S), len(S)), dtype = float), index = S, columns = S)
    for s in S:
        if S.index(s) % 100 == 0: 
            print("s: {}".format(s))
        s_npz = os.path.join('/data/project/3dith/pipelines/binned-difference-matrix-v2/result/', cohort, s+'.npz')
        s_diffmat = np.load(s_npz)['chr'+str(chrom)]
        s_diffmat_mask = np.load(s_npz)['chr'+str(chrom)+'_mask']
        s_diffmat_masked = s_diffmat[~s_diffmat_mask].T[~s_diffmat_mask].T #score 1
        s_diffmat_masked_inverse = 1/np.exp(s_diffmat_masked) #score 2
        pca = PCA(n_components = 2)
        s_pc1 = pca.fit_transform(s_diffmat_masked_inverse)[:,0]
        for s2 in S: #all-sample-pairwise
            if S.index(s) % 100 == 0 and S.index(s2) % 100 == 0:
                print("s2: {}".format(s2))
            s2_npz = os.path.join('/data/project/3dith/pipelines/binned-difference-matrix-v2/result/', cohort, s2 +'.npz')
            s2_diffmat = np.load(s2_npz)['chr'+str(chrom)]
            s2_diffmat_mask = np.load(s2_npz)['chr'+str(chrom)+'_mask']
            s2_diffmat_masked = s2_diffmat[~s2_diffmat_mask].T[~s2_diffmat_mask].T #score 1
            s2_diffmat_masked_inverse = 1/np.exp(s2_diffmat_masked) #score 2
            pca2 = PCA(n_components = 2)
            s2_pc1 = pca2.fit_transform(s2_diffmat_masked_inverse)[:,0]
            s_s_df.loc[s, s2] = distance.euclidean(s_pc1, s2_pc1)

    return s_s_df
'''


# In[28]:


'''
def sample_pairwise_distance_score3(cohort, distance):
    score_df = pd.read_pickle(os.path.join(SCORE_DIR, cohort, cohort+'_score3.pickle'))
    s_s_df = pd.DataFrame(np.zeros((score_df.shape[0], score_df.shape[0]), dtype = float), index = score_df.index.values, columns = score_df.index.values) # sample by sample
    for i in range(score_df.shape[0]): #iterate for all samples
        current_v = score_df.iloc[i,:]
        for j in range(score_df.shape[0]): #iterate for al samples (probe sample)
            current_probe_v = score_df.iloc[j,:]
            if distance == 'euclidean':
                s_s_df.iloc[i, j] = euclidean(current_v, current_probe_v)
            elif distance == 'cosine-similarity':
                s_s_df.iloc[i, j] = cosine_similarity(current_v, current_probe_v)
            else: #'jsd'
                s_s_df.iloc[i, j] = jensenshannon(current_v, current_probe_v)    
    return s_s_df
'''


# In[17]:


def sample_pairwise_distance_score7(cohort, distance):
    score_df = pd.read_csv(os.path.join(SCORE_DIR, cohort, 'score2-score4-cosine.csv'), index_col = 0)
    s_s_df = pd.DataFrame(np.zeros((score_df.shape[0], score_df.shape[0]), dtype = float), index = score_df.index.values, columns = score_df.index.values) # sample by sample
    for i in range(score_df.shape[0]): #iterate for all samples
        current_v = score_df.iloc[i,:]
        for j in range(score_df.shape[0]): #iterate for al samples (probe sample)
            current_probe_v = score_df.iloc[j,:]
            if distance == 'euclidean':
                s_s_df.iloc[i, j] = euclidean(current_v, current_probe_v)
            elif distance == 'cosine-similarity':
                s_s_df.iloc[i, j] = cosine_similarity(current_v, current_probe_v)
            else: #'jsd'
                s_s_df.iloc[i, j] = jensenshannon(current_v, current_probe_v)    
    return s_s_df




if __name__=='__main__':
    
    #for cohort in COHORTS:
    args = parse_arguments()
    #cohort = args.cohort
    print("cohort: {}".format(args.cohort))
    
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    
    if args.score not in ['score1', 'score2', 'score3', 'score4','score5','score6','score7']:
        raise ValueError
        
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    print("SAVEDIR: {}".format(SAVEDIR))
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    #d = os.path.join(os.getcwd(), cohort) #.npz file들이 있는 디렉토리
    if args.score == 'score1' or args.score == 'score2':
        for c in CHR_LIST:
            print('chr'+c)
            s_s_df = sample_pairwise_distance(d, S, c)
            #save_result(s_s_df, d, c)
    elif args.score == 'score3': #score3
        s_s_df = sample_pairwise_distance_score3(args.cohort, args.distance)
    elif args.score == 'score7':
        s_s_df = sample_pairwise_distance_score7(args.cohort, args.distance)
    
    
    fname = args.score + '_sample_sample_'+args.distance+'.csv'
    print("all-sample-pairwise distance file: {}".format(os.path.join(SAVEDIR, fname)))
    s_s_df.to_csv(os.path.join(SAVEDIR, fname), index = True)
    print("----")






