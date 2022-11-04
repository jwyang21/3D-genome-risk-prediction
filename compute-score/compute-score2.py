#!/usr/bin/env python
# coding: utf-8

# In[22]:


import numpy as np
import pandas as pd
import os
import sys
from scipy.integrate import simps
from numpy import trapz
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import jensenshannon
#from sklearn.metrics.pairwise import cosine_similarity
import argparse
from numpy import dot
from numpy.linalg import norm


# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[3]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]  #all-sample-PC1 계산할때 chrX, chrY는 길이가 짧아서 mask처리하면 아예 행렬이 없어지는 경우가 있어서 상염색체만 계산했음

PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result' # FIRE와의 intersection 고려하지 않음

NORMAL_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/normal-samples-pc1/result'

FIRE_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compare-450k-fire/result' # score2 PC1으로부터 FIRE-intersecting bins only 
# Thresholded 써야 함!!! 
#TCGA-ACC  TCGA-BLCA  TCGA-LIHC  TCGA-LUAD  TCGA-LUSC  TCGA-OV  TCGA-PAAD (FIRE, 450k 모두 있는 코호트들)만 있음 (FIRE_450K_COHORT)

FIRE_450K_COHORT = ['TCGA-BLCA', 'TCGA-LUAD', 'TCGA-ACC', 'TCGA-OV', 'TCGA-LIHC', 'TCGA-LUSC', 'TCGA-PAAD']

TUMOR_NORMAL_COHORT_FNAME ='/data/project/jeewon/research/3D-ITH/pipelines/cohort-composition/result/tumor_normal_cohort.npy'
TUMOR_NORMAL_COHORT = np.load(TUMOR_NORMAL_COHORT_FNAME)


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=1e6)
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #for score2 #euclidean, jsd, cosine-similarity
    return parser.parse_args()


# In[5]:


def cosine_similarity(a, b):
    c = dot(a, b) / (norm(a) * norm(b))
    return c


# In[6]:


def get_sample_list(directory):
    file_list = os.listdir(directory)
    T = []
    N = []
    for f in file_list:
        if 'TCGA' in f and f.endswith('.npz'):
            if f[13:15] in TUMOR_BARCODE:
                T.append(f.split('.')[0])
            elif f[13:15] in NORMAL_BARCODE:
                N.append(f.split('.')[0])
            else:
                raise Exception("Wrong TCGA barcode! This sample does not belong to neither tumor nor normal.")
    S = T + N
    print("len(tumor), len(normal), len(all): {}, {} and {}, respectively.".format(len(T), len(N), len(S)))
    return T, N, S


# In[7]:


''' 
# original
def import_binned_diffmat(directory, args.score, chrom): #directory should be DATA_DIR
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T
    if args.score == 'score1':
        return diffmat, diffmat_mask, diffmat_masked
    else: #score2 or score3
        return diffmat, diffmat_mask, 1/np.exp(diffmat_masked)
'''
# revised
def import_binned_diffmat(directory, chrom): #directory should be DATA_DIR
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T
    return diffmat_masked


# In[8]:


def compute_score1(m, t): #for all scores
    cnt = (m > t).sum().sum()
    return cnt


# In[15]:


def compute_score2(v1, v2, d):
    # v1과 v2 간의 거리를 계산 (하나는 각 샘플의 PC1, 다른 하나는 reference PC1)
    # input should be flattened 1d vector!
    if d == 'euclidean':
        return euclidean(v1, v2)
    elif d == 'cosine-similarity':
        #cosine_similarity(np.array([1, 2, 3]).reshape(1, -1),np.array([1, 2, 3]).reshape(1, -1) ).flatten()[0]  #how to compute cosine similarity from two 1d vectors
        return cosine_similarity(v1, v2)
    elif d == 'jsd':
        return jensenshannon(v1, v2) 
    else:
        pass


# In[10]:


def compute_score3(pc1_450k): #test for score3 only
    pc1_abs = np.array([abs(x) for x in pc1_450k])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area


# In[ ]:


if __name__ == '__main__':
    args = parse_arguments()
    
    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    DATA_DIR = os.path.join(BINNED_DIFFMAT_DIR, args.cohort) # raw binned diffmat dir of current cohort

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("Score: {}".format(args.score))
    if args.score not in ['score1', 'score2', 'score3']:
        raise Exception("Wrong score name. Score should be in ['score1', 'score2', 'score3']")  
        
    T, N, S = get_sample_list(DATA_DIR) 
    
    # 이 코호트의 모든 샘플들의 각 chromosome별 score를 저장. 
    all_score_df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), index = S, columns = CHR_LIST)
    for s in S: #repeat for every sample
        if S.index(s) % 50 == 0: #print samplename for everh 50-th sample.
            print("s: {}".format(s))
        score_list = [] #각 chromosome의 score가 순차적으로 저장되게 # 길이 22의 1-dimensional vector.
        for chrom in CHR_LIST: #repeat for every chromosome
            #raw_m, raw_m_mask, m = import_binned_diffmat(DATA_DIR, chrom)#
            m = import_binned_diffmat(DATA_DIR, chrom)
            #m_inv_exp = 1/np.exp(m)
            if args.score == 'score1':
                score_list.append(compute_score1(m, args.threshold))
                
            if args.score == 'score2':
                if args.reference == 'fire':
                    if args.cohort not in FIRE_450K_COHORT:
                        raise Exception("Since reference is FIRE, input cohort should be included in {}".format(FIRE_450K_COHORT))
                    # import FIRE-450k intersecting bins' PC1 vectors. (of current chrom)
                    globals()['{}_pc1_fire_intersection'.format(chrom)] = np.load(os.path.join(FIRE_PC1_DIR, args.cohort, s+'_'+chrom+'_thr.npz'))['pc1_450k'].flatten()
                    globals()['{}_fire_pc1'.format(chrom)] = np.load(os.path.join(FIRE_PC1_DIR, args.cohort, s+'_'+chrom+'_thr.npz'))['pc1_fire'].flatten()
                    score_list.append(compute_score2(globals()['{}_pc1_fire_intersection'.format(chrom)], globals()['{}_fire_pc1'.format(chrom)], args.distance)) 
                            #이 sample의 이 chromosome의 score 저장                
      
                elif args.reference == 'normal':
                    if args.cohort not in TUMOR_NORMAL_COHORT:
                        raise Exception("Since referene is normal, input cohort should be included in {}".format(TUMOR_NORMAL_COHORT))
                    # import all-bin PC1 of all samples and normal reference PC1 vectors.
                    globals()['{}_pc1'.format(chrom)] = np.load(os.path.join(PC1_DIR, args.cohort, s)+'_inv_exp.npz')[chrom].flatten()  # 모든 bin의 pc1
                    globals()['{}_normal_pc1'.format(chrom)] = np.load(os.path.join(NORMAL_PC1_DIR, args.cohort, 'normal_inv_exp_pc1.npz'))[chrom].flatten()  # 모든 bin의 pc1
                    score_list.append(compute_score2(globals()['{}_pc1'.format(chrom)], globals()['{}_normal_pc1'.format(chrom)], args.distance)) 
                            #이 sample의 이 chromosome의 score 저장
                else:
                    raise Exception("reference should be either fire or normal!")
      
            elif args.score == 'score3':
                globals()['{}_pc1'.format(chrom)] = np.load(os.path.join(PC1_DIR, args.cohort, s)+'_inv_exp.npz')[chrom]
                score_list.append(compute_score3(globals()['{}_pc1'.format(chrom)])) #이 sample의 이 chromosome의 score 저장
            #m, score, pc1_450k, chrom
            else:
                pass

        all_score_df.loc[s] = score_list #이 샘플의 각 chromosome에서의 score를, 이 sample에 해당하는 row에 저장.
    
    if args.score == 'score2':
        all_score_df_fname = os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'.pickle')
    else:
        all_score_df_fname = os.path.join(SAVEDIR, args.cohort+'_'+args.score+'.pickle')
    print("All_score_df fname: {}".format(all_score_df_fname))
    all_score_df.to_pickle(os.path.join(SAVEDIR, all_score_df_fname))


