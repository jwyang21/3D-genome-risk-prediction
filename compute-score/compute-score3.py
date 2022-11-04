#!/usr/bin/env python
# coding: utf-8

# In[1]:
# To-Do: compute score3 on original binned difference matrices. (all 22 autosomes)

import numpy as np
import pandas as pd
import os
import sys
from scipy.integrate import simps
from numpy import trapz
from scipy.spatial.distance import euclidean
from sklearn.metrics.pairwise import cosine_similarity
import argparse


# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[3]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]  #all-sample-PC1 계산할때 chrX, chrY는 길이가 짧아서 mask처리하면 아예 행렬이 없어지는 경우가 있어서 상염색체만 계산했음
PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result'
FIRE_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compare-450k-fire/result' #from PC1 of inv_exp_diffmat, FIRE-intersecting bins only # should use thresholded one!!


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #for score2
    return parser.parse_args()


# In[5]:

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

# revised
def import_binned_diffmat(directory, chrom): #directory should be DATA_DIR
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T
    return diffmat_masked

def compute_score1(m, t): #for all scores
    cnt = (m > t).sum().sum()
    return cnt


def compute_score2(pc1_450k, pc1_fire, pc1_450k_fire_intersection, pc1_450k_normal, ref, dist): #for all scores
    if ref == 'FIRE': #450k_PC1 중에서 FIRE와 겹치는 bin만 씀. 
        if dist == 'euclidean':
            return euclidean(pc1_450k_fire_intersection.values.flatten(), pc1_fire.values.flatten())
        elif dist == 'cosine_sim':
            #cosine_similarity(np.array([1, 2, 3]).reshape(1, -1),np.array([1, 2, 3]).reshape(1, -1) ).flatten()[0]  #how to compute cosine similarity from two 1d vectors
            return cosine_similarity(pc1_450k_fire_intersection.values.flatten().reshape(1, -1), pc1_fire.values.flatten().reshape(1, -1)).flatten()[0]
        else:
            pass

    elif ref == 'normal':  #normal을 reference로 쓰기 때문에 tumor sample에 대해서만 계산이 가능하다. #FIRE와의 교집합 고려하지 않고 PC1 그냥 그대로 씀.  
        if dist == 'euclidean':
            return euclidean(pc1_450k.values.flatten(), pc1_450k_normal.values.flatten())
        elif dist == 'cosine_sim':
            return cosine_similarity(pc1_450k.values.flatten().reshape(1, -1), pc1_450k_normal.values.flatten().reshape(1, -1)).flatten()[0]
        else:
            pass
    else:
        pass


def compute_score3(pc1_450k): #test for score3 only
    pc1_abs = np.array([abs(x) for x in pc1_450k])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area


if __name__ == '__main__':
    args = parse_arguments()
    
    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    DATA_DIR = os.path.join(BINNED_DIFFMAT_DIR, args.cohort) # binned diffmat dir of current cohort

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
        
    print("Score: {}".format(args.score))
    if args.score not in ['score1', 'score2', 'score3']:
        raise Exception("Wrong score name. Score should be in ['score1', 'score2', 'score3']")  
        
    T, N, S = get_sample_list(DATA_DIR) #
    
    all_score_df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), index = S, columns = CHR_LIST)
    for s in S: #repeat for every sample
        if S.index(s) % 50 == 0: #print samplename for everh 50-th sample.
            print("s: {}".format(s))
        score_list = [] #각 chromosome의 score가 순차적으로 저장되게
        for chrom in CHR_LIST: #repeat for every chromosome
            
            m = import_binned_diffmat(DATA_DIR, chrom)
            
            if args.score == 'score1':
                score_list.append(compute_score1(m, args.threshold))
                
            if args.score == 'score2':
                globals()['{}_pc1_fire_intersection'.format(chrom)] = np.load(os.path.join(FIRE_PC1_DIR, args.cohort, s+'_'+chrom+'_thr.npz'))['pc1_450k']
                globals()['{}_fire_pc1'.format(chrom)] = np.load(os.path.join(FIRE_PC1_DIR, args.cohort, s+'_'+chrom+'_thr.npz'))['pc1_fire']
                globals()['{}_pc1'.format(chrom)] = np.load(os.path.join(PC1_DIR, args.cohort, s)+'_inv_exp.npz')[chrom]  
                score_list.append(compute_score2(globals()['{}_pc1'.format(chrom)], globals()['{}_fire_pc1'.format(chrom)], 
                            globals()['{}_pc1_fire_intersection'.format(chrom)], pc1_450k_normal, args.reference, args.distance)) #이 sample의 이 chromosome의 score 저장                
            
            elif args.score == 'score3':
                globals()['{}_pc1'.format(chrom)] = np.load(os.path.join(PC1_DIR, args.cohort, s)+'_inv_exp.npz')[chrom]
                score_list.append(compute_score3(globals()['{}_pc1'.format(chrom)])) #이 sample의 이 chromosome의 score 저장
            
            else:
                pass


        all_score_df.loc[s] = score_list #이 샘플의 각 chromosome에서의 score를 저장. 
    
    # Save df as a pickle file.
    if args.score == 'score2':
        all_score_df_fname = os.path.join(SAVEDIR, args.cohort+'_'+args.score+'_'+args.reference+'_'+args.distance+'.pickle')
    else:
        all_score_df_fname = os.path.join(SAVEDIR, args.cohort+'_'+args.score+'.pickle')
    print("All_score_df fname: {}".format(all_score_df_fname))
    all_score_df.to_pickle(os.path.join(SAVEDIR, all_score_df_fname))
