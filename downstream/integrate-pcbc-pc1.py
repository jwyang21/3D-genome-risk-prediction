#!/usr/bin/env python
# coding: utf-8

# In[1]:


# integrate absolute value of PC1 of each pcbc sample -> average per chromosome -> pcbc SC reference
# use euclidean as a default distance. 
# pcbc 모든 샘플들 각각의 각 chromosome 별 PC1 그래프를 적분 -> 각 chromosome 별로 평균냄 -> 길이 22의 pcbc-reference vector.


# In[1]:


import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import argparse
from scipy.spatial.distance import euclidean
from scipy.integrate import simps


# In[2]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]  #all-sample-PC1 계산할때 chrX, chrY는 길이가 짧아서 mask처리하면 아예 행렬이 없어지는 경우가 있어서 상염색체만 계산했음
PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result'
FIRE_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compare-450k-fire/result' #from PC1 of inv_exp_diffmat, FIRE-intersecting bins only # should use thresholded one!!

PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#SC13-044_inv_exp.npz


# In[3]:


def get_pcbc_sample_list(directory):
    #command in main(): S= get_pcbc_sample_list(PCBC_PC1_DIR)
    files = os.listdir(directory) #PCBC_PC1_DIR includes both raw binned diffmat and inverse exponential binned diffmat
    S = []
    for f in files:
        if '_inv_exp' in f:
            S.append(f.split('_inv_exp')[0])
    return S


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort or PCBC', required = True)
    #parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)
    parser.add_argument('-f', '--flag', help = 'which type of binned diffmat you will use(raw, inv_exp)', default = 'inv_exp', required = False)

    return parser.parse_args()


# In[5]:


def import_pcbc_pc1(sample, chrom, flag):
    # import pre-computed PC1 of sample-of-interest
    
    if flag=='raw':
        fname = os.path.join(PCBC_PC1_DIR, sample+'.npz')
    elif flag=='inv':
        fname = os.path.join(PCBC_PC1_DIR, sample+'_inv_exp.npz')
    
    pc1 = np.load(fname)[chrom]

    return pc1


# In[6]:


def integrate_abs_pc1(pc1): 
    pc1_abs = np.array([abs(x) for x in pc1])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area


# In[7]:


def abs_integrate(S, flag, savedir, fname):
    #PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#SC13-044_inv_exp.npz
    abs_integration_df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), index = S, columns = CHR_LIST)
    for s in S:
        for chrom in CHR_LIST:
            pc1 = import_pcbc_pc1(s, chrom, flag)
            area = integrate_abs_pc1(pc1)
            abs_integration_df.loc[s][chrom] = area
    
    # save result
    abs_integration_df.to_csv(os.path.join(SAVEDIR, fname), index = True)
    print("abs_integration_df filename: {}".format(os.path.join(SAVEDIR, fname)))
    return abs_integration_df


# In[25]:


def get_pcbc_reference(abs_integration_df, savedir, fname):
    '''
    To-Do:
    - input: abs_integration_df #rows: all PCBC stem cell samples, columns: CHR_LIST (chr1 ~ chr22)
    - process: 각 chromosome 별 평균. 
    - output: pcbc_reference (길이 22짜리 1-d vector.)
    '''
    pcbc_reference = pd.DataFrame(abs_integration_df.mean().values, index = CHR_LIST)
    
    if len(pcbc_reference)!=len(CHR_LIST):
        raise ValueError
    
    pcbc_reference.to_csv(os.path.join(savedir, fname), index = True)
    print("pcbc_reference filename:", os.path.join(savedir, fname))
    
    return pcbc_reference    



if __name__ == '__main__':
    args = parse_arguments()
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort) #input 'pcbc' to args.cohort. 
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    #print("cohort: {}".format(args.cohort))
    #T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    #print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))
    
    S = get_pcbc_sample_list(PCBC_PC1_DIR)
    
    abs_integration_df = abs_integrate(S, args.flag, SAVEDIR, 'integrate-pcbc-abs-pc1.csv')
    pcbc_reference = get_pcbc_reference(abs_integration_df, SAVEDIR, 'pcbc-reference.csv')
    
    
    #save pcbc_reference as a csv file or pandas dataframe. 
    #score_df = compute_score(abs_integration_df, pcbc_reference)
    # abs_integration_df: row가 all pcbc stem cells, columns가 CHR_LIST인 df
    # pcbc_reference: 길이가 22인 1-d vector (one entry per one chromosome)

