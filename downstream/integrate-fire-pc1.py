#!/usr/bin/env python
# coding: utf-8

# In[3]:
# description: integrate FIRE PC1 values to get fire reference. (use for score9)

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
from sklearn.decomposition import PCA


# In[4]:


# FIRE의 PC1들 중 score2의 bin들과 겹치는 bin들만 찾아서 그래프 만든 후 적분 -> 각 chromosome 별로 평균 내면 길이 22의 normal-reference vector.


# In[40]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]  #all-sample-PC1 계산할때 chrX, chrY는 길이가 짧아서 mask처리하면 아예 행렬이 없어지는 경우가 있어서 상염색체만 계산했음
PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result'
FIRE_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compare-450k-fire/result' #from PC1 of inv_exp_diffmat, FIRE-intersecting bins only # should use thresholded one!!
CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
### int(CHR_LENGTH[CHR_LENGTH['Chromosome'] == 'chr1']['Total_length'].values[0])
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#SC13-044_inv_exp.npz
FIRE_COHORT = ['TCGA-BLCA','TCGA-LUAD','TCGA-ACC','TCGA-OV','TCGA-LIHC','TCGA-LUSC','TCGA-PAAD']
TCGA2FIRE = pd.read_csv('/data/project/jeewon/research/3D-ITH/data/tcga2fire.csv')#columns: TCGA, FIRE
### str(TCGA2FIRE[TCGA2FIRE['TCGA'] == 'TCGA-ACC']['FIRE'].values[0])
FIRE_PC1_FNAME = '/data/project/3dith/data/fire_pc1.csv'
FIRE_PC1 = pd.read_csv(FIRE_PC1_FNAME)
FIRE_BINSIZE = int(1e6)
# bin name format: chr1:1-1000000
FIRE_PC1['end'] = FIRE_PC1['start']+FIRE_BINSIZE-1


# In[ ]:





# In[49]:


FIRE_PC1[FIRE_PC1['LI'].notnull()]


# In[51]:


FIRE_PC1


# In[52]:


FIRE_PC1_LI =FIRE_PC1[FIRE_PC1.LI.notnull()]


# In[53]:


FIRE_PC1_LI[FIRE_PC1_LI['chr']==1]


# In[4]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    cohort_binned_diffmat_dir = os.path.join(BINNED_DIFFMAT_DIR, cohort)
    T = []
    N = []
    S = [] #all samples
    for l in os.listdir(cohort_binned_diffmat_dir):
        if l.startswith('TCGA') and l.endswith('.npz'):
            if l[13:15] in TUMOR_BARCODES:
                T.append(l.split('.')[0].strip())
            elif l[13:15] in NORMAL_BARCODES:
                N.append(l.split('.')[0].strip())
            else:
                pass
        else:
            pass
    S = T + N
    print("{}: tumor {}, normal {}, total {}".format(cohort, len(T), len(N), len(S)))

    return T, N, S


# In[6]:


def import_pc1(cohort, sample, chrom, flag):
    # import pre-computed PC1 of sample-of-interest
    if flag=='raw':
        fname = os.path.join(PC1_DIR, cohort, sample+'.npz')
    elif flag=='inv':
        fname = os.path.join(PC1_DIR, cohort, sample+'_inv_exp.npz')
    pc1 = np.load(fname)[chrom]

    return pc1


# In[65]:


def integrate_abs_pc1(pc1_450k): 
    pc1_abs = np.array([abs(x) for x in pc1_450k])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area


# In[ ]:


if __name__ == '__main__':
    
    RESULT_DIR = os.path.join(os.getcwd(), 'result')
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    print("RESULT_DIR: {}".format(RESULT_DIR))
    
    #initialize fire_reference df
    fire_reference = pd.DataFrame(np.zeros((len(FIRE_COHORT), len(CHR_LIST)), dtype = float), index = FIRE_COHORT, columns = CHR_LIST)
    
    for cohort in FIRE_COHORT: #이 cohrot는 TCGA-~~~ 형식

        cohort_fire = str(TCGA2FIRE[TCGA2FIRE['TCGA'] == cohort]['FIRE'].values[0])

        cohort_fire_data = FIRE_PC1[FIRE_PC1[cohort_fire].notnull()]

        for chrom in np.arange(1, 23):#chr1 ~ chr22
            current_fire_pc1 = cohort_fire_data[cohort_fire_data['chr']==chrom][cohort_fire].values

            fire_reference.loc[cohort][CHR_LIST[chrom-1]] = integrate_abs_pc1(current_fire_pc1)

    fname = os.path.join(RESULT_DIR, 'fire-cohort-reference.csv')

    fire_reference.to_csv(fname, index = True)
    print(fname)
    print("----\n")
    #이제 얘를 각 샘플의 score3과 비교해서 각 chromosome 별 값 간의 차이들 간 euclidean distance로 score9 계산

