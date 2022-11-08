#!/usr/bin/env python
# coding: utf-8

# In[14]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from collections import defaultdict
import matplotlib as mpl
import argparse


# In[5]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta')


# In[6]:


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
SCORE7_FNAME = 'score2-score4-cosine.csv'


## score2 parameters
#REFERENCE = 'normal'
#METRIC = 'euclidean'


# In[8]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r2', '--reference_score2', help = 'Reference for score2', type = str, default = 'normal', required = False) #for score2
    #parser.add_argument('-d2', '--distance_score2', help = 'Distance metric for score2', type = str, default = 'euclidean', required = False) #for score2
    #parser.add_argument('-r4', '--reference_score4', help = 'Reference for score4', type = str, default = 'all', required = False) #for score2
    #arser.add_argument('-d4', '--distance_score4', help = 'Distance metric for score4', type = str, default = 'euclidean', required = False) #for score2
    #s2r, s2d, s4r, s4d: normal, euclidean, all, euclidean
    return parser.parse_args()


# In[12]:


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


# In[13]:


def import_score7(cohort):
    score_fname = os.path.join(os.getcwd(), 'result', cohort, SCORE7_FNAME)
    raw_score_df = pd.read_csv(score_fname, index_col = 0)
    score = raw_score_df['cos_radian'].values.flatten()
    score_df = pd.DataFrame(score, index = raw_score_df.index.values, columns = ['score7'])
    return score_df


# In[ ]:


def plot_tumor_normal(score_df, fig_size, cohort, savedir, fig_name):
    # plot scores of tumors and normals with label.
    group  = ['Tumor' if int(x[13:15])<=9 else 'Normal' for x in score_df.index.values]
    score_df2 = score_df.copy()
    score_df2['group'] = group
    score_df2_tumor = score_df2[score_df2['group']=='Tumor'].copy()
    score_df2_normal = score_df2[score_df2['group']=='Normal'].copy()
    fig = plt.figure(figsize = (2*fig_size, fig_size))
    ax1 = fig.add_subplot(121)
    ax1.hist(score_df2_tumor.score7.values)
    ax1.set_title('score7 histogram (Tumor)')
    ax2 = fig.add_subplot(122)
    ax2.hist(score_df2_normal.score7.values)
    ax2.set_title('score7 histogram (Normal)')
    fig.suptitle(cohort)
    fname = os.path.join(savedir, fig_name)
    plt.savefig(fname)
    print(fname)
    return score_df2, score_df2_tumor, score_df2_normal


# In[ ]:


def tumor_normal_ttest(score_df2_tumor, score_df2_normal, savedir, cohort, file_name):
    pval = ttest_ind(score_df2_tumor.score7.values.flatten(), score_df2_normal.score7.values.flatten())[1]
    stat = ttest_ind(score_df2_tumor.score7.values.flatten(), score_df2_normal.score7.values.flatten())[0]
    fname = os.path.join(savedir, file_name)
    if os.path.exists(fname):
        with open(fname, "a") as file:#이어쓰기
            file.write(cohort)
            file.write("\n")
            file.write('pvalue\t')
            file.write(str(pval))
            file.write("\n")
            file.write('stat\t')
            file.write(str(stat))
            file.write("\n")
            file.write("----\n")
            file.close()
    else: #이 text file이 존재하지 않으면 새로 만듦. 
        with open(fname, "w") as file:#이어쓰기
            file.write(cohort)
            file.write("\n")
            file.write('pvalue\t')
            file.write(str(pval))
            file.write("\n")
            file.write('stat\t')
            file.write(str(stat))
            file.write("\n")
            file.write("----\n")
            file.close()  

    


# In[ ]:


if __name__=='__main__':

    args = parse_arguments()
    
    RESULTDIR = os.path.join(os.getcwd(), 'result')
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("cohort: {}".format(args.cohort))
    print("savedir: {}".format(SAVEDIR))
    
    print("1. import score7")
    score_df = import_score7(args.cohort)
    
    print("2. plot score7 of tumors and normals with label")
    score_df2, score_df2_tumor, score_df2_normal = plot_tumor_normal(score_df, 3, args.cohort, SAVEDIR, 'tumor-normal-score7-histogram.png') #add 'group' column to score_df
    
    print("3. t-test between tumors and normals' score7")
    tumor_normal_ttest(score_df2_tumor, score_df2_normal, RESULTDIR, args.cohort, 'tumor-normal-score7-ttest.txt') #모든 cohort들의 ttest 결과를 하나의 txt file로 모아서 저장.   

