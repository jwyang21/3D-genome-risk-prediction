#!/usr/bin/env python
# coding: utf-8

# In[37]:


import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import random
random.seed(2022)
np.random.seed(2022)
import math
import os
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import matplotlib as mpl


# ### reference
# https://www.adamsmith.haus/python/answers/how-to-calculate-the-angle-between-a-line-and-the-x-axis-in-python

# In[4]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score')


# In[38]:


# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))


# In[24]:


CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
#DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
#BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
# TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
#NORMAL_BARCODES = ['10', '11', '12', '13','14', '15', '16', '17', '18', '19']
## TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
#CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
#CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
#TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{cohort}/{sample}.npz or /{cohort}/{sample}_inv_exp.npz  #'chr1'
#PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/PCBC/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
#METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
#PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
#FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
#NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
#ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
#TCGA_SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/score3_simple_avg.pickle
#PCBC_SCORE3_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv'
SAVEDIR = os.path.join(os.getcwd(), 'result')
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
#ALL_COHORT_W_PCBC = ALL_COHORT.copy()
#ALL_COHORT_W_PCBC.append('PCBC')
print("SAVEDIR: {}".format(SAVEDIR))
#BINNED_DIFFMAT_BINS = '/data/project/jeewon/research/3D-ITH/pipelines/etc/binned-diffmat-bins' #pcbc_bins.npz #{TCGA_cohort}_diffmat_bins.npz
P_THRESHOLD = 5e-2
NORMAL_STEM_DISTANCE_MINMAX = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/normal-stem-distance-min-max.csv', index_col = 0)


# In[8]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort or PCBC', type = str, required = True)
    parser.add_argument('-s_type', '--score_type', help = 'which one you will use among PC1 fluctuation or averaged PC1 vector. avg_pc1 or pc1_fluctuation', required = True)
    parser.add_argument('-use_option', '--usage_option', help = 'use all samples or randomly picked samples. all or part', default = 'part', required = False) # all, part
    parser.add_argument('-n', '--normalize', help = 'whether you wnat to normalize score2 and score4 or not', required = True)#Y, N
    return parser.parse_args()


# score filenames: {score}_{score_type}_{usage_option}.csv
#     - score: score2 or score4
#     - score_type: pc1_fluctuation or avg_pc1
#     - usage_option: all or part

# In[9]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort+'_samples']
    if cohort=='PCBC':
        return samples.tolist()
    else: #TCGA cohort
        T = []
        N = []
        S = [] #all samples
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
        S = T + N
        print("{}: tumor {}, normal {}, total {}".format(cohort, len(T), len(N), len(S)))

        return T, N, S


# In[15]:


def compute_angle(x, y):
    # since x and y are both distances, all (x,y) shoud be located in the first quadrant
    if len(x) != len(y):
        raise ValueError
    slope = np.array([float(y[i]/x[i]) for i in range(len(x))])
    negative_slope_count = np.array([slope[i] < 0 for i in range(len(slope))]).sum()
    if negative_slope_count > 0:
        raise Exception("Negtive slope")
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


# In[75]:


def scatter_distances(normal_distance, stem_distance, final_df, T, N, cohort, score_type, figure_name):
    max_distance = max(np.max(normal_distance), np.max(stem_distance))
    fig = plt.figure(figsize = (3,3))
    ax1 = fig.add_subplot(111)
    ax1.scatter(final_df.loc[T].copy().normal_distance.values.flatten(), final_df.loc[T].copy().stem_distance.values.flatten(), label = 'Tumor')
    ax1.scatter(final_df.loc[N].copy().normal_distance.values.flatten(), final_df.loc[N].copy().stem_distance.values.flatten(), label = 'Normal')
    #ax1.plot(np.arange(0, round(max_distance)), np.arange(0, round(max_distance)) , label = 'y = x', color = 'k')
    ax1.set_xlabel('normal_distance')
    ax1.set_ylabel('stem_distance')
    ax1.set_title('Scatter plot of distances ({}, {})'.format(cohort, score_type))
    ax1.legend()
    fig.tight_layout()
    full_figname = os.path.join(SAVEDIR, figure_name)
    plt.savefig(full_figname)
    print("scatter plot: {}".format(full_figname))


if __name__ == '__main__':
    args = parse_arguments()
    
    if args.cohort not in NORMAL7_COHORT:
        sys.exit('This cohort has less than 7 normal samples. Cannot calculate normal_distance')
    else:
        T, N, S = get_sample_list(args.cohort) 

    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)


    # make SAVEDIR if it doesn't exist. 
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("1. normal_distance and stem_cell_distance")
    if args.score_type == 'avg_pc1':
        normal_distance = pd.read_csv(os.path.join(SAVEDIR, 'score2'+'_'+args.score_type+'_'+args.usage_option+'.csv'), index_col = 0).loc[S].copy().score2_avg_pc1.values.flatten()
        stem_distance = pd.read_csv(os.path.join(SAVEDIR, 'score4'+'_'+args.score_type+'_'+args.usage_option+'.csv'), index_col = 0).loc[S].copy().score4_avg_pc1.values.flatten()

    else:#args.score_type == 'pc1_fluctuation'
        normal_distance = pd.read_csv(os.path.join(SAVEDIR, 'score2'+'_'+args.score_type+'_'+args.usage_option+'.csv'), index_col = 0).loc[S].copy().values.flatten() #score2
        stem_distance = pd.read_csv(os.path.join(SAVEDIR, 'score4'+'_'+args.score_type+'_'+args.usage_option+'.csv'), index_col = 0).loc[S].copy().values.flatten() #score4
    
    if args.normalize=='Y':
        print("Normalize normal- and stem- distances.")
        normal_distance_max = float(NORMAL_STEM_DISTANCE_MINMAX.loc[args.cohort]['max_score2_'+args.score_type])
        stem_distance_max = float(NORMAL_STEM_DISTANCE_MINMAX.loc[args.cohort]['max_score4_'+args.score_type])
        normal_distance /= normal_distance_max
        stem_distance /= stem_distance_max
    
    # use raw score W/O any kind of scaling or normalization. 
    print("2. compute angle")
    radian_theta, degree_theta, cosine_radian = compute_angle(normal_distance, stem_distance) #input should be in order of x and y values
    final_df = pd.DataFrame(zip(normal_distance, stem_distance, radian_theta, degree_theta, cosine_radian), 
                            index = S, columns = ['normal_distance', 'stem_distance', 'angle_radian', 'angle_degree','cos_radian'])
                            #angle: (x, y) = (normal_distance, stem_distance)

    print("3. save result")
    normalize_yn = '_normalized' if args.normalize=='Y' else ''
    result_fname = os.path.join(SAVEDIR, 'stem-closeness_'+args.score_type+'_'+args.usage_option+normalize_yn+'.csv') 

    final_df.to_csv(result_fname)
    print('result file: {}'.format(result_fname))
    print("----")

    print("4. check whether tumor scores and normal scores differ significantly")
    score_df = pd.read_csv(result_fname, index_col = 0)
    tumor_mask = np.array([int(x[13:15])<=9 for x in score_df.index.values])
    tumor_score = score_df.iloc[tumor_mask,:].cos_radian.values.flatten()
    normal_score = score_df.iloc[~tumor_mask,:].cos_radian.values.flatten()
    print("---\nIndependent t-test between tumor and normal score")
    print(ttest_ind(tumor_score, normal_score))
    if ttest_ind(tumor_score, normal_score)[1]<=P_THRESHOLD:
        print("significant. p-value {} <= {}".format(ttest_ind(tumor_score, normal_score)[1], P_THRESHOLD))
    print("---\nTumor score mean and std")
    print("(mean, std) = ({}, {})".format(np.mean(tumor_score), np.std(tumor_score)))
    print("---\nNormal score mean and std")
    print("(mean, std) = ({}, {})".format(np.mean(normal_score), np.std(normal_score)))

    print("5. scatter plot (x: normal_distance, y: stem_distance")

    
    scatter_distances(normal_distance, stem_distance, final_df, T, N, args.cohort, args.score_type, 'scatter-normal-stem-distances_'+args.score_type+'_'+args.usage_option+normalize_yn+'.png')
    #scatter_distances(normal_distance, stem_distance, final_df, T, N, args.cohort, args.score_type, 'scatter-normal-stem-distances_'+args.score_type+'_'+args.usage_option+'.png')
    print('===')
    plt.clf()

