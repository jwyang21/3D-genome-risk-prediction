#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pandas as pd
import numpy as np
#from collections import defaultdict
#from scipy.integrate import simps
import os
import random
random.seed(2022)
np.random.seed(2022)
import sys

# # description
# - normal reference를 계산
# - avg_pc1 및 pc1_fluctuation의 두 가지 버전의 reference를 계산.
#     - avg_pc1: pc1 vector들의 평균 -> average pc1 vector를 만들고 sample-of-interest의 pc1 vector와 reference pc1 vector 간의 euclidean distance를 계산
#         - pc1 vector는 chromosome마다 존재하기 때문에, 한 샘플당 총 22개의 값이 생김 -> 22개 값들을 평균냄. -> scalar 값이 됨. 이게 normal_avg_pc1_distance
#     - pc1_fluctuation: 각 normal sample들의 각 chromosome의 pc1 vector를 절댓값 적분 -> 전체 normal sample들 중 절반만 골라서 각 chromosome 별로 값을 평균냄
#         - reference는 (1,22) shape의 vector가 됨.
#         - sample-of-interest도 (1, 22) shape의 vector를 가질 것이므로, 샘플의 벡터와 reference vector 간 euclidean distance 계산하면 scalar 값이 나옴. 이게 normal_pc1_fluctuation_distance
# - reference 계산 시 현재 cohort의 normal sample들 중 절반만 랜덤하게 샘플링.

# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-reference/')


# In[48]:


CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
# TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
NORMAL_BARCODES = ['10', '11', '12', '13','14', '15', '16', '17', '18', '19']
## TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
#CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
#CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{cohort}/{sample}.npz or /{cohort}/{sample}_inv_exp.npz  #'chr1'
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
TCGA_SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/score3_simple_avg.pickle
PCBC_SCORE3_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv'
SAVEDIR = os.path.join(os.getcwd(), 'result')
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
ALL_COHORT_W_PCBC = ALL_COHORT.copy()
ALL_COHORT_W_PCBC.append('PCBC')
print("SAVEDIR: {}".format(SAVEDIR))


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True) #'TCGA-{}' or 'PCBC'
    #parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)
    parser.add_argument('-r_type', '--reference_type', help = 'reference type. PCBC or TCGA', required = True) 
    parser.add_argument('-s_type', '--score_type', help = 'score type', required = True) #avg_pc1 #pc1_fluctuation
    return parser.parse_args()


# In[18]:


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


# In[6]:


def import_pc1(cohort, sample, chrom, flag, ref_type):
    # import pre-computed PC1 of sample-of-interest
    # ref_type: 'TCGA' or 'PCBC'. Type of reference you want to import. 
    # flag: 'raw' or 'inv'
    if ref_type=='TCGA':
        if flag=='raw':
            fname = os.path.join(TCGA_PC1_DIR, cohort, sample+'.npz')
        elif flag=='inv':
            fname = os.path.join(TCGA_PC1_DIR, cohort, sample+'_inv_exp.npz')
        else:
            pass

    elif ref_type == 'PCBC':
        if flag=='raw':
            fname = os.path.join(PCBC_PC1_DIR, sample+'.npz')
        elif flag=='inv':
            fname = os.path.join(PCBC_PC1_DIR, sample+'_inv_exp.npz')
        else:
            pass
        
    pc1 = np.load(fname)[chrom]
    return pc1




if __name__=='__main__':
    
    args = parse_arguments()
    
    print("cohort: {}".format(args.cohort))
    if args.cohort=='PCBC':
        S = get_sample_list(args.cohort)
    else: #TCGA cohort
        T, N, S = get_sample_list(args.cohort)

    # TCGA 코호트이rh 현재 코호트의 normal이 7개 이상이면
    if 'TCGA' in args.cohort and len(N) >= 7:
        # 현재 코호트의 normal들 중 절반을 랜덤샘플링해서 sample name list 만들기. 
        sampled_N = random.sample(N, len(N)//2)

    # PCBC이면
    elif args.cohort=='PCBC':
        sampled_N = random.sample(S, len(S)//2)
    else:
        sys.exit('Stop execution! Number of normal samples in this cohort is smaller than 7\n----') #stop execution.
    print("number of samples used to compute reference vector: {}".format(len(sampled_N)))
    
    # make saving directories if needed. 
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)

    print("SAVEDIR: {}".format(SAVEDIR))    
    
    # reference vector 계산. 
    # avg pc1 vector 구할 거면
    if args.score_type == 'avg_pc1':
        print("Calculate average pc1 vectors for each chromosome")
        # import_pc1 함수 써서, 위에서 구한 random-sampled normal들에 대한 각 chromosome별 pc1 vector를 import.
        for chrom in CHR_LIST:
            for s in S:              
                current_pc1 = import_pc1(args.cohort, s, chrom, 'inv', args.reference_type)
                if S.index(s)==0:#if this is the first sample
                    globals()['N_pc1_'+chrom] = current_pc1.copy()
                else:
                    globals()['N_pc1_'+chrom] = np.vstack((globals()['N_pc1_'+chrom], current_pc1))
            # chromosome 별로 random-sampled normal들의 PC1 vector들을 평균내서, chromosome 별 reference avg pc1 vector 만들기.
            

                
            globals()['N_pc1_'+chrom+'_ref'] = globals()['N_pc1_'+chrom].mean(axis=0)
        # output: 각 chromosome 별 averaged-pc1-vector
        # 저장
        save_fname = os.path.join(SAVEDIR, 'sampled_N_avg_pc1')
        np.savez(save_fname, chr1 = globals()['N_pc1_chr1_ref'], chr2 = globals()['N_pc1_chr2_ref'], chr3 = globals()['N_pc1_chr3_ref'], chr4 = globals()['N_pc1_chr4_ref'], 
                chr5 = globals()['N_pc1_chr5_ref'], chr6 = globals()['N_pc1_chr6_ref'], chr7 = globals()['N_pc1_chr7_ref'], chr8 = globals()['N_pc1_chr8_ref'], chr9 = globals()['N_pc1_chr9_ref'],
                chr10 = globals()['N_pc1_chr10_ref'], chr11 = globals()['N_pc1_chr11_ref'], chr12 = globals()['N_pc1_chr12_ref'], chr13 = globals()['N_pc1_chr13_ref'], chr14 = globals()['N_pc1_chr14_ref'],
                chr15 = globals()['N_pc1_chr15_ref'], chr16 = globals()['N_pc1_chr16_ref'], chr17 = globals()['N_pc1_chr17_ref'], chr18 = globals()['N_pc1_chr18_ref'], chr19 = globals()['N_pc1_chr19_ref'],
                chr20 = globals()['N_pc1_chr20_ref'], chr21 = globals()['N_pc1_chr21_ref'], chr22 = globals()['N_pc1_chr22_ref'])
        print("result file: {}".format(save_fname+'.npz'))
        print("----")


    # pc1 fluctuation 구할 거면
    elif args.score_type == 'pc1_fluctuation':
        print("Compute pc1 fluctuation for each chromosome")
        # score3 import (pc1 vector를 절댓값 적분한 것)
        if args.reference_type == 'TCGA':
            score3_fname = os.path.join(TCGA_SCORE3_DIR, args.cohort, 'score3.pickle')
            score_df = pd.read_pickle(score3_fname)#row: sample  #column: chromosome
        elif args.reference_type == 'PCBC':
            score_df = pd.read_csv(PCBC_SCORE3_FILE, index_col = 0)
        sampled_score_df = score_df.loc[sampled_N]
        # 각 chromosome 별로 random-sampled normal들의 값들을 average해서 1개의 scalar 값을 얻음
        globals()['N_pc1_fluctuation_avg'] = pd.DataFrame(sampled_score_df.mean(axis=0).values.flatten(), index = CHR_LIST, columns = ['avg_pc1_fluctuation'])
        # output: 현재 cohort에 대한 (1, num_autosome) 크기의 1d array
        # 저장  
        save_fname = os.path.join(SAVEDIR, 'sampled_N_avg_pc1_fluctuation.csv')
        globals()['N_pc1_fluctuation_avg'].to_csv(save_fname, index=True)
        print('result file: {}'.format(save_fname))
        print("----")

    else:
        pass

