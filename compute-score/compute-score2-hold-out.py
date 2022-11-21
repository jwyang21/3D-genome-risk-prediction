#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8


# - update on 22/11/21: reference 계산 시 전체 normal sample들을 다 쓰지 말고, normal들 중 50%를 hold out & normal sample 개수가 7개 이상인 cohort만 쓰기.
# - fix score2_reference and score2_distance as normal and euclidean, respectively.
# - use NORMAL7_COHORT

# In[22]:

# In[1]:


import numpy as np
import pandas as pd
import os
import sys
#from scipy.integrate import simps
#from numpy import trapz
from scipy.spatial.distance import euclidean
#from scipy.spatial.distance import jensenshannon
import argparse
#from numpy import dot
#from numpy.linalg import norm
import pickle
from scipy.stats import ttest_ind


# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# In[46]:


# global variables
CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
# TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
#NORMAL_BARCODES = ['10', '11', '12', '13','14', '15', '16', '17', '18', '19']
#CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
#CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{cohort}/{sample}.npz or /{cohort}/{sample}_inv_exp.npz  #'chr1'
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
#METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
#PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
#FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
#NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
#ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
#TCGA_SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/'#{cohort}/score3_simple_avg.pickle
#PCBC_SCORE3_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv
SAVEDIR = os.path.join(os.getcwd(), 'result')
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
REFERENCE_S2_S4_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-reference-v2/result/'#/{cohort}/sampled_N_avg_pc1_fluctuation.csv or sampled_N_avg_pc1.npz
print("SAVEDIR: {}".format(SAVEDIR))
P_THRESHOLD = 5e-2


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True)
    #parser.add_argument('-s4d', '--score4_distance', help = 'distance metric to be used in score4. euclidean, jsd, or cosine-similarity', defulat = 'euclidean', required = False)
    #parser.add_argument('-s4r', '--score4_reference', help = 'reference to be used in score4. all, sc, nonsc', default = 'all', required = False)
    parser.add_argument('-ref_type', '--reference_type', help = 'reference to be used in score2 or score4. TCGA or PCBC', required = True)
    parser.add_argument('-s_type', '--score_type', help = 'which one you will use among PC1 fluctuation or averaged PC1 vector. avg_pc1 or pc1_fluctuation', required = True)
    return parser.parse_args()


# In[6]:


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



# In[ ]:


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


# In[47]:


'''
# test in jupyter first

# initialize parameterws
#cohort, score2_reference, score2_distance, reference_type, score_type = 'TCGA-LIHC', 'normal', 'euclidean', 'TCGA', 'avg_pc1'
cohort, score2_reference, score2_distance, reference_type, score_type = 'TCGA-LIHC', 'normal', 'euclidean', 'TCGA', 'pc1_fluctuation'

if cohort not in NORMAL7_COHORT:
    sys.exit('This cohort has less than 7 normal samples. Cannot calculate score2')

print("cohort: {}".format(cohort))
SAVEDIR = os.path.join(os.getcwd(), 'result', cohort)

if not os.path.exists(os.path.join(os.getcwd(), 'result')):
    os.makedirs(os.path.join(os.getcwd(), 'result'))
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

T, N, S = get_sample_list(cohort) 


if score_type == 'avg_pc1':
    # 이 코호트의 모든 샘플들의 각 chromosome별 score 및 각 샘플마다 22개 score들의 평균값을 저장. 
    all_score_df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), index = S, columns = CHR_LIST)
    for s in S: #repeat for every sample
        score_list = [] #각 chromosome의 score가 순차적으로 저장되게 # 길이 22의 1-dimensional vector.
        for chrom in CHR_LIST: #repeat for every chromosome
            sample_pc1 = import_pc1(cohort, s, chrom, 'inv', reference_type).flatten()
            ref_pc1_fname = os.path.join(REFERENCE_S2_S4_DIR, cohort, 'sampled_N_avg_pc1.npz')
            ref_pc1 = np.load(ref_pc1_fname)[chrom].flatten()
            score = euclidean(sample_pc1, ref_pc1)
            score_list.append(score)
        all_score_df.loc[s] = score_list
    simple_avg = all_score_df.mean(axis = 1).values #rowmean
    all_score_df['score2_avg_pc1'] = simple_avg

elif score_type == 'pc1_fluctuation':
    all_score_df = pd.DataFrame(np.zeros((len(S), 1), dtype = float), index = S, columns = ['score2_pc1_fluctuation'])
    for s in S:
        sample_vector = pd.read_pickle('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-LIHC/score3.pickle').loc[s].values.flatten() #import score3 #길이 22의 1d vector
        ref_vector_fname = os.path.join(REFERENCE_S2_S4_DIR, cohort, 'sampled_N_avg_pc1_fluctuation.csv')
        ref_vector = pd.read_csv(ref_vector_fname, index_col = 0).values.flatten()
        score = euclidean(sample_vector, ref_vector)
        all_score_df.loc[s] = score

else: #score type is neither avg_pc1 nor pc1_fluctuation
    sys.exit('invalid score_type')


all_score_df_fname = os.path.join(SAVEDIR, 'score2_'+score_type+'.csv')#assume normal reference
all_score_df.to_csv(all_score_df_fname, index = True)
print("all_score_df (score2): {}".format(all_score_df_fname))
print("----")

print("Check if scores of tumor and normal differs significantly.")
if score_type == 'avg_pc1':
    tumor_score = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-LIHC/score2_avg_pc1.csv', index_col = 0).loc[T].score2_avg_pc1.values.flatten()
    normal_score = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-LIHC/score2_avg_pc1.csv', index_col = 0).loc[N].score2_avg_pc1.values.flatten()
elif score_type == 'pc1_fluctuation':
    tumor_score = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-LIHC/score2_pc1_fluctuation.csv', index_col = 0).loc[T].values.flatten()#score2_avg_pc1.values.flatten()
    normal_score = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-LIHC/score2_pc1_fluctuation.csv', index_col = 0).loc[N].values.flatten()#score2_avg_pc1.values.flatten()
else:
    sys.exit(0)
print("(1) Independent t-test (tumor score vs. normal score)")
print(ttest_ind(tumor_score, normal_score))
if ttest_ind(tumor_score, normal_score)[1] <= P_THRESHOLD:
    print("significant! p-value <= {}".format(P_THRESHOLD))
print("----")
print("(2) Mean and Std of tumor scores and normal score")
print("(tumor mean, tumor std): ({}, {})".format(np.mean(tumor_score), np.std(tumor_score)))
print("(normal mean, normal std): ({}, {})".format(np.mean(normal_score), np.std(normal_score)))
print("----")
'''


# In[ ]:


if __name__ == '__main__':
    args = parse_arguments()
    
    if 'TCGA' in args.cohort:
        if args.cohort not in NORMAL7_COHORT:
            sys.exit('This cohort has less than 7 normal samples. Cannot calculate score2')
        else:
            T, N, S = get_sample_list(args.cohort) 
    elif args.cohort=='PCBC':
        S = get_sample_list(args.cohort)

    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)


    if args.score_type == 'avg_pc1':
        # 이 코호트의 모든 샘플들의 각 chromosome별 score 및 각 샘플마다 22개 score들의 평균값을 저장. 
        all_score_df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), index = S, columns = CHR_LIST)
        for s in S: #repeat for every sample
            score_list = [] #각 chromosome의 score가 순차적으로 저장되게 # 길이 22의 1-dimensional vector.
            for chrom in CHR_LIST: #repeat for every chromosome
                sample_pc1 = import_pc1(args.cohort, s, chrom, 'inv', args.reference_type).flatten()#reference_type: 'TCGA' or 'PCBC'
                ref_pc1_fname = os.path.join(REFERENCE_S2_S4_DIR, args.cohort, 'sampled_N_avg_pc1.npz')
                ref_pc1 = np.load(ref_pc1_fname)[chrom].flatten()
                score = euclidean(sample_pc1, ref_pc1)
                score_list.append(score)
            all_score_df.loc[s] = score_list
        simple_avg = all_score_df.mean(axis = 1).values #rowmean
        all_score_df['score2_avg_pc1'] = simple_avg

    elif args.score_type == 'pc1_fluctuation':
        all_score_df = pd.DataFrame(np.zeros((len(S), 1), dtype = float), index = S, columns = ['score2_pc1_fluctuation'])
        for s in S:
            if 'TCGA' in args.cohort:
                sample_vector = pd.read_pickle(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/',args.cohort, 'score3.pickle')).loc[s].values.flatten() #import score3 #길이 22의 1d vector
            elif 'PCBC' in args.cohort:
                sample_vector = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv', index_col = 0).loc[s].values.flatten()
            else:
                sys.exit("wrong cohort.")
            ref_vector_fname = os.path.join(REFERENCE_S2_S4_DIR, args.cohort, 'sampled_N_avg_pc1_fluctuation.csv')
            ref_vector = pd.read_csv(ref_vector_fname, index_col = 0).values.flatten()
            score = euclidean(sample_vector, ref_vector)
            all_score_df.loc[s] = score

    else: #score type is neither avg_pc1 nor pc1_fluctuation
        sys.exit('invalid score_type')


    all_score_df_fname = os.path.join(SAVEDIR, 'score2_'+args.score_type+'.csv')#assume normal reference
    all_score_df.to_csv(all_score_df_fname, index = True)
    print("all_score_df (score2): {}".format(all_score_df_fname))
    print("----")
    
    if 'TCGA' in args.cohort: #conduct independent t-test and comparing mean/std value of tumor and normal scores only in case of TCGA cohorts!

        print("Check if scores of tumor and normal differs significantly.")
        if args.score_type == 'avg_pc1':
            tumor_score = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', args.cohort, 'score2_avg_pc1.csv'), index_col = 0).loc[T].score2_avg_pc1.values.flatten()
            normal_score = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', args.cohort, 'score2_avg_pc1.csv'), index_col = 0).loc[N].score2_avg_pc1.values.flatten()
        elif args.score_type == 'pc1_fluctuation':
            tumor_score = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', args.cohort, 'score2_pc1_fluctuation.csv'), index_col = 0).loc[T].values.flatten()#score2_avg_pc1.values.flatten()
            normal_score = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', args.cohort, 'score2_pc1_fluctuation.csv'), index_col = 0).loc[N].values.flatten()#score2_avg_pc1.values.flatten()

        else:
            sys.exit(0)
        print("(1) Independent t-test (tumor score vs. normal score)")
        print(ttest_ind(tumor_score, normal_score))
        if ttest_ind(tumor_score, normal_score)[1] <= P_THRESHOLD:
            print("significant! p-value <= {}".format(P_THRESHOLD))
        print("----")
        print("(2) Mean and Std of tumor scores and normal score")
        print("(tumor mean, tumor std): ({}, {})".format(np.mean(tumor_score), np.std(tumor_score)))
        print("(normal mean, normal std): ({}, {})".format(np.mean(normal_score), np.std(normal_score)))
        print("----")

