#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# command: bash 2_pcc-avg_beta-stem_closeness.sh > ../log/2_pcc-avg_beta-stem_closeness-ALL.log

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import argparse


# 각 cohort마다, 22개 autosome에 대한 target cpg type의 CpG probe들의 합집합들을 구함 -> 전체 beta value데이터에서, 각 샘플마다 이 CpG probe들의 데이터만 추출한 후 

# In[ ]:


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[ ]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


# In[ ]:


def parse_arguments():
    args = argparse.ArgumentParser()
    
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    #/data/project/3dith/pipelines/shelf_shore-pipeline
    
    args.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    # 일단 TCGA-BLCA
    
    args.add_argument('--score_fname', help = 'stem-closeness score file name, based on shelf_shore K-M result', type = str, default = '/data/project/3dith/data/cohort-1-best-score-km.csv', required = True)

    args.add_argument('--cpg_type', help = 'open/island/shelf_shore', type = str, required = True)
    return args.parse_args()


# In[ ]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = []
    else: #TCGA cohort
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S
    # usage: T, N, S = get_sample_list(args.cohort)


# In[ ]:


def get_cpg_list(chr_list, cpg_type):
    if cpg_type == 'shelf_shore':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', sep = '\t', header = None) #chrom // start // end // cpg_id

    elif cpg_type == 'island':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.island.sorted.bed', sep = '\t', header = None) #chrom // start // end // cpg_id

    elif cpg_type == 'shelf_shore':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.shelf_shore.sorted.bed', sep = '\t', header = None) #chrom / start / end / cpg_id
    else: #일단 shelf, shore, shelr_shore는 보류.
        raise Exception("Wrong cpg_type!")

    chr_list_mask = np.array([x in chr_list for x in df.iloc[:,0].values.flatten()])
    df_chr_list = df.iloc[chr_list_mask,:]

    for x in df_chr_list.iloc[:,0].values.flatten():
        assert x in chr_list

    total_list = df_chr_list.iloc[:,-1].values.flatten()
    return total_list


# In[ ]:


def get_avg_beta(cohort, cpg_list, S): # get average shelf_shore or resort CpG beta values (per sample)
    # 현재 cohort의 전체 beta 데이터 중에서, 입력받은 cpg_list들의 beta value들만 반환.
    beta_fname = '/data/project/3dith/data/450k_xena/'+cohort+'.HumanMethylation450.tsv' #이 파일로부터 beta value를 읽어옴. #cohort의 전체 beta value 값.  
    beta = pd.read_csv(beta_fname, sep = '\t', index_col = 0)
    
    beta_target_cpg_df = beta.loc[cpg_list].dropna() #beta 데이터: row가 cpg, column이 samplename
    print("num_target_cpg_after_dropna: {}".format(beta_target_cpg_df.shape[0]))
    samples_beta = beta_target_cpg_df.columns.values
    avg_target_cpg_beta = beta_target_cpg_df.mean().values #column mean
    
    avg_beta_target_cpg_df = pd.DataFrame(avg_target_cpg_beta, index = samples_beta, columns = ['avg_beta'])
    
    #avg_beta_target_cpg_df = pd.DataFrame(zip(samples_beta, avg_target_cpg_beta), columns = ['sample', 'avg_beta'])
    #print("avg_beta_target_cpg_df")#debug
    #print(avg_beta_target_cpg_df)#debug #index가 s
    return beta_target_cpg_df[S], avg_beta_target_cpg_df.loc[S] 


if __name__=='__main__':
    #for cohort in COHORTS:
    args = parse_arguments()
    
    os.chdir(args.working_dir)
    result_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'result')
     

    cohort_dir = os.path.join(result_dir, args.cohort)  #현재 cohort의 결과가 이 디렉토리에 저장돼야 함.
    if not os.path.exists(result_dir):
        os.makedirs(result_dir) 
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    print("cohort: {}".format(args.cohort))
    T, N, S = get_sample_list(args.cohort)

    # 먼저, 22개 염색체들의 target_cpg (shelf_shore OR resort CpG probes) 들의 합집합을 구하기 
    cpg_list = get_cpg_list(CHR_LIST, args.cpg_type) 
    print("num_{}_cpg (for 22 autosomes): {}".format(args.cpg_type, len(cpg_list)))

    # 현재 cohort의 beta value 전체를 불러들여서, 위에서 구한 CpG probe들만 남긴 후 개별 샘플의 평균 beta value 계산. 
    beta_target_cpg_df, avg_beta_target_cpg_df = get_avg_beta(args.cohort, cpg_list, S) #beta_target_cpg_df: column이 sample, index가 CpG -> 여기서 column mean 
                                                                            #-> avg_beta_target_cpg_df: index가 sample, column이 average target_cpg beta value.

    # beta_target_cpg_df: (num_{cpg_type}_CpGs, num_samples)
    # avg_beta_target_cpg_df: (num_samples, )
    beta_target_cpg_df_fname = args.cpg_type+'_all_samples_all_beta.csv'
    avg_beta_target_cpg_df_fname = args.cpg_type+'_all_samples_avg_beta.csv'

    beta_target_cpg_df_Tumor_fname = args.cpg_type+'_tumors_all_beta.csv'
    avg_beta_target_cpg_df_Tumor_fname = args.cpg_type+'_tumors_avg_beta.csv'

    beta_target_cpg_df.to_csv(os.path.join(cohort_dir, beta_target_cpg_df_fname))
    avg_beta_target_cpg_df.to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname))
    
    # test in jupyter
    if 'TCGA' in args.cohort:
        beta_target_cpg_df[T].copy().to_csv(os.path.join(cohort_dir, beta_target_cpg_df_Tumor_fname))
        avg_beta_target_cpg_df.loc[T].copy().to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_Tumor_fname))
        print(os.path.join(cohort_dir, beta_target_cpg_df_Tumor_fname))
        print(os.path.join(cohort_dir, avg_beta_target_cpg_df_Tumor_fname))
    else:
        pass
    
    # ESCA, HNSC는 score 값이 없음. kaplan-meier 결과가 안 좋기 때문. best score이 1개 event에 대해서만 significant하기 때문. 
    if 'ESCA' in args.cohort:
        sys.exit(0)
    elif 'HNSC' in args.cohort:
        sys.exit(0)
    else:
        # 현재 cohort의 각 sample의 score 값 import
        all_score_df = pd.read_csv(args.score_fname, index_col = 0)
        current_score_fname = str(all_score_df.loc[args.cohort][f'filename_{args.cpg_type}'])
        score_df = pd.read_csv(current_score_fname, index_col = 0)
        assert args.cohort in current_score_fname

        # score과 avg beta value 간 PCC
        ## 전체 samples
        print("---\nPCC(avg_beta({}), stem-closeness) in {}, using all samples (Tumor + Normal)".format(args.cpg_type, args.cohort))
        pcc, pvalue = pearsonr(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten())
        print("pcc: {}, p-value: {}".format(pcc, pvalue))

        print("---\nPCC(avg_beta({}), stem-closeness) in {}, using tumor samples".format(args.cpg_type, args.cohort))
        pcc, pvalue = pearsonr(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten())
        print("pcc: {}, p-value: {}".format(pcc, pvalue))

        # scatter (avg_beta, stem-closeness) #all samples
        fig = plt.figure(figsize = (3,3))
        ax = fig.add_subplot(111)
        ax.scatter(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten())
        ax.set_xlabel('average {} beta value'.format(args.cpg_type))
        ax.set_ylabel('stem-closeness')

        ax.set_title(f'{args.cohort} (all samples)')
        fig.tight_layout()
        plt.savefig(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.png'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.png'))
        figure_data = pd.DataFrame(zip(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten()), index = S, columns = ['avg_beta', 'stem_closeness'])
        figure_data.to_csv(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.csv'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.csv'))

        # scatter (avg_beta, stem-closeness) #tumor smaples
        fig = plt.figure(figsize = (3,3))
        ax = fig.add_subplot(111)
        ax.scatter(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten())
        ax.set_xlabel('average {} beta value'.format(args.cpg_type))
        ax.set_ylabel('stem-closeness')
        ax.set_title(f'{args.cohort} (tumor samples)')
        plt.savefig(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.png'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.png'))
        figure_data = pd.DataFrame(zip(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten()), index = T, columns = ['avg_beta', 'stem_closeness'])
        figure_data.to_csv(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.csv'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.csv'))

        plt.clf
