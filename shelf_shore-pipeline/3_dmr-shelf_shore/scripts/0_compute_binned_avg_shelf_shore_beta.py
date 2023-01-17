#!/usr/bin/env python
# coding: utf-8

# command: python3 0_compute_binned_avg_island_beta.py -w_dir /data/project/3dith/pipelines/island-pipeline/3_dmr-island --cpg_type island > ../log/0_compute_binned_avg_island_beta.log

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import argparse


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
#SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    args.add_argument('--cpg_type', help = 'CpG type. open/island/shelf/shore/shelf_shore', type = str, required = True)
    return args.parse_args()


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
    
def get_cpg_list(chr_list, cpg_type):
    if cpg_type == 'opensea':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', sep = '\t', header = None) #chrom // start // end // cpg_id
        '''
        # previous version
        total_list = np.array([])
        for chrom in chr_list: 
            #opensea version 
            # # opensesa cpg_list_fname: '/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+chrom+'_opensea_CpG.pickle'
            fname = '/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+chrom+'_opensea_CpG.pickle'
            cpgs = pd.read_pickle(fname).index.values
            total_list = np.union1d(total_list, cpgs)
        '''
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


def get_avg_beta(cohort, cpg_list, S):
    # 현재 cohort의 전체 beta 데이터 중에서, 입력받은 cpg_list들의 beta value들만 반환.
    beta_fname = '/data/project/3dith/data/450k_xena/'+cohort+'.HumanMethylation450.tsv' #이 파일로부터 beta value를 읽어옴. 
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

# jupyter notebook에서 test
#score='stem_closeness'#for test
#normalize='N'
#minmax='Y'
#score_type == ''
#args: cohort, score_type, S, stemness_type, normalize, minmax
if __name__ == '__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    SAVEDIR = os.path.join(os.getcwd(), 'result')
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    for cohort_ in NORMAL7_COHORT:#real
    #for cohort_ in NORMAL7_COHORT[:1]:#for test (debug)
        cohort=cohort_ #for test
        print("====")

        print("1. get sample name list")
        cohort_dir = os.path.join(SAVEDIR, cohort)  #현재 cohort의 결과가 이 디렉토리에 저장돼야 함. 
        print(cohort_dir)
        print("cohort: {}".format(cohort))
        print("cohort_dir: ", cohort_dir)
        if not os.path.exists(cohort_dir):
            os.makedirs(cohort_dir)
        T, N, S = get_sample_list(cohort) 


        #/data/project/3dith/data/450k_xena/TCGA-[XXXX].HumanMethylation450.tsv

        print("----\n2. find {} cpg probes of 22 autosomes".format(args.cpg_type))
        # 먼저, 22개 염색체들의 {cpg_type} CpG probe ID들의 합집합을 구하기 
        cpg_list = get_cpg_list(CHR_LIST, args.cpg_type) 
        print("num_{}_cpg: {}".format(args.cpg_type, len(cpg_list)))
        
        print("----\n3. import beta value corresponding to current cpg_type and average them")
        # 현재 cohort의 beta 값 불러들이기
        beta_target_cpg_df, avg_beta_target_cpg_df = get_avg_beta(cohort, cpg_list, S) #beta_target_cpg_df: column이 sample, index가 CpG -> 여기서 column mean 
                                                                                #-> avg_beta_target_cpg_df: index가 sample, column이 average target_cpg beta value.

        beta_target_cpg_df_fname = 'beta_'+args.cpg_type+'_df.csv'
        avg_beta_target_cpg_df_fname = 'avg_beta_'+args.cpg_type+'_df.csv'

        beta_target_cpg_df.to_csv(os.path.join(cohort_dir, beta_target_cpg_df_fname), index = True)
        avg_beta_target_cpg_df.to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname), index = True)
        print(os.path.join(cohort_dir, beta_target_cpg_df_fname))
        print(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname))
