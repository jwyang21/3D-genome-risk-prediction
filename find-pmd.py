#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
import sys
import argparse


# In[ ]:


# To-Do
## PMD 기준을 만족하는 genomic bin들 찾기 -> 여기에 위치한 cpg probe들의 목록을 probemap 사용해서 찾기 (hg19) -> cpg probe 목록 저장.


# - PMD definition
#     - 10kb 단위 (longer: 100kb)
#     - 최소 10개 이상의 methylated CpG dinucleotide 포함
#     - average meth.level <= 0.7
#     - Tumor에서만.
# - Non-PMD
#     - 전체 genome들 중 PMD를 제외한 지역 (여집합)
# - PMD-like
#     - Non-tumor에서.
#     - PMD와 동일한 기준 적용

# In[ ]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/find-pmd')


# In[ ]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
#BETA_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA-LIHC_hg19/TCGA.LIHC.sampleMap%2FHumanMethylation450'
#BETA = pd.read_csv(BETA_FNAME, sep = '\t', index_col = 0)
SAVEDIR = os.path.join(os.getcwd(), 'result')
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
#Target columns: 'CHR', 'MAPINFO', 'UCSC_CpG_Islands_Name', 'Relation_to_UCSC_CpG_Island’


# In[ ]:


# CPG_ANNOT의 'CHR' column에 data type이 int, str 섞여있어서 정리함.
CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()].copy()
tmp = [str(str(x).split('.')[0].strip()) for x in CPG_ANNOT2['CHR'].values.flatten()]
CPG_ANNOT2['CHR'] = tmp


# In[ ]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True) #어떤 score로 구한 PC1으로 sample 간 거리를 계산하는지.
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #which distance metric you want to use
    parser.add_argument('-b', '--binsize', help = 'binsize', type = int, required = True)
    return parser.parse_args()


# In[ ]:


def get_sample_list(directory):
#def get_sample_list(directory):
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
    print("Tumor: {}, Normal: {}, Total: {}".format(len(T), len(N), len(S)))
    return T, N, S


# In[85]:


def get_pmd(beta, binsize, chr_len):
    # input: 현재 코호트의 raw beta data
    # output: cpg column과 chromosome column으로 구성된 pmd-cpg df
    print("Find CpG probes which belong to PMD regions.")
    
    cpg_dictionary={}

    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()]
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 

    for chrom in CHR_LIST:#for real implementation
    #for chrom in CHR_LIST[:2]:#for test
        print(chrom)
        chrom_index = str(chrom.split('chr')[-1])
        
        tmp_list=np.array([], dtype=str)
        cpg_annot_chrom = CPG_ANNOT2[CPG_ANNOT2['CHR']==chrom_string][['Name', 'CHR', 'MAPINFO', 'UCSC_CpG_Islands_Name', 'Relation_to_UCSC_CpG_Island']]
        #display(cpg_annot_chrom.head(3))#debug
        chrlen = int(CHR_LENGTH[CHR_LENGTH['Chromosome']==chrom].Total_length.values[0])
        #print("chrlen: {}".format(chrlen))#debug
        n_bins = (chrlen // binsize) + 1
        #print("n_bins: {}".format(n_bins))#debug
        for i in range(n_bins):#for real implementation #현재 chromosome의 총 길이를 10kb 단위로 sliding window하면서
        #for i in range(10):#debug
            #print("----")#debug
            if i % 500 == 0:
                print('{}-th bin'.format(i+1))
            cur_window_start = i * binsize
            cur_window_end = (i+1) * binsize
            window_cpg_mask = [int(x)>=int(cur_window_start) and int(x)<int(cur_window_end) for x in cpg_annot_chrom.MAPINFO.values.flatten()]
            #print("current window: from {} to {}".format(cur_window_start, cur_window_end))#debug
            cur_cpg_annot_window = cpg_annot_chrom.iloc[window_cpg_mask, :]
            #display(cur_cpg_annot_window.head(3))#debug
            #print("cur_cpg_annot_window.shape: {}".format(cur_cpg_annot_window.shape))#debug
            #beta_extracted = beta.loc[np.intersect1d(cur_cpg_annot_window.Name.values.flatten(), beta.index.values.flatten())]
            intersecting_cpgs = np.intersect1d(cur_cpg_annot_window.Name.values.flatten(), beta.index.values.flatten())
            #print("intersecting cpgs: {}".format(intersecting_cpgs))#debug
        
            if len(intersecting_cpgs)>=10: #최소 10개 이상의 methylated CpG dinucleotide를 가져야 함.
                beta_extracted = beta.loc[intersecting_cpgs]
                #print("Mean DNA methylation: {}".format(beta_extracted.mean().mean()))#debug
                if beta_extracted.mean().mean() < 0.7: #average methylation level이 0.7보다 낮아야 함.
                    tmp_list = np.concatenate((tmp_list, intersecting_cpgs), axis = None)
                else:
                    pass
            else:
                pass
        cpg_dictionary[chrom] = tmp_list
        
    key_list = []
    value_list = []
    for i in range(len(list(cpg_dictionary.keys()))):
        current_key = list(cpg_dictionary.keys())[i]
        current_value = cpg_dictionary[current_key]
        for j in range(len(current_value)):
            value_list.append(current_value[j])
            key_list.append(current_key)
    if len(key_list) != len(value_list):
        raise ValueError
    pmd_cpg = pd.DataFrame(zip(key_list, value_list), columns = ['chrom', 'cpg'])
    
    return pmd_cpg


# In[ ]:


if __name__=='__main__':
    
    #for cohort in COHORTS:
    args = parse_arguments()
    cohort_dir = os.path.join(SAVEDIR, args.cohort)  #현재 cohort의 결과가 이 디렉토리에 저장돼야 함. 
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    print("cohort: {}".format(args.cohort))
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    
    print("Import beta value.")
    beta = pd.read_csv(os.path.join(METH_DIR, args.cohort+'.HumanMethylation450.tsv'), sep = '\t', index_col=0) #row: cpg, column: TCGA sample
    
    pmd_cpg = get_pmd(beta, args.binsize, CHR_LENGTH)  #tumor sample에서 pmd에 속하는 region에 위치한 cpg들의 list. 
    
    # save results
    fname = os.path.join(cohort_dir, 'pmd_cpg.csv')
    pmd_cpg.to_csv(fname, index=False)
    print(fname)
    print("---\n")

