#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8


# In[1]:


import pandas as pd
import numpy as np
import os
import sys
import argparse


# # To-Do
#  PMD 기준을 만족하는 genomic bin들 찾기 -> 여기에 위치한 cpg probe들의 목록을 probemap 사용해서 찾기 (hg19) -> cpg probe 목록 저장.

# # Definitions
# - PMD definition<br>
#     - 10kb 단위 (longer: 100kb)<br>
#     - 최소 10개 이상의 methylated CpG dinucleotide 포함<br>
#     - average meth.level <= 0.7<br>
#     - Tumor에서만.<br>
# - Non-PMD<br>
#     - 전체 genome들 중 PMD를 제외한 지역 (여집합)<br>
# - PMD-like<br>
#     - Non-tumor에서.<br>
#     - PMD와 동일한 기준 적용

# # output:
# - pmd bed file
# - pmd cpg list

# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/find-pmd')


# In[3]:


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
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples


# In[4]:


#CPG_ANNOT의 'CHR' column에 data type이 int, str 섞여있어서 정리함.
CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()].copy()
tmp = [str(str(x).split('.')[0].strip()) for x in CPG_ANNOT2['CHR'].values.flatten()]
CPG_ANNOT2['CHR'] = tmp


# In[5]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True) #어떤 score로 구한 PC1으로 sample 간 거리를 계산하는지.
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #which distance metric you want to use
    parser.add_argument('-b', '--binsize', help = 'binsize', type = int, required = False, default = int(1e6))
    return parser.parse_args()


# In[ ]:

# In[5]:


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


def get_pmd(beta, binsize, chr_len):
    # input: 현재 cohort의 tumor 또는 normal의 beta value 
    # output: pmd bedfile dataframe, pmd region에 속하는 cpg들의 목록.

    #print("Find CpG probes which belong to PMD regions.")
    pmd_cpg_dictionary={}#chrom 단위로 만들기. 
    non_pmd_cpg_dictionary={}
    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()] # '1', '10', '11', ...
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 

    pmd_df_data_chrom = []
    pmd_df_data_start = []
    pmd_df_data_end = []
    pmd_df_data_cpg = []
    non_pmd_df_data_chrom = []
    non_pmd_df_data_start = []
    non_pmd_df_data_end = []
    non_pmd_df_data_cpg = []

    for chrom in CHR_LIST:#for real 
    #for chrom in CHR_LIST[:1]:#for debug
        #print("===\nchrom: {}".format(chrom))
        chrom_index = str(chrom.split('chr')[-1])#'1', '10', '11', ...

        chrom_pmd_cpg_list = np.array([], dtype = str)
        chrom_non_pmd_cpg_list = np.array([], dtype = str)
        cpg_annot_chrom = CPG_ANNOT2[CPG_ANNOT2['CHR']==chrom_string][['Name', 'CHR', 'MAPINFO', 'UCSC_CpG_Islands_Name', 'Relation_to_UCSC_CpG_Island']]
        #print("cpg_annot_chrom")#debug
        #display(cpg_annot_chrom.head(3))#debug
        chrlen = int(CHR_LENGTH[CHR_LENGTH['Chromosome']==chrom].Total_length.values[0])
        #print("chrlen: {}".format(chrlen))#debug
        n_bins = (chrlen // binsize) + 1
        print("number of bins in {}: {}".format(chrom, n_bins))
        for i in range(n_bins):#for real #현재 chromosome의 총 길이를 10kb 단위로 sliding window하면서
        #for i in range(3):#debug for only the first 10 bins 
            if i % 500 == 0:
                print('---\nBin {}'.format(i+1))
            #print('---\n{}-th bin'.format(i+1))#debug
            current_window_start = i * binsize
            current_window_end = (i+1) * binsize

            #현재 보고 있는 window (genomic region) 내에 있는 cpg들을 골라내기 위한 boolean mask. 
            window_cpg_mask = [int(x)>=int(current_window_start) and int(x)<int(current_window_end) for x in cpg_annot_chrom.MAPINFO.values.flatten()] 

            #print("current window: from {} to {}".format(current_window_start, current_window_end))
            cur_cpg_annot_window = cpg_annot_chrom.iloc[window_cpg_mask, :]

            #print("cur_cpg_annot_window (CpGs in current window): {}".format(cur_cpg_annot_window.shape))
            #display(cur_cpg_annot_window.head(3))#debug

            # 현재 보고 있는 window 안의 cpg들의 목록 구하기
            intersecting_cpgs = np.intersect1d(cur_cpg_annot_window.Name.values.flatten(), beta.index.values.flatten())
            intersecting_cg_id_only = []
            for c in intersecting_cpgs:
                if c.startswith('cg'):
                    intersecting_cg_id_only.append(c)
            del(intersecting_cpgs)
            intersecting_cpgs = intersecting_cg_id_only.copy()
            del(intersecting_cg_id_only)

            #print("number of intersecting_cpgs (CpGs in current genomic region): {}".format(len(intersecting_cpgs)))#debug

            if len(intersecting_cpgs)>=10: # PMD가 되려면 최소 10개 이상의 methylated CpG dinucleotide를 가져야 함.
                beta_extracted = beta.loc[intersecting_cpgs]
                #print("Mean DNA methylation in current genomic region: {}".format(beta_extracted.mean().mean()))#debug
                if beta_extracted.mean().mean() < 0.7: #average methylation level이 0.7보다 낮아야 함.
                    pmd_df_data_chrom.append(chrom)
                    pmd_df_data_start.append(current_window_start)
                    pmd_df_data_end.append(current_window_end)
                    pmd_df_data_cpg.append(intersecting_cpgs)
                    chrom_pmd_cpg_list = np.concatenate((chrom_pmd_cpg_list, intersecting_cpgs), axis = None)
                    #print("---\nIdentified PMD")
                    #print("chrom: {}, current_window_start: {}, current_window_end: {}, intersecting_cpgs (last 5 elements): {}".format(chrom, current_window_start, current_window_end, intersecting_cpgs[-5:]))
                    #print("chrom_pmd_cpg_list (last 5 elements) : ", chrom_pmd_cpg_list[-5:])

                else:
                    non_pmd_df_data_chrom.append(chrom)
                    non_pmd_df_data_start.append(current_window_start)
                    non_pmd_df_data_end.append(current_window_end)
                    non_pmd_df_data_cpg.append(intersecting_cpgs)
                    chrom_non_pmd_cpg_list = np.concatenate((chrom_non_pmd_cpg_list, intersecting_cpgs), axis = None)
                    #print("---\nIdentified not-PMD")
                    #print("chrom: {}, current_window_start: {}, current_window_end: {}, intersecting_cpgs (last 5 elements): {}".format(chrom, current_window_start, current_window_end, intersecting_cpgs[-5:]))
                    #print("chrom_non_pmd_cpg_list (last 5 elements): ", current_non_pmd_cpg_list[-5:])

            else:
                non_pmd_df_data_chrom.append(chrom)
                non_pmd_df_data_start.append(current_window_start)
                non_pmd_df_data_end.append(current_window_end)
                non_pmd_df_data_cpg.append(intersecting_cpgs)
                chrom_non_pmd_cpg_list = np.concatenate((chrom_non_pmd_cpg_list, intersecting_cpgs), axis = None)
                #print("---\nIdentified not-PMD")
                #print("chrom: {}, current_window_start: {}, current_window_end: {}, intersecting_cpgs (last 5 elements): {}".format(chrom, current_window_start, current_window_end, intersecting_cpgs[-5:]))
                #print("chrom_non_pmd_cpg_list (last 5 elements): ", current_non_pmd_cpg_list[-5:])
        pmd_cpg_dictionary[chrom] = chrom_pmd_cpg_list
        non_pmd_cpg_dictionary[chrom] = chrom_non_pmd_cpg_list


    # make pmd_df
    pmd_df = pd.DataFrame(zip(pmd_df_data_chrom, pmd_df_data_start, pmd_df_data_end, pmd_df_data_cpg), columns = ['chrom', 'start', 'end', 'cpg'])

    # make non_pmd_df
    non_pmd_df = pd.DataFrame(zip(non_pmd_df_data_chrom, non_pmd_df_data_start, non_pmd_df_data_end, non_pmd_df_data_cpg), columns = ['chrom', 'start', 'end', 'cpg'])

    #print("---\npmd_df")#debug
    #display(pmd_df.head(3))#debug
    #print("---\nnon_pmd_df")#debug
    #display(non_pmd_df.head(3))#debug
    return pmd_df, non_pmd_df, pmd_cpg_dictionary, non_pmd_cpg_dictionary


# In[ ]:


if __name__=='__main__':
    
    #for cohort in COHORTS:
    args = parse_arguments()
    
    cohort_dir = os.path.join(SAVEDIR, args.cohort)  #현재 cohort의 결과가 이 디렉토리에 저장돼야 함. 
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    print("===\ncohort: {}".format(args.cohort))
    T, N, S = get_sample_list(args.cohort)
    
    print("Import beta value.")
    beta = pd.read_csv(os.path.join(METH_DIR, args.cohort+'.HumanMethylation450.tsv'), sep = '\t', index_col=0) #row: cpg, column: TCGA sample
    beta_tumor = beta[T].copy()
    beta_normal = beta[N].copy()
    
    print("Find PMD and non-PMD regions (in TUMOR samples only)")
    pmd_df, non_pmd_df, pmd_cpg_dictionary, non_pmd_cpg_dictionary = get_pmd(beta_tumor, args.binsize, CHR_LENGTH)
    
    print("Find PMD-like regions (in NORMAL samples only)")
    pmd_like_df, pmd_not_like_df, pmd_like_cpg_dictionary, pmd_not_like_cpg_dictionary = get_pmd(beta_normal, args.binsize, CHR_LENGTH)
    
    # make dataframe for bed file
    pmd_df_bed = pmd_df[['chrom', 'start', 'end']].copy()
    non_pmd_df_bed = non_pmd_df[['chrom', 'start', 'end']].copy()
    pmd_like_df_bed = pmd_like_df[['chrom', 'start', 'end']].copy()
    pmd_not_like_df_bed = pmd_not_like_df[['chrom', 'start', 'end']].copy()
        
    # save results
    pmd_df_fname = os.path.join(cohort_dir, 'tumor_pmd.csv')
    pmd_df_bed_fname = os.path.join(cohort_dir, 'tumor_pmd.bed')
    non_pmd_df_fname = os.path.join(cohort_dir, 'tumor_non_pmd.csv')
    non_pmd_df_bed_fname = os.path.join(cohort_dir, 'tumor_non_pmd.bed')
    pmd_cpg_dictionary_fname = os.path.join(cohort_dir, 'tumor_pmd_cpg_dictionary')#npz
    non_pmd_cpg_dictionary_fname = os.path.join(cohort_dir, 'tumor_non_pmd_cpg_dictionary')#npz
    pmd_like_df_fname = os.path.join(cohort_dir, 'normal_pmd-like.csv')
    pmd_like_df_bed_fname = os.path.join(cohort_dir, 'normal_pmd-like.bed')
    pmd_not_like_df_fname = os.path.join(cohort_dir, 'normal_pmd-not-like.csv')
    pmd_not_like_df_bed_fname = os.path.join(cohort_dir, 'normal_pmd-not-like.bed')
    pmd_like_cpg_dictionary_fname = os.path.join(cohort_dir, 'normal_pmd-like-cpg_dictipnary')#npz
    pmd_not_like_cpg_dictionary_fname = os.path.join(cohort_dir, 'normal_pmd-not-like-cpg_dictipnary')#npz
    
    
    pmd_df.to_csv(pmd_df_fname)
    non_pmd_df.to_csv(non_pmd_df_fname)
    np.savez(pmd_cpg_dictionary_fname, **pmd_cpg_dictionary)
    np.savez(non_pmd_cpg_dictionary_fname, **non_pmd_cpg_dictionary)
    pmd_like_df.to_csv(pmd_like_df_fname)
    pmd_not_like_df.to_csv(pmd_not_like_df_fname)
    np.savez(pmd_like_cpg_dictionary_fname, **pmd_like_cpg_dictionary)
    np.savez(pmd_not_like_cpg_dictionary_fname, **pmd_not_like_cpg_dictionary)
    pmd_df_bed.to_csv(pmd_df_bed_fname)
    non_pmd_df_bed.to_csv(non_pmd_df_bed_fname)
    pmd_like_df_bed.to_csv(pmd_like_df_bed_fname)
    pmd_not_like_df_bed.to_csv(pmd_not_like_df_bed_fname)
    
    print("result files: ")
    print("---Tumors: pmd, non-pmd")
    print('pmd_csv:', pmd_df_fname)
    print('pmd_bedfile:', pmd_df_bed_fname)
    print("non-pmd_csv:", non_pmd_df_fname)
    print("non-pmd_bedfile:", non_pmd_df_bed_fname)
    print("pmd_cpg_dictionary:", pmd_cpg_dictionary_fname+'.npz')
    print("non-pmd_cpg_dictionary:", non_pmd_cpg_dictionary_fname+'.npz')
    
    print("---\nNormals: pmd-like, pmd-not-like")
    print("pmd-like_csv:", pmd_like_df_fname)
    print("pmd-like_bedfile:",pmd_like_df_bed_fname)
    print("pmd-not-like_csv:",pmd_not_like_df_fname)
    print("pmd_not_like_bedfile:", pmd_not_like_df_bed_fname)
    print("pmd_like_cpg_dictionary:", pmd_like_cpg_dictionary_fname+'.npz')
    print("pmd_not_like_cpg_dictionary:",pmd_not_like_cpg_dictionary_fname+'.npz')

