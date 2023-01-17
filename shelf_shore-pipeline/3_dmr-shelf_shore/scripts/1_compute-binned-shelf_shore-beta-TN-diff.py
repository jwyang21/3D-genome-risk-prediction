#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import argparse


# In[5]:


# # pipeline
# key가 genomic binned regionn이고 value가 거기 위치한 cpg probe list인 딕셔너리 제작
# 
# 각 cohort마다:
#     - 각 sample마다;
#         - 각 chromosome마다:
#             - 각 bin마다:
#                 - 그 bin에 속한 cpg들의 beta 값들만 모아 평균내기
# 각 cohort마다:
#     - 각 chormosome마다:
#         - tumor들의 avg. beta 값을 평균냄
#         - normal들의 avg. beta 값을 평균냄
#         - (T-N)의 값 >= threshold인 genomic region들 찾아서 저장. 
#         
# 위에서 찾은 genomic region에 대해:
#     - cgc gene의 region과 겹치는 거 있는지 찾기
#     - 이 region에 위치한 gene들 찾기. 


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[8]:


#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
#BETA_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA-LIHC_hg19/TCGA.LIHC.sampleMap%2FHumanMethylation450'
#BETA = pd.read_csv(BETA_FNAME, sep = '\t', index_col = 0)

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
#SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
#CPG_ANNOT = pd.read_csv('/data/project/jeewon/research/3D-ITH/data/illumina/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
CHR_LENGTH.columns = ['chrom', 'len']

# refer to /data/project/jeewon/research/3D-ITH/binned_diff/snake/*probe*


#command: python3 1_compute-binned-island-beta-TN-diff.py -b $binsize -w_dir /data/project/3dith/pipelines/island-pipeline/3_dmr-island -c $cohort --cpg_metadata_fname /data/project/3dith/data/450k_metadata.island.sorted.bed --cpg_type island --dmr_type TN
# In[9]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--binsize', type=int, default=int(1e6), required = False)
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('-c', '--cohort', type = str, required = True)
    parser.add_argument('--cpg_metadata_fname', help = 'CpG metadata filename', type = str, required = True)#'/data/project/3dith/data/450k_metadata.{CPG_TYPE}.sorted.bed'
    parser.add_argument('--cpg_type', help = 'CpG type. opensea/island/shelf/shore/shelf_shore', type = str, required = True)
    parser.add_argument('--dmr_type', help = 'DMR type. TN (tumor-normal) or HL (high score - low score)', type = str, required = True)
    
   
    return parser.parse_args()


# In[10]:


def save_dictionary_to_csv_files(savedir, d, binsize, prefix=''):#use this in case when values in dictionary are dataframes.
    #ref: https://stackoverflow.com/questions/50786266/writing-dictionary-of-dataframes-to-file
    for k, v in d.items():
        print("key:{}".format(k))
        fname = os.path.join(savedir, prefix+str(k)+'_binsize_'+str(binsize)+'.csv')
        v.to_csv(fname)
        #print(fname)

        
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


def make_bin_target_type_cpg_dict(binsize):
    bin2cpg = {}
    
    for chrom in CHR_LIST:#real
    #for chrom in CHR_LIST[18:19]:#for test
        #current_chrom_df = CPG_METADATA[CPG_METADATA['chrom']==chrom].copy()
        chrom_mask = np.array([str(x)==chrom for x in CPG_METADATA['chrom'].values.flatten()])
        current_chrom_df = CPG_METADATA.iloc[chrom_mask,:].copy()
        
        for i in range(current_chrom_df.shape[0]):#real
        #for i in range(10):#test
            current_cpg = current_chrom_df.cpg.values[i]
            #print(current_cpg)
            current_cpg_mapinfo = int(current_chrom_df.end.values[i]) - 1
            #print(current_cpg_mapinfo)
            current_chrom_len = int(CHR_LENGTH[CHR_LENGTH['chrom']==chrom].len.values[0])
            current_chrom_num_bin = int(current_chrom_len // binsize)+1
            for j in range(current_chrom_num_bin):
                current_bin_start = int(j) * binsize
                current_bin_end = int((j+1) * binsize)
                
                if current_cpg_mapinfo >= current_bin_start and current_cpg_mapinfo < current_bin_end:
                    # use 0-based bin name since it should be compatible with gene, chrom, and regulatory features annotation file.
                    # (aforementioned files are all 0-based)
                    # for 0-based and 1-based system, refer to https://www.biostars.org/p/84686/
                    
                    #current_bin_name = chrom+':'+str(current_bin_start + 1)+'-'+str(current_bin_end)#former version. 1-based
                    current_bin_name = chrom+':'+str(current_bin_start)+'-'+str(current_bin_end)#current version. 0-based
                    #print(current_bin_start, current_bin_end, current_bin_name)
                    if current_bin_name not in list(bin2cpg.keys()):
                        bin2cpg[current_bin_name] = []
                    if current_cpg not in bin2cpg[current_bin_name]:
                        bin2cpg[current_bin_name].append(current_cpg)      
    return bin2cpg


# In[14]:


def compute_binned_avg_target_type_beta_df(S, bin2cpg, savedir, cpg_type):
    '''
    - 각 chromosome마다:
        - 각 bin마다:
            - 그 bin에 속한 cpg들의 beta 값들만 모아 평균내기
    '''
    all_chrom_results = {}

    for chrom in CHR_LIST:#real
    #for chrom in CHR_LIST[18:19]:#for test
        beta_target_type_df = pd.read_csv(os.path.join(savedir, 'beta_'+cpg_type+'_df.csv'), index_col = 0)
        current_chrom_target_type_bins = []
        for i in range(len(list(bin2cpg.keys()))):
            current_bin = list(bin2cpg.keys())[i]
            #if chrom in current_bin: #이렇게 하면 chr1의 DMR bin에 chr18:xxx-yyy 또는 chr19:xxx-yyy bin들도 포함됨.
            current_bin_chrom = str(current_bin.split(':')[0].strip())
            if current_bin_chrom == chrom:
                current_chrom_target_type_bins.append(current_bin)
        current_chrom_df = pd.DataFrame(np.zeros((len(current_chrom_target_type_bins), len(S)), dtype = float), index = current_chrom_target_type_bins, columns = S)
        for k in range(len(current_chrom_target_type_bins)):
            current_bin = current_chrom_target_type_bins[k]
            current_bin_cpg = bin2cpg[current_bin]
            intersecting_cpgs = np.intersect1d(beta_target_type_df.index.values, current_bin_cpg)
            current_bin_cpg_betas  = beta_target_type_df.loc[intersecting_cpgs]
            values = current_bin_cpg_betas.mean().values
            current_chrom_df.loc[current_bin] = values

        all_chrom_results[chrom] = current_chrom_df
            
    return all_chrom_results


def compute_binned_avg_target_type_BetaDiff(T, N, binned_avg_target_type_beta_dfs):
    '''
    각 chormosome마다:
        - tumor들의 avg. beta 값을 평균냄
        - normal들의 avg. beta 값을 평균냄
        - (T-N)의 값 >= threshold인 genomic region들 찾아서 저장. 
    '''
    binned_avg_target_type_BetaDiff_dfs = {}
    
    for chrom in CHR_LIST:#real
    #for chrom in CHR_LIST[18:19]:#test
        current_df = binned_avg_target_type_beta_dfs[chrom].copy()
        
        current_df_T = current_df[T].copy() #tumor samples only
        current_df_N = current_df[N].copy()
        assert current_df_T.shape[1] == len(T)
        assert current_df_N.shape[1] == len(N)
        
        current_df_T_rowmean = current_df_T.mean(axis=1).values
        current_df_N_rowmean = current_df_N.mean(axis=1).values
        
        current_df_T_N = (current_df_T_rowmean - current_df_N_rowmean).copy()
        binned_avg_target_type_BetaDiff_dfs[chrom] = current_df_T_N.flatten()
        binned_avg_target_type_BetaDiff_dfs[chrom+'_bins'] = current_df_T.index.values
        
    return binned_avg_target_type_BetaDiff_dfs


if __name__=='__main__':
    
    args = parse_arguments()
    os.chdir(args.working_dir) 
    
    CPG_METADATA = pd.read_csv(args.cpg_metadata_fname, header=None, sep = '\t')
    CPG_METADATA.columns = ['chrom', 'start', 'end', 'cpg']
    
    RESULTDIR = os.path.join(os.getcwd(), 'result')
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    print("cohort: {}".format(args.cohort))
    print("SAVEDIR: {}".format(SAVEDIR))
    
    T, N, S = get_sample_list(args.cohort)
    
    print("---\n1. Make bin_to_{}-cpg_dictionary".format(args.cpg_type))
    #key가 genomic binned regionn이고 value가 거기 위치한 cpg probe list인 딕셔너리 제작
    if not os.path.exists(os.path.join(RESULTDIR, 'bin_2_'+args.cpg_type+'_cpg_binsize_'+str(args.binsize)+'.npz')):
        bin2cpg = make_bin_target_type_cpg_dict(args.binsize)
        
        #save
        np.savez(os.path.join(RESULTDIR, 'bin_2_'+args.cpg_type+'_cpg_binsize_'+str(args.binsize)), **bin2cpg)
        print(os.path.join(RESULTDIR, 'bin_2_'+args.cpg_type+'_cpg_binsize_'+str(args.binsize)+'.npz'))
    else:
        bin2cpg = np.load(os.path.join(RESULTDIR, 'bin_2_'+args.cpg_type+'_cpg_binsize_'+str(args.binsize)+'.npz'))

    print("---\n2 Binned avg {} beta value".format(args.cpg_type))
    binned_avg_target_type_beta_dfs = compute_binned_avg_target_type_beta_df(S, bin2cpg, SAVEDIR, args.cpg_type) #dictionary #key: ['chr1', 'chr2', ...]
    save_dictionary_to_csv_files(SAVEDIR, binned_avg_target_type_beta_dfs, args.binsize, 'binned_avg_'+args.cpg_type+'_beta_')

    print("---\n3. Binned avg {} beta value diffmat".format(args.cpg_type))
    binned_avg_target_type_BetaDiff = compute_binned_avg_target_type_BetaDiff(T, N, binned_avg_target_type_beta_dfs)
    
    npz_fname = 'binned_avg_'+args.cpg_type+'_'+args.dmr_type+'_diff_binsize_'+str(args.binsize)
    np.savez(os.path.join(SAVEDIR, npz_fname), **binned_avg_target_type_BetaDiff)
    # filename example) binned_avg_opensea_TN_diff_binsize_1000000.npz
    print(os.path.join(SAVEDIR, npz_fname+'.npz'))
