#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
import argparse
import random


# In[ ]:


random.seed(2022)
np.random.seed(2022)


# In[ ]:


ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
NORMAL7_COHORT_PCBC = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD PCBC'.split(' ')


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


# In[ ]:


def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', type = str, required = True)

    args.add_argument('--cpg_type', help = 'opensea, island, shelf, shore, shelf_shore', type = str, required = True)
    
    args.add_argument('--tcga_bdm_dir', help = 'TCGA binned difference matrix directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/binned-difference-matrix-v2-{cpg_type}/result' #이 디렉토리 내의 cohort 디렉토리에 IEBDM들 있음.

    args.add_argument('--pcbc_bdm_dir', help = 'PCBC binned difference matrix directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/binned-difference-matrix-pcbc/result' #IEBDM들이 저장된 디렉토리
    
    args.add_argument('--hg19_chr_len', help = 'hg19 chromosome length file', type = str, required = True)
    #default: /data/project/3dith/data/hg19.fa.sizes'
    
    args.add_argument('-b', '--binsize', type = int, required = True, default = int(1e6))

    args.add_argument('--save_dir', type = str, required = True)
    #/data/project/3dith/data/bdm_bins
    
    return args.parse_args()


# In[ ]:


def import_bdm(directory, chrom, s): #directory should be data_dir
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T

    return diffmat_masked


# In[ ]:


# In[12]:
if __name__=='__main__':
    
    args = parse_arguments()
    os.chdir(args.working_dir)
    print("===cpg_type: ", args.cpg_type)

    if not os.path.exists(args.save_dir):
        os.makedirs(args.save_dir)
    if not os.path.exists(os.path.join(args.save_dir, args.cpg_type)):
        os.makedirs(os.path.join(args.save_dir, args.cpg_type))
    
    chr_len = pd.read_csv(args.hg19_chr_len, sep = '\t', header = None, index_col = 0, names = ['len'])

    for cohort in NORMAL7_COHORT_PCBC:
        
        globals()[cohort+'_diffmat_bins'] = {}
        
        T, N, S = get_sample_list(cohort) 
        
        if 'TCGA' in cohort:
            bdm_dir = args.tcga_bdm_dir
            data_dir = os.path.join(bdm_dir, cohort)
        else:#PCBC
            bdm_dir = args.pcbc_bdm_dir
            data_dir = bdm_dir
            
        for chrom in CHR_LIST:
            bdm = import_bdm(data_dir, chrom, S[0])
            n_bin = int(chr_len.loc[chrom].values // args.binsize)+1
            bin_index = [chrom+':'+str(int(i * args.binsize))+'-'+str(int((i+1) * args.binsize)) for i in range(n_bin)]
            diffmat_mask = np.load(os.path.join(data_dir, S[0]+'.npz'))[chrom+'_mask']
            globals()[cohort+'_diffmat_bins'][chrom+'_bins'] = np.array(bin_index)[~diffmat_mask].tolist()
            
        result_fname = cohort+'_diffmat_bins'#npz
        
        np.savez(os.path.join(args.save_dir, args.cpg_type, result_fname), **globals()[cohort+'_diffmat_bins'])
        
        print("{}: {}".format(cohort, os.path.join(args.save_dir, args.cpg_type, result_fname)))        