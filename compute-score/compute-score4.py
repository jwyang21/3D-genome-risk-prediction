#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import argparse
from scipy.spatial.distance import euclidean


# In[51]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/')


# ### score4: reference로 pcbc 사용. 전반적인 계산 방법은 score2와 동일

# In[3]:


BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
PCBC_BINNED_DIFFMAT_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/pcbc-compute-binned-diffmat/result/'

TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] + ['chrX','chrY'] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/'
CHR_LENGTH_FNAME = '/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv'
CHR_LENGTH = pd.read_csv(CHR_LENGTH_FNAME)


# In[ ]:





# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    parser.add_argument('-b', '--binsize', type=int, default=1e6)
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    #parser.add_arguemnt('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_arguemnt('-t', '--threshold', help = 'Threshold for score1', type = float, required = True)
    #parser.add_arguemnt('-r', '--reference', help = 'Reference for score2', type = str, required = True)
    #parser.add_arguemnt('-d', '--distance', help = 'Distance metric for score2', type = str, required = True)
    #parser.add_argument('-cr', '--chrom', help = 'chromosome', type = str, required = True)
    return parser.parse_args()


# In[5]:


def get_sample_list(directory):
    file_list = os.listdir(directory)
    T = []
    N = []
    for f in file_list:
        if 'TCGA' in f and f.endswith('.npz'):
            if f[13:15] in TUMOR_BARCODE:
                T.append(f.split('.')[0])
            elif f[13:15] in NORMAL_BARCODE:
                N.append(f.split('.')[0])
            else:
                raise Exception("Wrong TCGA barcode! This sample does not belong to neither tumor nor normal.")
    S = T + N
    print("len(tumor), len(normal), len(all): {}, {} and {}, respectively.".format(len(T), len(N), len(S)))
    return T, N, S


# In[35]:


def import_pcbc_bins(chrom):
    
    return np.load(os.path.join(PCBC_BINNED_DIFFMAT_DIR, 'H9_A', chrom+'_binned_diffmat.npz'), allow_pickle = True)['bins']


# In[36]:


def import_pc1_ref(chrom, bin_index):
    pc1 = np.load(os.path.join(PC1_DIR, 'pcbc', 'ref_all.npz'), allow_pickle = True)[chrom]
    return pd.DataFrame(pc1, index = bin_index, columns = ['pc1'])


# In[37]:


def import_450k_bins(chrom, binsize, cohort, sample):
    n_bin = int(int(CHR_LENGTH[CHR_LENGTH['Chromosome']==chrom].Total_length) // binsize + 1)
    bin_index = [chrom+':'+str(int(i * binsize))+'-'+str(int((i+1) * binsize)) for i in range(n_bin)]
    diffmat_mask = np.load(os.path.join(BINNED_DIFFMAT_DIR, cohort, sample+'.npz'), allow_pickle = True)[chrom+'_mask']
    
    return np.array(bin_index)[~diffmat_mask].tolist()


# In[38]:


def import_pc1_450k(cohort, sample, chrom, bin_index):
    f = os.path.join(PC1_DIR, cohort, sample+'_inv_exp.npz')
    pc1 = np.load(f, allow_pickle = True)[chrom]
    
    return pd.DataFrame(pc1, index = bin_index, columns = ['pc1'])







# In[ ]:





# In[ ]:


if __name__ == '__main__':   
    args = parse_arguments()
    
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    df = pd.DataFrame(np.zeros((len(S), len(CHR_LIST)), dtype = float), columns = CHR_LIST, index = S)
    
    for chrom in CHR_LIST:
        pcbc_bins = import_pcbc_bins(chrom)
        pc1_ref = import_pc1_ref(chrom, pcbc_bins) #pcbc ref pc1를 import #flattened vector #df (index: bin names, column: pc1)
        
        for s in S:
            if S.index(s) == 0: # import 450k pc1 bin names if this is the first sample
                pc1_450k_bins = import_450k_bins(chrom, args.binsize, args.cohort, s) 
                
            pc1_450k = import_pc1_450k(args.cohort, s, chrom, pc1_450k_bins) # pc1 of current chromosome #flattened vector #df (index: bin names, column: pc1)
            
            # extract intersecting bins only from pc1_450k and pc1_ref
            pc1_450k_ = pc1_450k.loc[np.intersect1d(pc1_450k.index.values, pc1_ref.index.values)]
            pc1_ref_ = pc1_ref.loc[np.intersect1d(pc1_450k.index.values, pc1_ref.index.values)]
            
            df.loc[s][chrom] = euclidean(pc1_450k_.pc1.values.flatten(), pc1_ref_.pc1.values.flatten())
            
    df['simple_avg'] = df.mean(axis = 1).values
    df.to_pickle(os.path.join(SAVEDIR, 'score4_euclidean.pickle'))
    print(os.path.join(SAVEDIR, 'score4_euclidean.pickle')) # column이 chromosome, row가 해당 cohort의 sample names


# In[ ]:




