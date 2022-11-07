#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# To-Do: compute pc1 from original binned difference matrix (of each chromosome)

import numpy as np
import pandas as pd
import os
import sys
from sklearn.decomposition import PCA
import argparse


# diffmat의 PC1, inverse exponential diffmat의 PC1 계산

# In[ ]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/')


# In[ ]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] + ['chrX','chrY'] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
#FIRE_DIR = '/data/project/3dith/data/'
#FIRE_AB_FNAME = 'fire_a_b_labels.csv'
#FIRE_PC1_FNAME = 'fire_pc1.csv'
#COHORT_TCGA_FIRE_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA_FIRE_intersection.csv'
#COHORT_TCGA_FIRE = pd.read_csv(COHORT_TCGA_FIRE_FNAME)
BETA_PMD_DIFFMAT_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/beta-pmd-binned-diffmat/result/'#TCGA-ACC/TCGA-OR-A5KW-01_chr1.npz # value, colnames, rownames
PMD_DIFFMAT_ITEMS = ['value', 'colnames', 'rownames']


# In[ ]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=1e6)
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    #parser.add_arguemnt('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA_cohort', type = str, required = True)
    #parser.add_arguemnt('-t', '--threshold', help = 'Threshold for score1', type = float, required = True)
    #parser.add_arguemnt('-r', '--reference', help = 'Reference for score2', type = str, required = True)
    #parser.add_arguemnt('-d', '--distance', help = 'Distance metric for score2', type = str, required = True)
    parser.add_argument('-cr', '--chrom', help = 'chromosome(chr1~chr22)', type = str, required = True)
    return parser.parse_args()


# In[ ]:


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


# In[ ]:


'''
def import_binned_diffmat(directory, chrom): #directory should be DATA_DIR
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T
    #if args.score == 'score1':
        #return diffmat, diffmat_mask, diffmat_masked
    #else: #score2 or score3
        #return diffmat, diffmat_mask, 1/np.exp(diffmat_masked)
    return diffmat_masked
'''


# In[ ]:


def import_pmd_binned_diffmat(cohort, sample, chrom): #chrom should be in format of 'chr1', 'chr2', ...
    
    npz_fname = sample+'_'+chrom+'.npz'
    f = os.path.join(BETA_PMD_DIFFMAT_DIR, cohort, npz_fname)
    value = np.load(f, allow_pickle=True)['value']
    colnames = np.load(f, allow_pickle=True)['colnames']
    rownames = np.load(f, allow_pickle=True)['rownames']
    df = pd.DataFrame(value, index = rownames, columns = colnames) 
    
    return df


# In[ ]:


def pc1(m):
    
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    pc1 = pc[:,0]

    return pc1


# In[ ]:





# In[ ]:




# In[ ]:


# BETA_PMD_DIFFMAT으로부터 PC1 계산. 

if __name__ == '__main__':
    args = parse_arguments()
    
    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, 'pmd-diffmat')
    DATA_DIR = os.path.join(BINNED_DIFFMAT_DIR, args.cohort) # binned diffmat dir of current cohort #used to get list of all sample names. 

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(os.path.join(os.getcwd(), 'result', args.cohort)):
        os.makedirs(os.path.join(os.getcwd(), 'result', args.cohort))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
        
    print("SAVEDIR: {}".format(SAVEDIR))

    T, N, S = get_sample_list(DATA_DIR) 

    for s in S:#real implementation
    #for s in S[:2]: #for test
        if S.index(s) % 100 == 0:
            print("s: {}".format(s))
        print(args.chrom)
        m = import_pmd_binned_diffmat(args.cohort, s, args.chrom)
        globals()['{}_pc1'.format(args.chrom)] = pc1(m)
        globals()['{}_inv_exp_pc1'.format(args.chrom)] = pc1(1/np.exp(m))
        
        
        # save temporary result. (only for chr1, of which the calculation is completed.)    
        ## run for chromosome 1 only. 
        npz_fname1 = os.path.join(SAVEDIR, s) #diffmat의 PC1
        npz_fname2 = os.path.join(SAVEDIR, s + '_inv_exp') #diffmat의 PC1
        
        print("npz_fname: {}".format(npz_fname1+'.npz'))
        np.savez(npz_fname1, chr1 = chr1_pc1) #chrX = chrX_pc1, chrY = chrY_pc1
        
        npz_fname2 = os.path.join(SAVEDIR, s + '_inv_exp') # ivnerse diffmat의 PC1
        print("npz_fname: {}".format(npz_fname2 + '.npz'))
        np.savez(npz_fname2, chr1 = chr1_inv_exp_pc1)
    
        '''
        # save result of all chromosomes
        print("npz_fname: {}".format(npz_fname1+'.npz'))
                np.savez(npz_fname1, chr1 = chr1_pc1, chr2 = chr2_pc1, chr3 = chr3_pc1, chr4 = chr4_pc1, chr5 = chr5_pc1, chr6 = chr6_pc1,
                chr7 = chr7_pc1, chr8 = chr8_pc1, chr9 = chr9_pc1, chr10 = chr10_pc1, chr11 = chr11_pc1, chr12 = chr12_pc1, chr13 = chr13_pc1,
                chr14 = chr14_pc1, chr15 = chr15_pc1, chr16 = chr16_pc1, chr17 = chr17_pc1, chr18 = chr18_pc1, chr19 = chr19_pc1,
                    chr20 = chr20_pc1, chr21 = chr21_pc1, chr22 = chr22_pc1) #chrX = chrX_pc1, chrY = chrY_pc1
        
        npz_fname2 = os.path.join(SAVEDIR, s + '_inv_exp') # ivnerse diffmat의 PC1
        print("npz_fname: {}".format(npz_fname2 + '.npz'))
                np.savez(npz_fname2, chr1 = chr1_inv_exp_pc1, chr2 = chr2_inv_exp_pc1, chr3 = chr3_inv_exp_pc1, chr4 = chr4_inv_exp_pc1, chr5 = chr5_inv_exp_pc1, chr6 = chr6_inv_exp_pc1,
                chr7 = chr7_inv_exp_pc1, chr8 = chr8_inv_exp_pc1, chr9 = chr9_inv_exp_pc1, chr10 = chr10_inv_exp_pc1, chr11 = chr11_inv_exp_pc1, chr12 = chr12_inv_exp_pc1, 
                    chr13 = chr13_inv_exp_pc1, chr14 = chr14_inv_exp_pc1, chr15 = chr15_inv_exp_pc1, chr16 = chr16_inv_exp_pc1, chr17 = chr17_inv_exp_pc1, 
                    chr18 = chr18_inv_exp_pc1, chr19 = chr19_inv_exp_pc1, chr20 = chr20_inv_exp_pc1, chr21 = chr21_inv_exp_pc1, chr22 = chr22_inv_exp_pc1)
        '''
    


# In[ ]:


# ORIGINAL BINNED DIFFMAT으로부터 PC1 추출. (version 1 code)
## ORIGINAL BINNED DIFFMAT filenames: /data/project/3dith/pipelines/binned-difference-matrix-v2/result/[TCGA-COHORT]/[BARCODE].npz

'''
if __name__ == '__main__':
    args = parse_arguments()
    
    print("cohort: {}".format(args.cohort))
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    DATA_DIR = os.path.join(BINNED_DIFFMAT_DIR, args.cohort) # binned diffmat dir of current cohort

    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    T, N, S = get_sample_list(DATA_DIR) 
    
    for s in S:
        if S.index(s) % 100 == 0:
            print("s: {}".format(s))
        for chrom in CHR_LIST:
            print(chrom)
            m = import_binned_diffmat(DATA_DIR, chrom)
            globals()['{}_pc1'.format(chrom)] = pc1(m)
            globals()['{}_inv_exp_pc1'.format(chrom)] = pc1(1/np.exp(m))
        
        npz_fname1 = os.path.join(SAVEDIR, s) #diffmat의 PC1
        print("npz_fname: {}".format(npz_fname1+'.npz'))
                np.savez(npz_fname1, chr1 = chr1_pc1, chr2 = chr2_pc1, chr3 = chr3_pc1, chr4 = chr4_pc1, chr5 = chr5_pc1, chr6 = chr6_pc1,
                chr7 = chr7_pc1, chr8 = chr8_pc1, chr9 = chr9_pc1, chr10 = chr10_pc1, chr11 = chr11_pc1, chr12 = chr12_pc1, chr13 = chr13_pc1,
                chr14 = chr14_pc1, chr15 = chr15_pc1, chr16 = chr16_pc1, chr17 = chr17_pc1, chr18 = chr18_pc1, chr19 = chr19_pc1,
                 chr20 = chr20_pc1, chr21 = chr21_pc1, chr22 = chr22_pc1) #chrX = chrX_pc1, chrY = chrY_pc1
        
        npz_fname2 = os.path.join(SAVEDIR, s + '_inv_exp') #diffmat의 PC1
        print("npz_fname: {}".format(npz_fname2 + '.npz'))
                np.savez(npz_fname2, chr1 = chr1_inv_exp_pc1, chr2 = chr2_inv_exp_pc1, chr3 = chr3_inv_exp_pc1, chr4 = chr4_inv_exp_pc1, chr5 = chr5_inv_exp_pc1, chr6 = chr6_inv_exp_pc1,
                chr7 = chr7_inv_exp_pc1, chr8 = chr8_inv_exp_pc1, chr9 = chr9_inv_exp_pc1, chr10 = chr10_inv_exp_pc1, chr11 = chr11_inv_exp_pc1, chr12 = chr12_inv_exp_pc1, 
                 chr13 = chr13_inv_exp_pc1, chr14 = chr14_inv_exp_pc1, chr15 = chr15_inv_exp_pc1, chr16 = chr16_inv_exp_pc1, chr17 = chr17_inv_exp_pc1, 
                 chr18 = chr18_inv_exp_pc1, chr19 = chr19_inv_exp_pc1, chr20 = chr20_inv_exp_pc1, chr21 = chr21_inv_exp_pc1, chr22 = chr22_inv_exp_pc1)
                #chrX = chrX_inv_exp_pc1, chrY = chrY_inv_exp_pc1
'''
