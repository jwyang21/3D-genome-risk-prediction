#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import sys

DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
NORMAL_BARCODES = ['11', '12', '13','14', '15', '16', '17', '18', '19']
#Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{TCGA-cohort}/{sample}.npz or #/{TCGA-cohort}/{sample}_inv_exp.npz #'chr1'
CHROM_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/hg19.fa.sizes', sep = '\t', header=None)
CHROM_LENGTH.columns = ['chr', 'len']
FIRE_PC1 = pd.read_csv('/data/project/3dith/data/fire_pc1.csv') #FIRE_PC1 'end' column contains NA, so preprocessing is needed.
FIRE_BINSIZE = int(1e6)
FIRE_PC1['end'] = [FIRE_PC1['start'].values[i]+FIRE_BINSIZE-1 for i in range(FIRE_PC1.shape[0])]
FIRE_AB = pd.read_csv('/data/project/3dith/data/fire_a_b_labels.csv')
FIRE2TCGA = pd.read_csv(os.path.join(DATA_DIR, 'tcga-fire-cohorts.csv'))
PROBEMAP_CGC = pd.read_csv(os.path.join(DATA_DIR, 'gene-expr-hg19-xena', 'probemap-cgc-sampled.csv'))
PROBEMAP_CGC.index = [PROBEMAP_CGC['chrom'].values[i]+':'+str(PROBEMAP_CGC['chromStart_binned'].values[i])+'-'+str(PROBEMAP_CGC['chromEnd_binned'].values[i]-1) for i in range(PROBEMAP_CGC.shape[0])]
SAVEDIR = os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/ab-agreement/result')


# In[4]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True)
    parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)
    parser.add_argument('-f', '--flag', help = 'binned diffmat type(raw, inv)', required=True)
    return parser.parse_args()


# In[5]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    cohort_binned_diffmat_dir = os.path.join(BINNED_DIFFMAT_DIR, cohort)
    T = []
    N = []
    S = [] #all samples
    for l in os.listdir(cohort_binned_diffmat_dir):
        if l.startswith('TCGA') and l.endswith('.npz'):
            if l[13:15] in TUMOR_BARCODES:
                T.append(l.split('.')[0].strip())
            elif l[13:15] in NORMAL_BARCODES:
                N.append(l.split('.')[0].strip())
            else:
                pass
        else:
            pass
    S = T + N
    print("{}: tumor {}, normal {}, total {}".format(cohort, len(T), len(N), len(S)))
    
    return T, N, S


# In[6]:


def import_binned_diffmat(cohort, sample, chrom, binsize):
    # return a binned diffmat
    fname = os.path.join(BINNED_DIFFMAT_DIR, cohort, sample+'.npz')
    raw_diffmat = np.load(fname)[chrom]
    raw_mask = np.load(fname)['{}_mask'.format(chrom)]
    diffmat_masked = raw_diffmat[~raw_mask].T[~raw_mask].T
    
    # return bins with mask applied.
    chrom_len = int(CHROM_LENGTH[CHROM_LENGTH['chr']==chrom].len.values[0])
    n_bins = int((chrom_len//binsize) + 1)
    bins = np.array([chrom+':'+str((i * binsize) + 1)+'-'+str((i+1)*binsize) for i in range(n_bins)])
    #return diffmat_masked, 1/np.exp(diffmat_masked)
    return diffmat_masked, 1/np.exp(diffmat_masked), bins[~raw_mask]


# In[7]:


def import_pc1(cohort, sample, chrom, flag):
    # import pre-computed PC1 of sample-of-interest
    if flag=='raw':
        fname = os.path.join(PC1_DIR, cohort, sample+'.npz')
    elif flag=='inv': #inverse exponential
        fname = os.path.join(PC1_DIR, cohort, sample+'_inv_exp.npz')
    pc1 = np.load(fname)[chrom]

    return pc1.flatten()


# In[8]:


def import_fire(tcga_cohort, flag, chrom, fire_cohort):
    # import FIRE A,B or PC1 data
    cols = ['chr', 'start', 'end']
    
    cols.append(fire_cohort)

    if flag=='ab':
        df = FIRE_AB[cols].dropna()
        df = df[df['chr']==int(chrom.split('chr')[-1])]
        df.index = ['chr'+str(df['chr'].values[i])+':'+str(df['start'].values[i])+'-'+str(df['end'].values[i]) for i in range(df.shape[0])]  
        
    elif flag=='pc1':
        df = FIRE_PC1[cols].dropna()
        df = df[df['chr']==int(chrom.split('chr')[-1])]
        df.index = ['chr'+str(df['chr'].values[i])+':'+str(df['start'].values[i])+'-'+str(df['end'].values[i]) for i in range(df.shape[0])]  
    
    else:
        pass

    return df


