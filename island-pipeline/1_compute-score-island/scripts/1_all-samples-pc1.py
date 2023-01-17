#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import sys
from sklearn.decomposition import PCA
import argparse

'''
# cohort는 TCGA, reference type은 PCBC -> stem distance of each TCGA sample
# cohort는 TCGA, reference type도 TCGA -> normal distance of each TCGA sample
# cohort는 PCBC, reference type도 PCBC -> stem distance of each PCBC sample
# cohort는 PCBC, reference type은 TCGA -> normal distance of each PCBC sample.
# distance 계산할 때 simple_avg하거나 weighted_avg (weighted by chromosome length / all autosome length)하거나. 
'''

# To Do: 각 sample의 각 chromosome의 inverse exponential diffmat의 PC1 계산

# GLOBAL VARIABLES
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='island'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    #default: 'data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}'

    parser.add_argument('--tcga_bdm_dir', help = 'TCGA binned difference matrix directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/binned-difference-matrix-v2-{CPG_TYPE}/result' #이 디렉토리 내의 cohort 디렉토리에 IEBDM들 있음.

    parser.add_argument('--pcbc_bdm_dir', help = 'PCBC binned difference matrix directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/binned-difference-matrix-pcbc/result' #IEBDM들이 저장된 디렉토리

    parser.add_argument('-r_dir', '--result_dir', help = 'result directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result' 

    parser.add_argument('-m_type', '--matrix_type', help = 'from which matrix you will compute PC1. bdm, iebdm, or all.', type = str, required = True)
    #default: 'iebdm'
    
    parser.add_argument('--cohort', help = 'cohort. TCGA-{} or PCBC', type = str, required = True)
    #usage: TCGA-xxxx or PCBC
    return parser.parse_args()

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

def import_bdm(directory, chrom, s): #directory should be data_dir
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T

    return diffmat_masked

def pc1(m):
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    #print(pc)
    pc1 = pc[:,0]
    #print(pc1)
    #print('-----')
    return pc1

if __name__ == '__main__':
    args = parse_arguments()
    
    os.chdir(args.working_dir)

    if 'TCGA' in args.cohort:
        bdm_dir = args.tcga_bdm_dir
        data_dir = os.path.join(bdm_dir, args.cohort) #binned diffmat dir of current cohort
    else:
        bdm_dir = args.pcbc_bdm_dir
        data_dir = bdm_dir #binned diffmat dir of current cohort

    print("===\ncohort: {}".format(args.cohort))
    
    cohort_dir = os.path.join(args.result_dir, args.cohort) #cohort result directory #result/{CPG_TYPE}/TCGA-BRCA

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)

    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)

    if not os.path.exists(os.path.join(cohort_dir, 'pc1')):
        os.makedirs(os.path.join(cohort_dir, 'pc1'))
    
    T, N, S = get_sample_list(args.cohort) 
    
    for s in S:
        print("s: {}".format(s))
        
        if args.matrix_type == 'bdm':
            bdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = import_bdm(data_dir, chrom, s)
                bdm_pc1[current_key] = pc1(m)
            # 현재 샘플의 22개 chromosome 각각에 대한 PC1 벡터들을 한 개의 npz file로 저장. 
            npz_fname1 = os.path.join(cohort_dir,'pc1', s) #bdm pc1 #cohort_dir: cohort result directory. 
            print("npz_fname: {}".format(npz_fname1+'.npz'))
            np.savez(npz_fname1, **bdm_pc1)
        
        elif args.matrix_type == 'iebdm':
            iebdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = 1/np.exp(import_bdm(data_dir, chrom, s))
                iebdm_pc1[current_key] = pc1(m)
            # 현재 샘플의 22개 chromosome 각각에 대한 PC1 벡터들을 한 개의 npz file로 저장.
            npz_fname2 = os.path.join(cohort_dir, 'pc1', s + '_inv_exp') #iebdm PC1
            print("npz_fname: {}".format(npz_fname2 + '.npz'))
            np.savez(npz_fname2, **iebdm_pc1)
        
        elif args.matrix_type == 'all':
            bdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = import_bdm(data_dir, chrom, s)
                bdm_pc1[current_key] = pc1(m)
            # 현재 샘플의 22개 chromosome 각각에 대한 PC1 벡터들을 한 개의 npz file로 저장. 
            npz_fname1 = os.path.join(cohort_dir, 'pc1', s) #bdm pc1 #cohort_dir: cohort result directory. 
            print("npz_fname: {}".format(npz_fname1+'.npz'))
            np.savez(npz_fname1, **bdm_pc1)

            iebdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = 1/np.exp(import_bdm(data_dir, chrom, s))
                iebdm_pc1[current_key] = pc1(m)
            # 현재 샘플의 22개 chromosome 각각에 대한 PC1 벡터들을 한 개의 npz file로 저장.
            npz_fname2 = os.path.join(cohort_dir, 'pc1', s + '_inv_exp') #iebdm PC1
            print("npz_fname: {}".format(npz_fname2 + '.npz'))
            np.savez(npz_fname2, **iebdm_pc1)
