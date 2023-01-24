import numpy as np
import pandas as pd
import os
import sys
from sklearn.decomposition import PCA
import argparse

CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='opensea'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('--tcga_bdm_dir', help = 'TCGA binned difference matrix directory', type = str, required = True)
    parser.add_argument('--pcbc_bdm_dir', help = 'PCBC binned difference matrix directory', type = str, required = True)
    parser.add_argument('-r_dir', '--result_dir', help = 'result directory', type = str, required = True)
    parser.add_argument('-m_type', '--matrix_type', help = 'from which matrix you will compute PC1. bdm, iebdm, or all.', type = str, required = True)
    parser.add_argument('--cohort', help = 'cohort. TCGA-{} or PCBC', type = str, required = True)
    return parser.parse_args()

def get_sample_list(cohort):
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = [] 
    else: 
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: 
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S

def import_bdm(directory, chrom, s): 
    f = os.path.join(directory, s+'.npz')
    diffmat = np.load(f)[chrom]
    diffmat_mask = np.load(f)[chrom+'_mask']
    diffmat_masked = diffmat[~diffmat_mask].T[~diffmat_mask].T

    return diffmat_masked

def pc1(m):
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    pc1 = pc[:,0]
    return pc1

if __name__ == '__main__':
    args = parse_arguments()
    
    os.chdir(args.working_dir)

    if 'TCGA' in args.cohort:
        bdm_dir = args.tcga_bdm_dir
        data_dir = os.path.join(bdm_dir, args.cohort)
    else:
        bdm_dir = args.pcbc_bdm_dir
        data_dir = bdm_dir 

    print("===\ncohort: {}".format(args.cohort))
    
    cohort_dir = os.path.join(args.result_dir, args.cohort) 

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
            
            npz_fname1 = os.path.join(cohort_dir,'pc1', s) 
            print("npz_fname: {}".format(npz_fname1+'.npz'))
            np.savez(npz_fname1, **bdm_pc1)
        
        elif args.matrix_type == 'iebdm':
            iebdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = 1/np.exp(import_bdm(data_dir, chrom, s))
                iebdm_pc1[current_key] = pc1(m)
            
            npz_fname2 = os.path.join(cohort_dir, 'pc1', s + '_inv_exp') 
            print("npz_fname: {}".format(npz_fname2 + '.npz'))
            np.savez(npz_fname2, **iebdm_pc1)
        
        elif args.matrix_type == 'all':
            bdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = import_bdm(data_dir, chrom, s)
                bdm_pc1[current_key] = pc1(m)
           
            npz_fname1 = os.path.join(cohort_dir, 'pc1', s) 
            print("npz_fname: {}".format(npz_fname1+'.npz'))
            np.savez(npz_fname1, **bdm_pc1)

            iebdm_pc1 = {}
            for chrom in CHR_LIST:
                current_key = chrom+'_pc1'
                m = 1/np.exp(import_bdm(data_dir, chrom, s))
                iebdm_pc1[current_key] = pc1(m)
           
            npz_fname2 = os.path.join(cohort_dir, 'pc1', s + '_inv_exp') 
            print("npz_fname: {}".format(npz_fname2 + '.npz'))
            np.savez(npz_fname2, **iebdm_pc1)
