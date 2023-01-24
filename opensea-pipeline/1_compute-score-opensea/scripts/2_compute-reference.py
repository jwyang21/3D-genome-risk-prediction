import argparse
import pandas as pd
import numpy as np
import os
import random
import sys
from scipy.integrate import simps
from numpy import trapz
random.seed(2022)
np.random.seed(2022)

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('--pc1_upper_dir', help = 'directory where pc1 directories are in', type = str, required = True)
    parser.add_argument('-r_dir', '--result_dir', help = 'result directory', type = str, required = True)
    parser.add_argument('--score_type', help = 'pc1-avg or pc1-fluctuation.', type = str, default = 'pc1-avg', required = True)
    parser.add_argument('--usage_option', help = 'whether you will use all or part of the available samples to compute reference pc1.', default = 'half', type = str, required = True)
    parser.add_argument('-m_type', '--matrix_type', help = 'from which matrix you will compute PC1. bdm, iebdm, or all.', type = str, required = True)
    parser.add_argument('--cohort', help = 'cohort. TCGA-{} or PCBC', type = str, required = True)
    parser.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
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

def standardize(v):
    return (v - v.mean()) / v.std()

def pick_random_samples(cohort, T, N, S, cohort_dir, fname):
    if 'TCGA' in cohort:
        sampled_N = random.sample(N, len(N)//2)
        excluded_N =list(set(N)-set(sampled_N))
    elif cohort=='PCBC':
        sampled_N = random.sample(S, len(S)//2) 
        excluded_N =list(set(S)-set(sampled_N))
    print("number of randomly picked normal samples: {}".format(len(sampled_N)))
    print("Number of excluded normal samples: {}".format(len(excluded_N)))

    np.savez(fname, picked = np.array(sampled_N), excluded = np.array(excluded_N))
    print("names of randomly-picked or excluded samples: {}".format(fname+'.npz'))
    return sampled_N, excluded_N

def import_pc1(pc1_dir, sample, chrom, standardize_flag, flag='inv'):

    if flag=='raw':
        fname = os.path.join(pc1_dir, sample+'.npz')
    elif flag=='inv':
        fname = os.path.join(pc1_dir, sample+'_inv_exp.npz')
    else:
        pass
    key_ = chrom+'_pc1'
    pc1 = np.load(fname)[key_]
    if standardize_flag=='Y':
        pc1 = standardize(pc1)
    return pc1

def compute_reference_avg_pc1(cohort, sampled_N, cohort_dir, pc1_dir, pc1_flag, m_type, score_type, usage_option, standardize_flag):
    reference_pc1 = {}
    standardize_flag_yn = '_standardized' if standardize_flag=='Y' else ''
    for chrom in CHR_LIST:
        for s in sampled_N:           
            current_pc1 = import_pc1(pc1_dir, s, chrom, standardize_flag, pc1_flag)
            if sampled_N.index(s)==0:
                globals()['N_pc1_'+chrom] = current_pc1.copy()
            else:
                globals()['N_pc1_'+chrom] = np.vstack((globals()['N_pc1_'+chrom], current_pc1)) 
        reference_pc1[chrom] = globals()['N_pc1_'+chrom].mean(axis=0) 
        del(globals()['N_pc1_'+chrom])
    if 'TCGA' in cohort:
        fname = 'normal-reference_'+m_type+'_'+score_type+'_'+usage_option+standardize_flag_yn
        save_fname = os.path.join(cohort_dir, fname)
    else: 
        fname = 'stem-reference_'+m_type+'_'+score_type+'_'+usage_option+standardize_flag_yn
        save_fname = os.path.join(cohort_dir, fname)
    return reference_pc1, save_fname

def compute_pc1_fluctuation(pc1): 
    pc1_abs = np.array([abs(x) for x in pc1])
    abs_area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return abs_area 

if __name__=='__main__':
    
    args = parse_arguments()

    T, N, S = get_sample_list(args.cohort)

    cohort_dir = os.path.join(args.result_dir, args.cohort)
    pc1_dir = os.path.join(args.pc1_upper_dir, args.cohort, 'pc1')
    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir) 

    pc1_flag = 'inv' if args.matrix_type == 'iebdm' else 'raw'

    if args.usage_option == 'half': 
        samples_fname = os.path.join(cohort_dir, 'picked-not-picked-samples')
        if os.path.exists(samples_fname):
            f = np.load(samples_fname)
            sampled_N, excluded_N = f['picked'], f['excluded']
        else:
            fname = os.path.join(cohort_dir, 'picked-not-picked-samples')
            sampled_N, excluded_N = pick_random_samples(args.cohort, T, N, S, cohort_dir, fname)

    if args.score_type == 'pc1-avg':
        reference_pc1, save_fname = compute_reference_avg_pc1(args.cohort, sampled_N, cohort_dir, pc1_dir, pc1_flag, 
            args.matrix_type, args.score_type, args.usage_option, args.standardize)


    elif args.score_type == 'pc1-fluctuation':
        pass

    np.savez(save_fname, **reference_pc1)
    print("result file: {}".format(save_fname+'.npz'))
    print("----") 


