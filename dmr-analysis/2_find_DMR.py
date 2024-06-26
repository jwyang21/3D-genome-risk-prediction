#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import argparse
import pickle

CHR_LIST = [f'chr{i}' for i in range(1, 23)]
SCORE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('-b', '--binsize', help = 'bnsize', default = int(1e6), type = int, required = True)
    parser.add_argument('--cpg_type', help = 'cpg type. opensea/island/shelf/shore/shelf_shore', type = str, required = True)
    parser.add_argument('-c', '--cohort', type = str, required = True) 
    parser.add_argument('--dmr_type', help = 'DMR type. TN (tumor-normal), HL (high score - low score) risk_HL', type = str, required = True)
    parser.add_argument('--event', help = 'survival event', type = str, required = True)
    parser.add_argument('--fold', help = 'fold number', type = str, required = True, default = 'fold1')
    parser.add_argument('--version', help = 'feature version', type = str, required = True, default = 'v7.1')
    parser.add_argument('--lr', help = 'learning rate. 0.001, 0.0005, 0.0001', type = float, required = True, default = 0.001)
    
    return parser.parse_args()

def get_risk_HL_samples(cohort, version, cpg_type, fold, event, lr, s2l): 
    risk_dir = os.path.join('/data/project/3dith/result', cpg_type, cohort, version)
    risk_group_fname = f'{fold}_{event}_lr_{lr}_test_df_with_risk.csv'
    risk_df = pd.read_csv(os.path.join(risk_dir, risk_group_fname), index_col = 0)
    H = risk_df[risk_df['risk_group']=='High'].copy().index.values
    L = risk_df[risk_df['risk_group']=='Low'].copy().index.values
    T = risk_df.index.values
    
    H_out = [s2l[x] for x in H]
    L_out = [s2l[x] for x in L]
    T_out = [s2l[x] for x in T]
    
    return H_out, L_out, T_out

if __name__ == '__main__':
    
    args = parse_arguments()
    os.chdir(args.working_dir)
    
    s2l_fname = f'/data/project/3dith/data/{args.cohort}_s2l.pickle'
    with open(s2l_fname, 'rb') as f:
        s2l = pickle.load(f)
    f.close()

    result_dir = os.path.join(os.getcwd(), 'result')
    print("result_dir: {}".format(result_dir))
    
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
    
    colname_ = f'{args.dmr_type}_diff'
    mean_std_df_data = [] 
    chrom_list = []
    mean_std_df_column = ['chrom', f'{args.cohort}_{args.event}_{args.fold}']
    print("===\ncohort: ", args.cohort)    
    fname = os.path.join(SAVEDIR, f'binned_avg_{args.cpg_type}_{args.dmr_type}_diff_binsize_{args.binsize}.npz')
    f = np.load(fname, allow_pickle = True)
    globals()[f'{args.cohort}_DMR'] = {}
    
    for chrom in CHR_LIST:
        print("---\n"+chrom)
        v = f[chrom].copy()
        b = f[chrom+'_bins'].copy()
        current_chrom_df = pd.DataFrame(v, index = b, columns = [colname_])
        current_chrom_df.dropna(inplace = True)
        current_mean = np.mean(current_chrom_df[colname_].values.flatten())
        current_std = np.std(current_chrom_df[colname_].values.flatten())
        chrom_list.append(chrom)
        mean_std_df_data.append((current_mean - current_std))
        mean_std_mask = current_chrom_df[colname_].values.flatten() < (current_mean - current_std)
        globals()[f'{args.cohort}_DMR'][chrom+'_mean_std_bins'] = []
        mean_std_bins = (current_chrom_df.index.values * mean_std_mask)
        
        mean_std_bins_cnt = 0
        for x in mean_std_bins:
            if x != '':
                globals()[f'{args.cohort}_DMR'][chrom+'_mean_std_bins'].append(x)
                mean_std_bins_cnt += 1
                
        assert mean_std_bins_cnt == len(globals()[f'{args.cohort}_DMR'][chrom+'_mean_std_bins'])
        print("proportion_mean_std_bins: {}".format(mean_std_bins_cnt / current_chrom_df.shape[0]))

    print("---\nSave DMR bins of current cohort (per each chrom).")
    result_fname = os.path.join(SAVEDIR, 'DMR_binsize_'+str(args.binsize))
    print("DMR result file: {}".format(result_fname+'.npz'))
    np.savez(result_fname, **globals()[f'{args.cohort}_DMR'])
   
    mean_std_df = pd.DataFrame(zip(chrom_list, mean_std_df_data))
    mean_std_df.columns = mean_std_df_column
    print("===\nSave result of all cohorts (shape: num_chrom, num_cohorts).")
    print(os.path.join(SAVEDIR, f'chrom_cohort_mean_std_{args.dmr_type}_diff_binsize_{args.binsize}.csv'))
    mean_std_df.to_csv(os.path.join(SAVEDIR, f'chrom_cohort_mean_std_{args.dmr_type}_diff_binsize_{args.binsize}.csv'), index = True)
