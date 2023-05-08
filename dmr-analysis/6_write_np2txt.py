#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import sys
import pickle
import argparse

SCORE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
THRESHOLD = ['mean', 'mean_std', 'mean_2std']
SMALL_CATEGORY_FNAME = '/data/project/3dith/data/etc/dmr-feature-small-category.npz'
if os.path.exists(SMALL_CATEGORY_FNAME):
    SMALL_CATEGORY = np.load(SMALL_CATEGORY_FNAME, allow_pickle = True)
else:
    SMALL_CATEGORY = {}
    SMALL_CATEGORY['GENE'] = ['gene','transcript']
    SMALL_CATEGORY['REG'] = ['open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', 'enhancer', 'promoter',                             'promoter_flanking_region']
    SMALL_CATEGORY['EPI'] = ['18_Quies', '13_Het', '17_ReprPCWk', '16_ReprPC', '14_TssBiv', '2_TssFlnk',     '12_ZNF/Rpts', '11_EnhWk', '1_TssA', '6_TxWk', '5_Tx', '9_EnhA1', '7_EnhG1',     '4_TssFlnkD' ,'15_EnhBiv', '10_EnhA2', '3_TssFlnkU', '8_EnhG2']
    np.savez(SMALL_CATEGORY_FNAME,**SMALL_CATEGORY)
cohort2eid = pd.read_csv('/data/project/3dith/data/etc/cohort2eid.txt', sep = '\t', header = None)
cohort2eid.columns = ['cohort', 'eid']
eid_cohorts = cohort2eid.cohort.values
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', help = 'DMR threshold. mean, mean_std, mean_2std', type = str, required = True, default = 'mean_std')
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('--dmr_type', type = str, help = 'TN or HL', required = True)
    parser.add_argument('--event', help = 'survival event', type = str, required = True)
    parser.add_argument('--fold', help = 'fold number', type = str, required = True, default = 'fold1')
    parser.add_argument('--version', help = 'feature version', type = str, required = True, default = 'v7.1')
    parser.add_argument('--lr', help = 'learning rate. 0.001, 0.0005, 0.0001', type = float, required = True, default = 0.001)    
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    return parser.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')

    npz_fname = 'DMR_GENE_features_threshold_'+args.threshold+'_id_len.npz'
    full_npz_fname = os.path.join(SAVEDIR, npz_fname)
    current_ensg = np.load(full_npz_fname)['ENSG_all']
    current_enst = np.load(full_npz_fname)['ENST_all']
    result_ensg_fname = os.path.join(SAVEDIR, 'ENSG_'+args.threshold+'.txt')
    result_enst_fname = os.path.join(SAVEDIR, 'ENST_'+args.threshold+'.txt')

    f = open(result_ensg_fname, 'w')
    for i in range(len(current_ensg)):
        f.write(current_ensg[i])
        if i != len(current_ensg)-1:
            f.write("\n")
    f.close()
    print(result_ensg_fname)

    f = open(result_enst_fname, 'w')
    for i in range(len(current_enst)):
        f.write(current_enst[i])
        if i != len(current_enst)-1:
            f.write("\n")
    f.close()
    print(result_enst_fname)
