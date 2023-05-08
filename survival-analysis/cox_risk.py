#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import pandas as pd
from lifelines import CoxPHFitter
import matplotlib as mpl
mpl.use('AGG')
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

unavailable = ['[Unknown]', '[Not Available]', '[Not Evaluated]', '[Not Applicable]', '[Discrepancy]'] 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--survival_type', required=True)
    parser.add_argument('--result_dir', required=True)
    parser.add_argument('--result_file', required=True)
    parser.add_argument('--use_score', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_score_group', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_avg_beta', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_beta_group', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_risk', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_risk_group', type = str, help = 'Y or N', required = True)    
    parser.add_argument('--score_file', type = str, default = '', required=False) 
    parser.add_argument('--avg_beta_file', type = str, default = '', required=False)     
    parser.add_argument('--risk_file', type = str, default = '', required = False)
    parser.add_argument('--clinical_file', default='/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv') #short barcode
    
    return parser.parse_args()

def print_args(args):
    for arg in vars(args):
        print(f'{arg}: {getattr(args, arg)}')

def mkdir(dir_name):
    if not os.path.exists(dir_name):
        os.system(f'mkdir -p {dir_name}')
        print(f'\nDirectory {dir_name} is created.')

def drop_dupl(score_df):
    score_df['barcode'] = [x[:12] for x in score_df.index.values]
    score_df2 = score_df.drop_duplicates(subset = ['barcode']).copy() 
    score_df2.index = score_df2.barcode.values
    score_df2.drop(['barcode'], axis = 1, inplace=True)
    return score_df2

def add_column_to_clinical_df(clinical_df, score_df):
    intersection = np.intersect1d(clinical_df.index.values, score_df.index.values)
    clinical_df2 = clinical_df.loc[intersection].copy().sort_index()
    score_df2 = score_df.loc[intersection].copy().sort_index() # clinical과 겹치는 애들만 남김
    merged = pd.merge(clinical_df2, score_df2, left_index=True, right_index=True)
    return merged  

def prep_clinical(args):    
    clinical = pd.read_csv(args.clinical_file)
    usecols = ['bcr_patient_barcode', 'age_at_initial_pathologic_diagnosis', 'gender', f'{args.survival_type}', f'{args.survival_type}.time']        
    clinical_cohort = clinical[clinical['type'] == args.cohort.split('-')[1]][usecols]
    clinical_cohort.dropna(inplace=True)
    clinical_cohort.index = clinical_cohort['bcr_patient_barcode'].values.flatten()

    if args.use_score=='Y':
        assert args.score_file != ''
        stem_closeness = pd.read_csv(args.score_file, index_col=0)[['stem_closeness']]
        stem_closeness = drop_dupl(stem_closeness)
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, stem_closeness)
        
    if args.use_score_group == 'Y':
        assert args.score_file != ''
        score_group = pd.read_csv(args.score_file, index_col = 0)[['group']]
        score_group.columns = ['score group']
        score_group = drop_dupl(score_group)
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, score_group)
        
    if args.use_avg_beta == 'Y':
        assert args.avg_beta_file != ''
        avg_beta = pd.read_csv(args.avg_beta_file, index_col = 0)[['avg_beta']]
        avg_beta.columns = ['avg beta']
        avg_beta = drop_dupl(avg_beta)
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, avg_beta)
        
    if args.use_beta_group == 'Y':
        assert args.avg_beta_file != ''
        beta_group = pd.read_csv(args.avg_beta_file, index_col = 0)[['group']]
        beta_group.columns = ['beta group']
        beta_group = drop_dupl(beta_group)
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, beta_group)
    
    if args.use_risk == 'Y':
        assert args.risk_file != ''
        risk = pd.read_csv(args.risk_file, index_col = 0)[['risk']]
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, risk)
    
    if args.use_risk_group == 'Y':
        assert args.risk_file != ''
        risk_group_dict = {'High': 1, 'Low': 0}
        risk_group = pd.read_csv(args.risk_file, index_col = 0)[['risk_group']]
        risk_group['risk_group'] = risk_group['risk_group'].map(risk_group_dict) 
        risk_group.columns = ['risk group']
        clinical_cohort = add_column_to_clinical_df(clinical_cohort, risk_group)
    
    return clinical_cohort

def perform_CoxPH(args, clinical):    
    gender_dict = {'MALE': 0, 'FEMALE': 1}

    clinical['gender'] = clinical['gender'].map(gender_dict) 
    clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 
    
    if args.use_beta_group == 'Y':
        clinical.rename(columns = {'beta group':'beta group=High'},inplace=True)
    if args.use_score_group == 'Y':
        clinical.rename(columns = {'score group':'score group=High'},inplace=True)
    if args.use_risk_group == 'Y':
        clinical.rename(columns = {'risk group': 'risk group=High'}, inplace = True) 
        
    clinical.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'},inplace=True)
    clinical.drop(columns='bcr_patient_barcode', inplace=True)        
    events = clinical[f'{args.survival_type}'].astype(bool)
    
    low_var_col = []
    for col in clinical.columns:
        if col != f'{args.survival_type}':
            var1 = clinical.loc[events, f'{col}'].var()
            var2 = clinical.loc[events, f'{col}'].var()
            if var1 == 0 or var2 == 0:
                if col not in low_var_col:
                    low_var_col.append(col)
    if len(low_var_col) > 0:
        print("columns with low variance when conditioned on death event: {}".format(low_var_col))
        clinical.drop(low_var_col, axis = 1, inplace = True)
    print("final covariates: {}".format(clinical.columns))
    
    cph = CoxPHFitter()
    cph.fit(df=clinical, duration_col=f'{args.survival_type}.time', event_col=f'{args.survival_type}')

    result_file = f'{args.result_dir}/{args.result_file}'
    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
    cph.summary.to_csv(f'{args.result_dir}/{args.result_file}.csv')
    print(f'\nResult is saved to {args.result_dir}/{args.result_file}.csv.')

    plt.figure(figsize=(1.5,1.5))
    cph.plot()
    plt.savefig(f'{args.result_dir}/{args.result_file}.png', dpi=350, bbox_inches='tight')
    print(f'\nFigure is saved to {args.result_dir}/{args.result_file}.png.')
    
def main():
    args = parse_args()
    print_args(args)    
    clinical = prep_clinical(args)
    perform_CoxPH(args, clinical)

if __name__ == '__main__':
    main()
