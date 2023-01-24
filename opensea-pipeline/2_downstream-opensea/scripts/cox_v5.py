import argparse
import os
import pandas as pd
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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
    parser.add_argument('--score_file', type = str, default = '', required=False)
    parser.add_argument('--score_group_file', type = str, default = '', required = False)
    parser.add_argument('--avg_beta_file', type = str, default = '', required=False)
    parser.add_argument('--beta_group_file', type = str, default = '', required = False)
    parser.add_argument('--clinical_file', default='/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv')
    
    return parser.parse_args()

def print_args(args):
    for arg in vars(args):
        print(f'{arg}: {getattr(args, arg)}')

def mkdir(dir_name):
    if not os.path.exists(dir_name):
        os.system(f'mkdir -p {dir_name}')
        print(f'\nDirectory {dir_name} is created.')

def prep_clinical(args):
    clinical = pd.read_csv(args.clinical_file)
    usecols = ['bcr_patient_barcode', 'age_at_initial_pathologic_diagnosis', 'gender', f'{args.survival_type}', f'{args.survival_type}.time']
        
    clinical = clinical[clinical['type'] == args.cohort.split('-')[1]][usecols]
    clinical.dropna(inplace=True)
    
    if args.use_score=='Y':
        assert args.score_file != ''
        stem_closeness = pd.read_csv(args.score_file, index_col=0)[['cos_radian']]
        stem_closeness.columns = ['stem closeness']
    avg_beta = pd.read_csv(args.avg_beta_file, index_col=0)
    avg_beta.columns = ['avg beta']
    if args.use_beta_group == 'Y':
        assert args.beta_group_file != ''
        beta_group = pd.read_csv(args.beta_group_file, index_col = 0)[['group_name']]
        beta_group.columns = ['beta group']
    if args.use_score_group == 'Y':
        assert args.score_group_file != ''
        score_group = pd.read_csv(args.score_group_file, index_col = 0)[['group_name']]
        score_group.columns = ['score group']
        
    sample_one_idx = []
    sample_one_barcode = []
    sample_more_idx = []
    sample_more_barcodes = []
    
    for idx in clinical['bcr_patient_barcode'].index:
        row = clinical.loc[idx]
        sample_short = row['bcr_patient_barcode']
        sample_long_list = avg_beta.index[avg_beta.index.str.startswith(sample_short)].to_list()
        if len(sample_long_list) == 1:
            sample_one_idx.append(idx)
            sample_one_barcode.append(sample_long_list[0])
        else:
            sample_more_idx.extend([idx] * len(sample_long_list))
            sample_more_barcodes.extend(sample_long_list)
            
    clinical = pd.concat([clinical.loc[sample_one_idx], clinical.loc[sample_more_idx]], axis=0)
    clinical.index = sample_one_barcode + sample_more_barcodes
    if args.use_score=='Y':
        clinical = clinical.merge(stem_closeness, left_index=True, right_index=True)
    if args.use_avg_beta=='Y':
        clinical = clinical.merge(avg_beta, left_index=True, right_index=True)
    if args.use_beta_group == 'Y':
        clinical = clinical.merge(beta_group, left_index=True, right_index=True)
    if args.use_score_group == 'Y':
        clinical = clinical.merge(score_group, left_index=True, right_index=True)
    
    return clinical

def perform_CoxPH(args, clinical): 
    gender_dict = {'MALE': 0, 'FEMALE': 1}

    clinical['gender'] = clinical['gender'].map(gender_dict) 
    clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 
    group_dict = {'Low': 0, 'High': 1}

    if args.use_beta_group == 'Y':
        clinical['beta group'] = clinical['beta group'].map(group_dict) 
        clinical.rename(columns = {'beta group':'beta group=High'},inplace=True)
    if args.use_score_group == 'Y':
        clinical['score group'] = clinical['score group'].map(group_dict) 
        clinical.rename(columns = {'score group':'score group=High'},inplace=True)
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
