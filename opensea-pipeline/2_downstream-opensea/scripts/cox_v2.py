#!/usr/bin/env python
# coding: utf-8

# In[1]:

#v3: age, gender, score_group
import argparse
import os
import pandas as pd
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


# In[2]:


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)


# In[3]:


unavailable = ['[Unknown]', '[Not Available]', '[Not Evaluated]', '[Not Applicable]', '[Discrepancy]'] 
# items that should be treated equal to NaN 
# not needed when using only age and gender from TCGA clinical data.


# In[ ]:


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--survival_type', required=True)
    parser.add_argument('--result_dir', required=True)
    parser.add_argument('--result_file', required=True)

    parser.add_argument('--use_score', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_score_group', type = str, help = 'Y or N', required = True)
    parser.add_argument('--use_avg_beta', type = str, help = 'Y or N', required = True)
    
    parser.add_argument('--score_file', type = str, default = '', required=False)
    parser.add_argument('--score_group_file', type = str, default = '', required = False)
    parser.add_argument('--avg_beta_file', type = str, default = '', required=True)
    
    parser.add_argument('--clinical_file', default='/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv')
    
    return parser.parse_args()


# In[ ]:


def print_args(args):
    for arg in vars(args):
        print(f'{arg}: {getattr(args, arg)}')


# In[ ]:


def mkdir(dir_name):
    if not os.path.exists(dir_name):
        os.system(f'mkdir -p {dir_name}')
        print(f'\nDirectory {dir_name} is created.')


# In[ ]:


def prep_clinical(args):
    clinical = pd.read_csv(args.clinical_file)
    usecols = ['bcr_patient_barcode', 'age_at_initial_pathologic_diagnosis', 'gender', f'{args.survival_type}', f'{args.survival_type}.time']
        
    clinical = clinical[clinical['type'] == args.cohort.split('-')[1]][usecols]
    clinical.dropna(inplace=True)
    
    if args.use_score=='Y':
        assert args.score_file != ''
        stem_closeness = pd.read_csv(args.score_file, index_col=0)[['cos_radian']]
    if args.use_score_group == 'Y':
        assert args.score_group_file != ''
        score_group = pd.read_csv(args.score_group_file, index_col = 0)[['group_name']]
        score_group.columns = ['score group']
    if args.use_avg_beta == 'Y':
        assert args.avg_beta_file != ''
        avg_beta = pd.read_csv(args.avg_beta_file, index_col=0)
        avg_beta.columns = ['avg beta']
        
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
    if args.use_score_group == 'Y':
        clinical = clinical.merge(score_group, left_index=True, right_index=True)
    #print('clinical:')
    #print(clinical.head(3))
    
    return clinical


# In[ ]:


def perform_CoxPH(args, clinical):
    # gender는 categorical 변수이므로 one-hot vector로 변경. 
    
    gender_dict = {'MALE': 0, 'FEMALE': 1}

    clinical['gender'] = clinical['gender'].map(gender_dict) 
    clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) #'gender' column 이름을 'gender=FEMALE'로 바꾸기.

    # group은 categorical 변수이므로 one-hot vector로 변경. 
    group_dict = {'Low': 0, 'High': 1}
    if args.use_score_group == 'Y':
        clinical['score group'] = clinical['score group'].map(group_dict) 
    clinical.rename(columns = {'score group':'score group=High'},inplace=True)#'group' column 이름을 'group=High'로 바꾸기.
    
    clinical.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'},inplace=True)#'age_at_initial_pathologic_diagnosis' column 이름을 'age'로 바꾸기.

    clinical.drop(columns='bcr_patient_barcode', inplace=True)

        
    events = clinical[f'{args.survival_type}'].astype(bool)
    
    # drop columns with low variance conditioned on the death event. 
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
    '''
    # related warning log.
    # the warning below occurred in TCGA-BRCA, DFI. 
    /data/project/jeewon/miniconda3/envs/3dith/lib/python3.10/site-packages/lifelines/utils/__init__.py:1122: ConvergenceWarning: Column gender=FEMALE have very low variance when conditioned on death event present or not. This may harm convergence. This could be a form of 'complete separation'. For example, try the following code:

    >>> events = df['DFI'].astype(bool)
    >>> print(df.loc[events, 'gender=FEMALE'].var())
    >>> print(df.loc[~events, 'gender=FEMALE'].var())

    A very low variance means that the column gender=FEMALE completely determines whether a subject dies or not. See https://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression.
    '''
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


# ---

# # test in jupyter

# In[4]:


'''
cohort = 'TCGA-BLCA'
survival_type = 'OS'
version = '2'
result_dir = f'/data/project/jeewon/cox_v{version}/{cohort}'
result_file = survival_type

use_score = 'N'
use_score_group = 'Y'
use_avg_beta= 'Y'

score_file = '/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BLCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-8.csv'
score_group_file = '/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/sc-group.csv'
avg_beta_file = '/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/opensea_tumors_avg_beta.csv'

clinical_file = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
'''


# # prep_clinical

# In[5]:


'''
#def prep_clinical(args):
clinical = pd.read_csv(clinical_file)

usecols = ['bcr_patient_barcode', 'age_at_initial_pathologic_diagnosis', 'gender', f'{survival_type}', f'{survival_type}.time']
clinical = clinical[clinical['type'] == cohort.split('-')[1]][usecols]
clinical.dropna(inplace=True)

if use_score=='Y':
    assert score_file != ''
    stem_closeness = pd.read_csv(score_file, index_col=0)[['cos_radian']]
if use_score_group == 'Y':
    assert score_group_file != ''
    score_group = pd.read_csv(score_group_file, index_col = 0)[['group_name']]
    score_group.columns = ['group']
if use_avg_beta == 'Y':
    assert avg_beta_file != ''
    avg_beta = pd.read_csv(avg_beta_file, index_col=0)
    avg_beta.columns = ['avg beta']

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
if use_score=='Y':
    clinical = clinical.merge(stem_closeness, left_index=True, right_index=True)
if use_avg_beta=='Y':
    clinical = clinical.merge(avg_beta, left_index=True, right_index=True)
if use_score_group == 'Y':
    clinical = clinical.merge(score_group, left_index=True, right_index=True)
'''


# In[6]:


#display(clinical)


# ## perform_coxPH

# In[7]:


'''
#def perform_CoxPH(args, clinical):
# gender는 categorical 변수이므로 one-hot vector로 변경. 
gender_dict = {'MALE': 0, 'FEMALE': 1}

clinical['gender'] = clinical['gender'].map(gender_dict) 
clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) #'gender' column 이름을 'gender=FEMALE'로 바꾸기.

# group은 categorical 변수이므로 one-hot vector로 변경. 
group_dict = {'Low': 0, 'High': 1}

clinical['group'] = clinical['group'].map(group_dict) 
clinical.rename(columns = {'group':'group=High'},inplace=True)#'group' column 이름을 'group=High'로 바꾸기.

clinical.drop(columns='bcr_patient_barcode', inplace=True)
'''


# In[8]:


#display(clinical)


# In[7]:


'''
cph = CoxPHFitter()
cph.fit(df=clinical, duration_col=f'{survival_type}.time', event_col=f'{survival_type}')

result_file = f'{result_dir}/{result_file}'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
cph.summary.to_csv(f'{result_file}.csv')
print(f'\nResult is saved to {result_file}.csv.')

plt.figure(figsize=(1.5,1.5))
cph.plot()
plt.savefig(f'{result_file}.png', dpi=350, bbox_inches='tight')
print(f'\nFigure is saved to {result_file}.png.')
'''


# In[ ]:


def main():
    args = parse_args()
    print_args(args)
    
    clinical = prep_clinical(args)
    perform_CoxPH(args, clinical)


# In[ ]:


if __name__ == '__main__':
    main()

