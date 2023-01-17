import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

plt.rc('font', family='FreeSans', size=7)
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cohort', required=True)
    parser.add_argument('--survival_type', required=True)
    parser.add_argument('--result_dir', required=True)
    parser.add_argument('--result_file', required=True)

    parser.add_argument('--score_file', required=True)
    parser.add_argument('--avg_beta_file', required=True)

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

    stem_closeness = pd.read_csv(args.score_file, index_col=0)[['cos_radian']]
    avg_beta = pd.read_csv(args.avg_beta_file, index_col=0)

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

    clinical = clinical.merge(stem_closeness, left_index=True, right_index=True)
    clinical = clinical.merge(avg_beta, left_index=True, right_index=True)

    return clinical

def perform_CoxPH(args, clinical):
    gender_dict = {'MALE': 0, 'FEMALE': 1}

    clinical['gender'] = clinical['gender'].map(gender_dict)
    clinical.drop(columns='bcr_patient_barcode', inplace=True)
    clinical.columns = ['age', 'gender=FEMALE', f'{args.survival_type}', f'{args.survival_type}.time', 'avg beta', 'stem closeness']
    
    if args.cohort == 'TCGA-UCEC':
        clinical.drop(columns='gender=FEMALE', inplace=True)
        
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
    '''
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
    mkdir(args.result_dir)
    cph.summary.to_csv(f'{result_file}.csv')
    print(f'\nResult is saved to {result_file}.csv.')

    plt.figure(figsize=(1.5,1.5))
    cph.plot()
    plt.savefig(f'{result_file}.png', dpi=350, bbox_inches='tight')
    print(f'\nFigure is saved to {result_file}.png.')

def main():
    args = parse_args()
    print_args(args)

    clinical = prep_clinical(args)
    perform_CoxPH(args, clinical)

if __name__ == '__main__':
    main()
