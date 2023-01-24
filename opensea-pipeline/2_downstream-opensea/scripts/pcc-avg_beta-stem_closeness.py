import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import argparse

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    args.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    args.add_argument('--score_fname', help = 'stem-closeness score file name, based on opensea K-M result', type = str, default = '/data/project/3dith/data/cohort-1-best-score-km.csv', required = True)
    args.add_argument('--cpg_type', help = 'open/island/shelf_shore', type = str, required = True)
    return args.parse_args()

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

def get_cpg_list(chr_list, cpg_type):
    if cpg_type == 'opensea':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', sep = '\t', header = None) 

    elif cpg_type == 'island':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.island.sorted.bed', sep = '\t', header = None) 

    elif cpg_type == 'shelf_shore':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.shelf_shore.sorted.bed', sep = '\t', header = None) 
    else: 
        raise Exception("Wrong cpg_type!")

    chr_list_mask = np.array([x in chr_list for x in df.iloc[:,0].values.flatten()])
    df_chr_list = df.iloc[chr_list_mask,:]

    for x in df_chr_list.iloc[:,0].values.flatten():
        assert x in chr_list

    total_list = df_chr_list.iloc[:,-1].values.flatten()
    return total_list

def get_avg_beta(cohort, cpg_list, S): 
    beta_fname = '/data/project/3dith/data/450k_xena/'+cohort+'.HumanMethylation450.tsv' 
    beta = pd.read_csv(beta_fname, sep = '\t', index_col = 0)
    beta_target_cpg_df = beta.loc[cpg_list].dropna() 
    print("num_target_cpg_after_dropna: {}".format(beta_target_cpg_df.shape[0]))
    samples_beta = beta_target_cpg_df.columns.values
    avg_target_cpg_beta = beta_target_cpg_df.mean().values 
    avg_beta_target_cpg_df = pd.DataFrame(avg_target_cpg_beta, index = samples_beta, columns = ['avg_beta'])

    return beta_target_cpg_df[S], avg_beta_target_cpg_df.loc[S] 


if __name__=='__main__':
    args = parse_arguments()
    
    os.chdir(args.working_dir)
    result_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'result')
    cohort_dir = os.path.join(result_dir, args.cohort)  
    if not os.path.exists(result_dir):
        os.makedirs(result_dir) 
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    print("cohort: {}".format(args.cohort))
    T, N, S = get_sample_list(args.cohort)

    cpg_list = get_cpg_list(CHR_LIST, args.cpg_type) 
    print("num_{}_cpg (for 22 autosomes): {}".format(args.cpg_type, len(cpg_list)))

    beta_target_cpg_df, avg_beta_target_cpg_df = get_avg_beta(args.cohort, cpg_list, S) 
    beta_target_cpg_df_fname = args.cpg_type+'_all_samples_all_beta.csv'
    avg_beta_target_cpg_df_fname = args.cpg_type+'_all_samples_avg_beta.csv'

    beta_target_cpg_df_Tumor_fname = args.cpg_type+'_tumors_all_beta.csv'
    avg_beta_target_cpg_df_Tumor_fname = args.cpg_type+'_tumors_avg_beta.csv'

    beta_target_cpg_df.to_csv(os.path.join(cohort_dir, beta_target_cpg_df_fname))
    avg_beta_target_cpg_df.to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname))

    if 'TCGA' in args.cohort:
        beta_target_cpg_df[T].copy().to_csv(os.path.join(cohort_dir, beta_target_cpg_df_Tumor_fname))
        avg_beta_target_cpg_df.loc[T].copy().to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_Tumor_fname))
        print(os.path.join(cohort_dir, beta_target_cpg_df_Tumor_fname))
        print(os.path.join(cohort_dir, avg_beta_target_cpg_df_Tumor_fname))
    else:
        pass

    if 'ESCA' in args.cohort:
        sys.exit(0)
    elif 'HNSC' in args.cohort:
        sys.exit(0)
    else:
        all_score_df = pd.read_csv(args.score_fname, index_col = 0)
        current_score_fname = str(all_score_df.loc[args.cohort][f'filename_{args.cpg_type}'])
        score_df = pd.read_csv(current_score_fname, index_col = 0)
        assert args.cohort in current_score_fname

        print("---\nPCC(avg_beta({}), stem-closeness) in {}, using all samples (Tumor + Normal)".format(args.cpg_type, args.cohort))
        pcc, pvalue = pearsonr(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten())
        print("pcc: {}, p-value: {}".format(pcc, pvalue))

        print("---\nPCC(avg_beta({}), stem-closeness) in {}, using tumor samples".format(args.cpg_type, args.cohort))
        pcc, pvalue = pearsonr(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten())
        print("pcc: {}, p-value: {}".format(pcc, pvalue))

        fig = plt.figure(figsize = (3,3))
        ax = fig.add_subplot(111)
        ax.scatter(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten())
        ax.set_xlabel('average {} beta value'.format(args.cpg_type))
        ax.set_ylabel('stem-closeness')

        ax.set_title(f'{args.cohort} (all samples)')
        fig.tight_layout()
        plt.savefig(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.png'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.png'))
        figure_data = pd.DataFrame(zip(avg_beta_target_cpg_df.loc[S].values.flatten(), score_df.loc[S].cos_radian.values.flatten()), index = S, columns = ['avg_beta', 'stem_closeness'])
        figure_data.to_csv(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.csv'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-ALL_samples.csv'))

        fig = plt.figure(figsize = (3,3))
        ax = fig.add_subplot(111)
        ax.scatter(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten())
        ax.set_xlabel('average {} beta value'.format(args.cpg_type))
        ax.set_ylabel('stem-closeness')
        ax.set_title(f'{args.cohort} (tumor samples)')
        plt.savefig(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.png'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.png'))
        figure_data = pd.DataFrame(zip(avg_beta_target_cpg_df.loc[T].values.flatten(), score_df.loc[T].cos_radian.values.flatten()), index = T, columns = ['avg_beta', 'stem_closeness'])
        figure_data.to_csv(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.csv'))
        print(os.path.join(cohort_dir, 'scatter-avg_beta-stem_closeness-TUMOR_samples.csv'))

        plt.clf
