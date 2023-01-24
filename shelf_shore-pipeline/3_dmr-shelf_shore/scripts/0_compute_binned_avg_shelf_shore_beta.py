import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
import argparse

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    args.add_argument('--cpg_type', help = 'CpG type. open/island/shelf/shore/shelf_shore', type = str, required = True)
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

if __name__ == '__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    SAVEDIR = os.path.join(os.getcwd(), 'result')
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    for cohort_ in NORMAL7_COHORT:
        cohort=cohort_ 
        print("====")

        print("1. get sample name list")
        cohort_dir = os.path.join(SAVEDIR, cohort)
        print(cohort_dir)
        print("cohort: {}".format(cohort))
        print("cohort_dir: ", cohort_dir)
        if not os.path.exists(cohort_dir):
            os.makedirs(cohort_dir)
        T, N, S = get_sample_list(cohort) 

        print("----\n2. find {} cpg probes of 22 autosomes".format(args.cpg_type))
        cpg_list = get_cpg_list(CHR_LIST, args.cpg_type) 
        print("num_{}_cpg: {}".format(args.cpg_type, len(cpg_list)))
        print("----\n3. import beta value corresponding to current cpg_type and average them")
        beta_target_cpg_df, avg_beta_target_cpg_df = get_avg_beta(cohort, cpg_list, S) 
        beta_target_cpg_df_fname = 'beta_'+args.cpg_type+'_df.csv'
        avg_beta_target_cpg_df_fname = 'avg_beta_'+args.cpg_type+'_df.csv'
        beta_target_cpg_df.to_csv(os.path.join(cohort_dir, beta_target_cpg_df_fname), index = True)
        avg_beta_target_cpg_df.to_csv(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname), index = True)
        print(os.path.join(cohort_dir, beta_target_cpg_df_fname))
        print(os.path.join(cohort_dir, avg_beta_target_cpg_df_fname))
