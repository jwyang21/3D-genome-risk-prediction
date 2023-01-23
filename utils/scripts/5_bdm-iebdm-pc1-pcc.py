#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
import os
import glob
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

mpl.rcParams['figure.dpi'] = 300 #150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)


# In[3]:

# command example: python3 230118_bdm-iebdm-pc1.py --cpg_type opensea > ../log/230118_bdm-iebdm-pc1_opensea.log
cohorts = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[4]:


#cpg_type = 'opensea'


# In[8]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]


# In[37]:


'''
result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/1_compute-score-{cpg_type}/result/bdm-iebdm-pc1'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
'''


# In[32]:


SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples


# In[ ]:


def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', type = str, required = True, help = 'opensea, island, shelf_shore')
    return args.parse_args()


# In[ ]:


if __name__=='__main__':
    args = parse_arguments()
    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/1_compute-score-{args.cpg_type}/result/bdm-iebdm-pc1'
    print(f'result_dir: {result_dir}')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
        
    for c in cohorts:
        cohort_dir = os.path.join(result_dir, c)
        if not os.path.exists(cohort_dir):
            os.makedirs(cohort_dir)
        print(f'===\nCohort:{c}')
        cohort_sample_names = np.load(SAMPLE_NAME_FILE)[f'{c}_samples']
        pc1_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/1_compute-score-{args.cpg_type}/result/{c}/pc1'
        all_bdm_samples_flist = []
        all_iebdm_samples_flist = []

        for s in cohort_sample_names:
            s_list = glob.glob(os.path.join(pc1_dir, f'*{s}*'))
            for f in s_list:
                if f.endswith('.npz'):
                    if f.endswith('_inv_exp.npz'):
                        all_iebdm_samples_flist.append(f)
                    else:
                        all_bdm_samples_flist.append(f)


        assert len(all_bdm_samples_flist) == len(all_iebdm_samples_flist)
        for i in range(len(all_bdm_samples_flist)):
            assert os.path.basename(all_bdm_samples_flist[i]).split('.npz')[0] in all_iebdm_samples_flist[i]
        print(f"num_samples: {len(all_bdm_samples_flist)}")

        for chrom in CHR_LIST:
            print(f"---\n{chrom}")
            pcc_df = pd.DataFrame(index=cohort_sample_names, columns = ['pcc'])
            pcc_df_fname = f'{chrom}_bdm-iebdm-pcc.csv'
            chrom_key = f'{chrom}_pc1'
            for i in range(len(all_bdm_samples_flist)):
                bdm = all_bdm_samples_flist[i]
                iebdm = all_iebdm_samples_flist[i]
                sample_name = os.path.basename(bdm).split('.npz')[0]
                assert sample_name in iebdm

                bdm_pc1 = np.load(bdm)[chrom_key]
                iebdm_pc1 = np.load(iebdm)[chrom_key]
                
                pcc = pearsonr(bdm_pc1, iebdm_pc1)[0]
                pval = pearsonr(bdm_pc1, iebdm_pc1)[1]

                if pval <= 5e-2:
                    pcc_df.loc[f'{sample_name}']['pcc'] = pcc

            pcc_df.to_csv(os.path.join(cohort_dir, pcc_df_fname))
            print(os.path.join(cohort_dir, pcc_df_fname))
            if pcc_df.shape[0] == len(cohort_sample_names):
                print(f"ALL samples' {chrom} bdm-iebdm-pcc values are saved.")
                
    avg_pcc_df = pd.DataFrame(index = cohorts, columns = CHR_LIST)
    fig = plt.figure(figsize = (2*5, 2*3))
    #for i in range(1):
    for i in range(len(cohorts)):
        c = cohorts[i]
        ax = fig.add_subplot(3, 5, i+1)
        #fig = plt.figure(figsize = (3*6, 3*4))
        #globals()[f'{cohort}'] = []
        for chrom in CHR_LIST:
            #ax = fig.add_subplot(4, 6, i+1)
            df_name = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/1_compute-score-{args.cpg_type}/result/bdm-iebdm-pc1/{c}/{chrom}_bdm-iebdm-pcc.csv'
            df = pd.read_csv(df_name, index_col = 0)
            if df.isna().sum().sum() > 0:
                print(f'{df.isna().sum().sum()} samples are dropped in ({c}, {chrom}, {args.cpg_type}) because of insignificnat pcc with p value > 0.05')
                print('---\ndropped sample:')
                print(df.iloc[df.isna().values.flatten(),:].index.values)
                print("---")
            pcc = df.pcc.dropna().values.flatten()#drop insignificant pcc, whose pval > 5e-2
            #if df.isna().sum().sum():
            #    print("\n\n\n\n\n"+c+' '+chrom)
            #pcc = [abs(x) for x in pcc]
            #lobals()[f'{cohort}'][chrom] = pcc
            avg_pcc_df.loc[c][chrom] = np.mean(pcc)
        ax.hist(avg_pcc_df.loc[c].values.flatten())
        ax.set_title(c)
    fig.tight_layout()
    plt.savefig(os.path.join(result_dir, 'pcc-bdm-iebdm-ALL-cohorts.png'))
    print(os.path.join(result_dir, 'pcc-bdm-iebdm-ALL-cohorts.png'))
    avg_pcc_df.to_csv(os.path.join(result_dir, 'pcc-bdm-iebdm-ALL-cohorts.csv'))
    print(os.path.join(result_dir, 'pcc-bdm-iebdm-ALL-cohorts.csv'))


