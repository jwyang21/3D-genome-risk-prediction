#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[ ]:


NORMAL7_COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
# # scatter (average opensea beta, stem_closeness)
# remove ESCA and HNSC since score is not finalized in these cohorts. 
NORMAL7_COHORTS.remove('TCGA-ESCA')
NORMAL7_COHORTS.remove('TCGA-HNSC')

si_rna_types = ['mRNAsi', 'EREG-mRNAsi']
si_dna_types = ['mDNAsi', 'EREG-mDNAsi', 'DMPsi', 'ENHsi']
all_si_types = si_rna_types + si_dna_types

all_df_columns = ['pcc_'+x for x in all_si_types]
all_df =  pd.DataFrame(np.zeros((len(NORMAL7_COHORTS), len(all_df_columns))), index = NORMAL7_COHORTS, columns = all_df_columns)


# In[72]:

#command: python3 4_pcc-sc-si.py --cpg_type opensea --si_fname /data/project/3dith/data/stemness-index/1-s2.0-S0092867418303581-mmc1.xlsx
def parse_arguments():
    args = argparse.ArgumentParser()
    
    args.add_argument('--cpg_type', help = 'island, opensea, shelf_shore', type = str, required = True)
    
    args.add_argument('--si_fname', help = 'stemness index filename, written in absolute path', type = str, required = True)
    #default: /data/project/3dith/data/stemness-index/1-s2.0-S0092867418303581-mmc1.xlsx
    
    args.add_argument('--score_fname', help = 'stem-closeness score file name, based on opensea K-M result', type = str, default = '/data/project/3dith/data/cohort-1-best-score-km.csv', required = True)
    
    return args.parse_args()  


# In[ ]:


if __name__=='__main__':
    args = parse_arguments()
    
    working_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type)
    os.chdir(working_dir)
    
    print("0. assign save directory")
    save_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'result', 'compare-si-sc')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    print("save_dir: {}".format(save_dir))
    
    print("1. import stemness index")
    si_rna = pd.read_excel(args.si_fname, index_col = 0, sheet_name = 'StemnessScores_RNAexp')
    si_dna = pd.read_excel(args.si_fname, index_col = 0, sheet_name = 'StemnessScores_DNAmeth')
    
    print("2. preprocess stemness index")
    si_rna_index = [x[:15] for x in si_rna.index.values]
    si_dna_index = [x[:15] for x in si_dna.index.values]
    si_rna['index'] = si_rna_index
    si_dna['index'] = si_dna_index
    si_rna.drop_duplicates(subset = 'index', inplace = True)
    si_dna.drop_duplicates(subset = 'index', inplace = True)
    si_rna.index = si_rna['index'].values
    si_dna.index = si_dna['index'].values
    si_rna.drop(['index'], axis = 1, inplace = True)
    si_dna.drop(['index'], axis = 1, inplace = True)
    
    print(f'2. Import best {args.cpg_type} score of each cohort')
    #scores_fname = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'result', 'cohort-1-best-score-km.csv')
    scores_cohort = pd.read_csv(args.score_fname, index_col = 0)
    
    print("3. compare stem closeness (sc) and stemness index (si)")
    
    for si_type in all_si_types:
        print("===\nsi_type: {}".format(si_type))
        fig = plt.figure(figsize = (2*7, 2*2)) #13 cohorts-of-interest

        for i in range(len(NORMAL7_COHORTS)):
            cohort = NORMAL7_COHORTS[i]
            print("---\ncohort: {}".format(cohort))

            sc_filename = scores_cohort.loc[cohort][f'filename_{args.cpg_type}']
            sc = pd.read_csv(sc_filename, index_col = 0)

            if si_type in si_rna_types:
                si_df = si_rna.copy()
            else:
                si_df = si_dna.copy()

            si_df = si_df[si_df['cancer.type']==cohort.split('TCGA-')[1]]

            intersecting_sample = np.intersect1d(si_df.index.values, sc.index.values)

            sc = sc.loc[intersecting_sample]
            si_df = si_df.loc[intersecting_sample]

            sc_si_merged = pd.merge(si_df, sc, left_index = True, right_index = True)

            sc_si_merged = sc_si_merged[[sc.columns[-1], si_type]]

            pcc = pearsonr(sc_si_merged.cos_radian.values, sc_si_merged[si_type].values.flatten())[0]
            all_df.loc[cohort][f'pcc_{si_type}'] = pcc
            #display(sc_si_merged)   
            ax = fig.add_subplot(2, 7, i+1)
            ax.scatter(sc_si_merged.cos_radian.values, sc_si_merged[si_type].values.flatten())
            ax.set_xlabel('stem closeness')
            ax.set_ylabel(si_type)
            num_samples = sc_si_merged.shape[0]
            ax.set_title(f'{cohort} (n = {num_samples})',fontsize = 8)
            #ax.text(0, 0, str(pcc))
            #plt.scatter(sc_si_merged.cos_radian.values, sc_si_merged.mRNAsi.values)
            sc_si_merged_fname = f'scatter-{si_type}-stem_closeness-{cohort}-{num_samples}.csv'
            sc_si_merged.to_csv(os.path.join(save_dir, sc_si_merged_fname))
            print(os.path.join(save_dir, sc_si_merged_fname))

        plt.suptitle(f'scatter(stem closeness, {si_type})', fontsize = 11)
        fig.tight_layout()

        fig_name = f'scatter-{si_type}-stem_closeness.png'
        plt.savefig(os.path.join(save_dir, fig_name))
        print("---")
        print(os.path.join(save_dir, fig_name))
    all_df.to_csv(os.path.join(save_dir, 'pcc-sc-si.csv'))
    print(os.path.join(save_dir, 'pcc-sc-si.csv'))
    
    print("4. histogram of pcc(stem closeness, stemness index) for each stemness index")
    fig = plt.figure(figsize = (1.5*all_df.shape[1], 1.5))
    for i in range(all_df.shape[1]):
        ax = fig.add_subplot(1, all_df.shape[1], i+1)
        ax.hist(all_df.iloc[:,i].values.flatten())
        si_type = all_df.columns[i].split('pcc_')[1]

        ax.set_title(si_type)
    suptitle_ = 'histogram of pcc(stem closeness, stemness index)'
    plt.suptitle(suptitle_)
    fig.tight_layout()
    fig_title = 'pcc-sc-si.png'
    plt.savefig(os.path.join(save_dir, fig_title))
    print(os.path.join(save_dir, fig_title))
    
    plt.clf()
