#!/usr/bin/env python
# coding: utf-8

# In[93]:


import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr
import glob
from sklearn.decomposition import PCA
import argparse

'''
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
'''

# tumor sample 수가 너무 많아서 그림 저장 시 에러가 발생 -> 맨 처음 돌릴 때는 pcc 값들만 저장 후, 다 돌리고 나서 pcc 값 확인해서 가장 성능 좋은 idnividual을 backtrack해서 그 individual만 따로 plot하자.

# ## cancer hi-c (3div) 데이터로 만든 corr.mat에서 뽑은 PC1과 450K IEBDM에서 뽑은 PC1을 비교.
# - 450K PC1으로 cancer hi-c PC1이 재현이 잘 돼야 좋음.
CHR_LIST = [f'chr{i}' for i in range(1, 23)]
#matching_df = pd.read_csv('/data/project/3dith/data/etc/3div-cancer_hic-tcga-matching.csv', index_col = 0)
#hic_corrmat_dir = '/data/project/3dith/data/hic_corrmat'
#hg19_len_fname = '/data/project/jeewon/research/reference/hg19.fa.sizes'
#binsize = int(1e6)


def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'cpg type. opensea, island, shelf_shore', type = str, required = True)
    args.add_argument('--matching_df_fname', help = '3div and TCGA matching result', type = str, required = False, default = '/data/project/3dith/data/etc/3div-cancer_hic-tcga-matching.csv')
    args.add_argument('--hic_corrmat_dir', help = 'Hi-C correlation matrix directory', type = str, required = True, default = '/data/project/3dith/data/hic_corrmat')
    args.add_argument('--hg19_len_fname', help = 'chromosome sizes in hg19', type = str, required = False, default = '/data/project/jeewon/research/reference/hg19.fa.sizes')
    args.add_argument('--binsize', help = 'BDM binsize', type = int, required = True, default = int(1e6))
    args.add_argument('--cohort', help = 'TCGA cohort', type = str, required = True)
    return args.parse_args()

def pc1(m):
    #print("pc1(m)")
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    #print(pc)
    pc1 = pc[:,0]
    #print(pc1)
    #print('-----')
    return pc1

def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()

if __name__=='__main__':
    args = parse_arguments()
    #matching_df = pd.read_csv(args.matching_df_fname, index_col = 0)#여러 cohorts 동시에
    matching_df = pd.read_csv(args.matching_df_fname, index_col = 0).loc[[args.cohort]].copy()#한개 cohort만.
    bdm_bins_dir = f'/data/project/3dith/data/bdm_bins/{args.cpg_type}'#/{cohort}_diffmat_bins.npz. #item: chr{i}_bins
    hg19_len = pd.read_csv(args.hg19_len_fname, sep = '\t', index_col = 0, header = None)
    hg19_len.columns = ['len']
    hic_corrmat_pc1_dir = os.path.join(args.hic_corrmat_dir, 'pc1')

    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/compare-hic-450k-pc1/tumor-individual'

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    for chrom in CHR_LIST:#real
    #for chrom in CHR_LIST[:1]:#debug
        print(f"===\n{chrom}")
        for i in range(matching_df.shape[0]):#real
        #for i in range(1):#debug
            cohort = matching_df.index.values[i]
            cohort_result_dir = os.path.join(result_dir, cohort)
            if not os.path.exists(cohort_result_dir):
                os.makedirs(cohort_result_dir)

            current_key = cohort.split('-')[1]+'_'+chrom

            globals()[f'{current_key}_all'] = {}
            globals()[f'{current_key}_all']['pcc'] = []
            globals()[f'{current_key}_all']['pval'] = []
            tissue = matching_df.tissue.values[i]
            print(f"---\n{cohort}, {tissue}")
            cohort_bdm_bin_fname = os.path.join(bdm_bins_dir, f'{cohort}_diffmat_bins.npz')
            cohort_bdm_bin = np.load(cohort_bdm_bin_fname)
            cohort_bdm_chrom_bin = cohort_bdm_bin[f'{chrom}_bins']

            n_total_chrom_bins = ( int(hg19_len.loc[f'{chrom}']['len']) // args.binsize) + 1
            n_total_chrom_bin_names = np.array([f'{chrom}:{str(i * args.binsize)}-{str((i+1) * args.binsize)}' for i in range(n_total_chrom_bins)])
            bdm_bin_mask = np.array([x in cohort_bdm_chrom_bin for x in n_total_chrom_bin_names])

            hic_corrmat_fname = os.path.join(args.hic_corrmat_dir, f'{tissue}.hg19_1mb.{chrom}.corrmat.npy')
            hic_corrmat = np.load(hic_corrmat_fname)
            hic_corrmat_masked = hic_corrmat[bdm_bin_mask].T[bdm_bin_mask].T.copy()
            hic_corrmat_pc1 = pc1(hic_corrmat_masked)

            hic_corrmat_pc1_basename = os.path.basename(hic_corrmat_fname).split('.npy')[0]+'.pc1.npy'
            hic_corrmat_pc1_fname = os.path.join(hic_corrmat_pc1_dir, hic_corrmat_pc1_basename)

            #np.save(hic_corrmat_pc1_fname, hic_corrmat_pc1) #save raw PC1
            hic_corrmat_pc1 = standardize(hic_corrmat_pc1)
            
            individual_450k_pc1_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/1_compute-score-{args.cpg_type}/result/{cohort}/pc1'
            individual_pc1_fnames = glob.glob(os.path.join(individual_450k_pc1_dir, f'*_inv_exp.npz'))
            individual_pc1_fname_basename = np.array([os.path.basename(x) for x in individual_pc1_fnames])
            individual_T_pc1_fnames = []
            
            for i in range(len(individual_pc1_fnames)):
                f = individual_pc1_fnames[i]
                b = individual_pc1_fname_basename[i]
                if int(b[13:15])<=9:
                    individual_T_pc1_fnames.append(f)
            
            individual_T_pc1_fnames = np.array(individual_T_pc1_fnames)
            
            #repr_450k_pc1_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/repr_vectors/{cohort}'
            #repr_T_vector_pc1_fnames = glob.glob(os.path.join(repr_450k_pc1_dir, f'repr_T_vector_{chrom}_*.npy')) 
                #이 안에는 (각 repr vector를 정의하는 데 필요한 샘플 수, number of bdm bins in this chromosome of this cohort) shape의 행렬이 있음.
                #column mean해서 쓰기.
            for j in range(len(individual_T_pc1_fnames)):
                current_fname = individual_T_pc1_fnames[j]
                current_basename = os.path.basename(current_fname)
                current_sample = current_basename.split('_inv_exp')[0]
                
                assert len(current_sample) == 15
                
                print(current_fname)
                
                individual_T_pc1_all = np.load(current_fname)
                chrom_key = f'{chrom}_pc1'
                individual_T_pc1 = individual_T_pc1_all[chrom_key]
                individual_T_pc1 = standardize(individual_T_pc1)
                
                pcc = pearsonr(hic_corrmat_pc1, individual_T_pc1)[0]
                pval = pearsonr(hic_corrmat_pc1, individual_T_pc1)[1]
                globals()[f'{current_key}_all']['pcc'].append(pcc)
                globals()[f'{current_key}_all']['pval'].append(pval)

                
                '''
                # plot inidvidual PC1 and hi-c PC1.
                #fig = plt.figure(figsize = (6, 0.75))               
                fig = plt.figure(figsize = (8, 1.75))
                ax = fig.add_subplot(111)
                ax.axhline(0, lw=0.75, ls='--', c='0.8')
                
                y1 = hic_corrmat_pc1.copy()
                ax.plot(smoothen(y1, window=3), lw=1, c='0.5', label='Hi-C ({})'.format(tissue))
                
                y2 = individual_T_pc1.copy()
                #ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimations ({}, averaged)'.format(cohort))#repr vector
                ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimation ({}, 450K)'.format(cohort))#individual vector

                for direction in ['top', 'right', 'bottom']:
                    ax.spines[direction].set_visible(False)
                ax.spines['left'].set_position(('outward', 3))
                ax.spines['left'].set_position(('outward', 3))

                assert len(y1)==len(y2)
                ax.set_xlim([0, len(y1)])
                ax.set_ylim([-2.5, 2.5])

                ax.set_xlabel('Position ({})'.format(chrom))
                ax.set_ylabel('PC1\n(Standardized)')

                ax.set_xticks([])
                mb_length = 1 / len(y1)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.1, -0.1], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.1, '10Mb', transform=ax.transAxes, ha='left', va='center')
                ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.4))
                
                fig_basename = f'individual_T_vector_{chrom}_{current_sample}-w-hic.png'
                fig.tight_layout()
                plt.savefig(os.path.join(cohort_result_dir, fig_basename), dpi = 300)
                plt.clf()
                '''
            current_result_all = pd.DataFrame.from_dict(globals()[f'{current_key}_all'])
            current_result_all.to_csv(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'), index = False)
            print(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'))
