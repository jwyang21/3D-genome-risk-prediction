import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from sklearn.decomposition import PCA
import argparse

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

CHR_LIST = [f'chr{i}' for i in range(1, 23)]

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'cpg type. opensea, island, shelf_shore', type = str, required = True)
    args.add_argument('--matching_df_fname', help = '3div and TCGA matching result', type = str, required = False, default = '/data/project/3dith/data/etc/3div-normal_hic-tcga-matching.csv')
    args.add_argument('--hic_corrmat_dir', help = 'Hi-C correlation matrix directory', type = str, required = True, default = '/data/project/3dith/data/hic_corrmat')
    args.add_argument('--hg19_len_fname', help = 'chromosome sizes in hg19', type = str, required = False, default = '/data/project/jeewon/research/reference/hg19.fa.sizes')
    args.add_argument('--binsize', help = 'BDM binsize', type = int, required = True, default = int(1e6))
    args.add_argument('--cohort', help = 'TCGA cohort', type = str, required = True)
    args.add_argument('--matrix_type', help = 'bdm of iebdm', type = str, required = True)
    return args.parse_args()

def pc1(m):
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    pc1 = pc[:,0]
    return pc1

def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()

if __name__=='__main__':
    args = parse_arguments()
    matching_df = pd.read_csv(args.matching_df_fname, index_col = 0).loc[args.cohort].copy()
    bdm_bins_dir = f'/data/project/3dith/data/bdm_bins/{args.cpg_type}'#/{cohort}_diffmat_bins.npz. #item: chr{i}_bins
    hg19_len = pd.read_csv(args.hg19_len_fname, sep = '\t', index_col = 0, header = None)
    hg19_len.columns = ['len']
    hic_corrmat_pc1_dir = os.path.join(args.hic_corrmat_dir, 'pc1')

    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/compare-hic-450k-pc1-{args.matrix_type}/normal-repr/3div'

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    for chrom in CHR_LIST:
        print(f"===\n{chrom}")
        for i in range(matching_df.shape[0]):
            cohort = matching_df.index.values[i]
            cohort_result_dir = os.path.join(result_dir, cohort+'_'+str(i+1))
            if not os.path.exists(cohort_result_dir):
                os.makedirs(cohort_result_dir)

            current_key = cohort.split('-')[1]+'_'+chrom
            globals()[f'{current_key}_repr'] = {}
            globals()[f'{current_key}_repr']['pcc'] = []
            globals()[f'{current_key}_repr']['pval'] = []
            
            tissue = matching_df.tissue.values[i]
            print(f"---\n{cohort}, {tissue}")
            cohort_bdm_bin_fname = os.path.join(bdm_bins_dir, f'{cohort}_diffmat_bins.npz')
            cohort_bdm_bin = np.load(cohort_bdm_bin_fname)
            cohort_bdm_chrom_bin = cohort_bdm_bin[f'{chrom}_bins']

            n_total_chrom_bins = ( int(hg19_len.loc[f'{chrom}']['len']) // args.binsize) + 1
            n_total_chrom_bin_names = np.array([f'{chrom}:{str(i * args.binsize)}-{str((i+1) * args.binsize)}' for i in range(n_total_chrom_bins)])
            bdm_bin_mask = np.array([x in cohort_bdm_chrom_bin for x in n_total_chrom_bin_names])

            hic_corrmat_fname = os.path.join(args.hic_corrmat_dir, f'{tissue}{str(i+1)}.hg19_1mb.{chrom}.corrmat.npy')
            hic_corrmat = np.load(hic_corrmat_fname)
            hic_corrmat_masked = hic_corrmat[bdm_bin_mask].T[bdm_bin_mask].T.copy()
            hic_corrmat_pc1 = pc1(hic_corrmat_masked)

            hic_corrmat_pc1_basename = os.path.basename(hic_corrmat_fname).split('.npy')[0]+'.pc1.npy'
            hic_corrmat_pc1_fname = os.path.join(hic_corrmat_pc1_dir, hic_corrmat_pc1_basename)

            np.save(hic_corrmat_pc1_fname, hic_corrmat_pc1) 
            hic_corrmat_pc1 = standardize(hic_corrmat_pc1)

            repr_450k_pc1_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/repr_vectors-{args.matrix_type}/{cohort}'
            repr_N_vector_pc1_fnames = glob.glob(os.path.join(repr_450k_pc1_dir, f'repr_N_vector_{chrom}_*.npy')) 

            repr_samples = []
            for j in range(len(repr_N_vector_pc1_fnames)):
                current_fname = repr_N_vector_pc1_fnames[j]
                print(current_fname)
                repr_N_pc1_10_vectors = np.load(current_fname)
                repr_N_pc1 = np.load(current_fname).mean(axis = 0) 
                assert len(hic_corrmat_pc1) == len(repr_N_pc1)
                repr_samples.append(os.path.basename(current_fname).split('.np')[0])

                repr_N_pc1 = standardize(repr_N_pc1)
                pcc = pearsonr(hic_corrmat_pc1, repr_N_pc1)[0]
                pval = pearsonr(hic_corrmat_pc1, repr_N_pc1)[1]
                globals()[f'{current_key}_repr']['pcc'].append(pcc)
                globals()[f'{current_key}_repr']['pval'].append(pval)

                fig = plt.figure(figsize = (8, 1.75))
                ax = fig.add_subplot(111)
                ax.axhline(0, lw=0.75, ls='--', c='0.8')
                y1 = hic_corrmat_pc1.copy()
                ax.plot(smoothen(y1, window=3), lw=1, c='0.5', label='Hi-C ({})'.format(tissue))
                y2 = repr_N_pc1.copy()
                ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimations ({}, 450K, averaged)'.format(cohort))
                for k in range(len(repr_N_pc1_10_vectors)):
                    y3 = repr_N_pc1_10_vectors[k].copy()
                    y3 = standardize(y3)
                    ax.plot(smoothen(y3, window=3), c='C3', lw=0.75, alpha=0.2)

                for direction in ['top', 'right', 'bottom']:
                    ax.spines[direction].set_visible(False)
                ax.spines['left'].set_position(('outward', 3))
                ax.spines['left'].set_position(('outward', 3))

                assert len(y1)==len(y2) and len(y2)==len(y3)
                ax.set_xlim([0, len(y1)])
                ax.set_ylim([-2.5, 2.5])

                ax.set_xlabel('Position ({})'.format(chrom))
                ax.set_ylabel('PC1\n(Standardized)')

                ax.set_xticks([])
                mb_length = 1 / len(y1)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.1, -0.1], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.1, '10Mb', transform=ax.transAxes, ha='left', va='center')
                ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.4))
                fig_basename = f'repr_N_vector_{chrom}_{str(j)}-w-hic.png'
                fig.tight_layout()
                plt.savefig(os.path.join(cohort_result_dir, fig_basename), dpi = 300)
                plt.clf()

            current_result_repr = pd.DataFrame.from_dict(globals()[f'{current_key}_repr'])
            
            current_result_repr.index = repr_samples

            current_result_repr.to_csv(os.path.join(cohort_result_dir, f'{chrom}_repr_standardized.csv'), index = True)

            print(os.path.join(cohort_result_dir, f'{chrom}_repr_standardized.csv'))

