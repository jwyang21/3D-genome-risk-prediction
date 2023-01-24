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
normal7_cohort = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'cpg type. opensea, island, shelf_shore', type = str, required = True)
    args.add_argument('--hg19_len_fname', help = 'chromosome sizes in hg19', type = str, required = False, default = '/data/project/jeewon/research/reference/hg19.fa.sizes')
    args.add_argument('--binsize', help = 'BDM binsize', type = int, required = True, default = int(1e6))
    args.add_argument('--tcga_fire_cohort_fname', help = 'description of TCGA cohorts intersecting with FIRE', type = str, required = False, default = '/data/project/3dith/data/etc/tcga-fire-cohorts.csv')
    args.add_argument('--hic_pc1_fname', help = 'Hi-C PC1 fname', type = str, required = False, default = '/data/project/3dith/data/fire_pc1.csv')
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

if __name__ == '__main__':
    args = parse_arguments()
    tcga_fire_cohort = pd.read_csv(args.tcga_fire_cohort_fname, index_col = 0)
    fire_normal7_cohort = np.intersect1d(tcga_fire_cohort.index.values, normal7_cohort)
    bdm_bins_dir = f'/data/project/3dith/data/bdm_bins/{args.cpg_type}'
    hg19_len = pd.read_csv(args.hg19_len_fname, sep = '\t', index_col = 0, header = None)
    hg19_len.columns = ['len']
    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/compare-hic-450k-pc1-{args.matrix_type}/normal-individual/fire'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    fire_pc1 = pd.read_csv(args.hic_pc1_fname)
    
    chroms = np.array([int(x) for x in fire_pc1.chr.values])
    fire_pc1['chr'] = chroms
    fire_indices = np.array([f'chr{fire_pc1.chr.values[i]}:{int(fire_pc1.start.values[i])-1}-{int(fire_pc1.start.values[i])+args.binsize-1}' for i in range(fire_pc1.shape[0])])
    fire_pc1.index = fire_indices


    for chrom in CHR_LIST:
        print(f"===\n{chrom}")

        for i in range(len(fire_normal7_cohort)):
            cohort = fire_normal7_cohort[i]

            cohort_FIRE_key = tcga_fire_cohort.loc[cohort]['FIRE']

            cohort_result_dir = os.path.join(result_dir, cohort)
            if not os.path.exists(cohort_result_dir):
                os.makedirs(cohort_result_dir)

            current_key = cohort.split('-')[1]+'_'+chrom

            globals()[f'{current_key}_all'] = {}
            globals()[f'{current_key}_all']['pcc'] = []
            globals()[f'{current_key}_all']['pval'] = []

            tissue = tcga_fire_cohort.loc[cohort]['FIRE_tissue'].lower()
            print(f"---\n{cohort}, {tissue}")

            cohort_bdm_bin_fname = os.path.join(bdm_bins_dir, f'{cohort}_diffmat_bins.npz')
            cohort_bdm_bin = np.load(cohort_bdm_bin_fname)
            cohort_bdm_chrom_bin = cohort_bdm_bin[f'{chrom}_bins']

            tmp = fire_pc1[['chr',cohort_FIRE_key]]
            tmp = tmp[tmp['chr'] == int(chrom.split('chr')[-1])].copy()
            tmp = tmp.dropna()
            tmp.drop(['chr'], axis = 1, inplace = True)

            fire_tcga_intersecting_bins = np.intersect1d(tmp.index.values, cohort_bdm_chrom_bin)
 
            hic_corrmat_pc1 = tmp.loc[fire_tcga_intersecting_bins][cohort_FIRE_key].values.flatten()
            hic_corrmat_pc1 = hic_corrmat_pc1 * (-1)
            hic_corrmat_pc1 = standardize(hic_corrmat_pc1)

            individual_450k_pc1_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/1_compute-score-{args.cpg_type}/result/{cohort}/pc1'
            
            if args.matrix_type == 'bdm':               
                all_npzs = glob.glob(os.path.join(individual_450k_pc1_dir, f'*.npz')) 
                individual_pc1_fnames = []
                for i in range(len(all_npzs)):
                    if 'inv_exp' not in all_npzs[i]:
                        individual_pc1_fnames.append(all_npzs[i])
            elif args.matrix_type == 'iebdm':
                individual_pc1_fnames = glob.glob(os.path.join(individual_450k_pc1_dir, f'*_inv_exp.npz'))
            else:
                pass
            
            individual_pc1_fname_basename = np.array([os.path.basename(x) for x in individual_pc1_fnames])
            individual_N_pc1_fnames = []

            for i in range(len(individual_pc1_fnames)):
                f = individual_pc1_fnames[i]
                b = individual_pc1_fname_basename[i]
                if int(b[13:15])>9:
                    individual_N_pc1_fnames.append(f)

            individual_N_pc1_fnames = np.array(individual_N_pc1_fnames)

            individual_samplenames = []
            for j in range(len(individual_N_pc1_fnames)):
                current_fname = individual_N_pc1_fnames[j]
                current_basename = os.path.basename(current_fname)
                if args.matrix_type == 'bdm':
                    current_sample = current_basename.split('.np')[0]
                elif args.matrix_type == 'iebdm':
                    current_sample = current_basename.split('_inv_exp')[0]
                else:
                    pass
                assert len(current_sample) == 15
                
                individual_samplenames.append(current_sample)
                print(current_fname)
                
                
                individual_N_pc1_all = np.load(current_fname)
                chrom_key = f'{chrom}_pc1'

                tmp = individual_N_pc1_all[chrom_key]
                tmp_df = pd.DataFrame(tmp)

                tmp_df.index = cohort_bdm_chrom_bin
                tmp = tmp_df.loc[fire_tcga_intersecting_bins].values.flatten()
                individual_N_pc1 = standardize(tmp)

                pcc = pearsonr(hic_corrmat_pc1, individual_N_pc1)[0]
                pval = pearsonr(hic_corrmat_pc1, individual_N_pc1)[1]
                globals()[f'{current_key}_all']['pcc'].append(pcc)
                globals()[f'{current_key}_all']['pval'].append(pval)

                fig = plt.figure(figsize = (8, 1.75))
                ax = fig.add_subplot(111)
                ax.axhline(0, lw=0.75, ls='--', c='0.8')
                y1 = hic_corrmat_pc1.copy()
                ax.plot(smoothen(y1, window=3), lw=1, c='0.5', label='Hi-C ({})'.format(tissue))
                
                y2 = individual_N_pc1.copy()
                ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimation ({}, 450K)'.format(cohort))

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
                fig_basename = f'individual_N_vector_{chrom}_{current_sample}-w-hic.png'
                fig.tight_layout()
                plt.savefig(os.path.join(cohort_result_dir, fig_basename), dpi = 300)
                print(os.path.join(cohort_result_dir, fig_basename))
                plt.clf()

            current_result_all = pd.DataFrame.from_dict(globals()[f'{current_key}_all'])
            current_result_all.index = individual_samplenames
            current_result_all.to_csv(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'), index = True)
            print(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'))
