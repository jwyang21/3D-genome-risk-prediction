import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr
import glob
from sklearn.decomposition import PCA
import argparse

CHR_LIST = [f'chr{i}' for i in range(1, 23)]

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'cpg type. opensea, island, shelf_shore', type = str, required = True)
    args.add_argument('--matching_df_fname', help = '3div and TCGA matching result', type = str, required = False, default = '/data/project/3dith/data/etc/3div-cancer_hic-tcga-matching.csv')
    args.add_argument('--hic_corrmat_dir', help = 'Hi-C correlation matrix directory', type = str, required = True, default = '/data/project/3dith/data/hic_corrmat')
    args.add_argument('--hg19_len_fname', help = 'chromosome sizes in hg19', type = str, required = False, default = '/data/project/jeewon/research/reference/hg19.fa.sizes')
    args.add_argument('--binsize', help = 'BDM binsize', type = int, required = True, default = int(1e6))
    args.add_argument('--cohort', help = 'TCGA cohort', type = str, required = True)
    args.add_argument('--matrix_type', help = 'bdm or iebdm', type = str, required = True)
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
    matching_df = pd.read_csv(args.matching_df_fname, index_col = 0).loc[[args.cohort]].copy()
    bdm_bins_dir = f'/data/project/3dith/data/bdm_bins/{args.cpg_type}'
    hg19_len = pd.read_csv(args.hg19_len_fname, sep = '\t', index_col = 0, header = None)
    hg19_len.columns = ['len']
    hic_corrmat_pc1_dir = os.path.join(args.hic_corrmat_dir, 'pc1')

    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/compare-hic-450k-pc1-{args.matrix_type}/tumor-individual'

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    for chrom in CHR_LIST:
        print(f"===\n{chrom}")
        for i in range(matching_df.shape[0]):
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
            individual_T_pc1_fnames = []
            
            for i in range(len(individual_pc1_fnames)):
                f = individual_pc1_fnames[i]
                b = individual_pc1_fname_basename[i]
                if int(b[13:15])<=9:
                    individual_T_pc1_fnames.append(f)
            
            individual_T_pc1_fnames = np.array(individual_T_pc1_fnames)
            individual_samplename = []
            for j in range(len(individual_T_pc1_fnames)):
                
                current_fname = individual_T_pc1_fnames[j]
                current_basename = os.path.basename(current_fname)
                if args.matrix_type == 'bdm':
                    current_sample = current_basename.split('.np')[0]
                elif args.matrix_type == 'iebdm':
                    current_sample = current_basename.split('_inv_exp')[0]
                
                assert len(current_sample) == 15
                
                individual_samplename.append(current_sample)
                print(current_fname)
                
                individual_T_pc1_all = np.load(current_fname)
                chrom_key = f'{chrom}_pc1'
                individual_T_pc1 = individual_T_pc1_all[chrom_key]
                individual_T_pc1 = standardize(individual_T_pc1)
                
                pcc = pearsonr(hic_corrmat_pc1, individual_T_pc1)[0]
                pval = pearsonr(hic_corrmat_pc1, individual_T_pc1)[1]
                globals()[f'{current_key}_all']['pcc'].append(pcc)
                globals()[f'{current_key}_all']['pval'].append(pval)

            current_result_all = pd.DataFrame.from_dict(globals()[f'{current_key}_all'])
            current_result_all.index = individual_samplename
            current_result_all.to_csv(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'), index = True)
            print(os.path.join(cohort_result_dir, f'{chrom}_all_standardized.csv'))
