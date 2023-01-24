import pandas as pd
import numpy as np
import os
import argparse

CHR_LIST = [f'chr{i}' for i in range(1, 23)]
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('-b', '--binsize', help = 'bnsize', default = int(1e6), type = int, required = True)
    parser.add_argument('--cpg_type', help = 'cpg type. opensea/island/shelf/shore/shelf_shore', type = str, required = True)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    result_dir = os.path.join(os.getcwd(), 'result')
    print("result_dir: {}".format(result_dir))

    mean_std_df = pd.DataFrame(np.zeros((len(CHR_LIST), len(NORMAL7_COHORT)), dtype = float), index = CHR_LIST, columns = NORMAL7_COHORT)

    for cohort in NORMAL7_COHORT:
        print("===\ncohort: ", cohort)
        cohort_dir = os.path.join(result_dir, cohort)
        print("cohort_dir: {}".format(cohort_dir))
        globals()[cohort+'_DMR'] = {} 
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        if not os.path.exists(cohort_dir):
            os.makedirs(cohort_dir)
        fname = os.path.join(cohort_dir, 'binned_avg_'+args.cpg_type+'_TN_diff_binsize_'+str(args.binsize)+'.npz')
        f = np.load(fname, allow_pickle = True)

        for chrom in CHR_LIST:
            print("---\n"+chrom)
            v = f[chrom].copy()
            b = f[chrom+'_bins'].copy()
            current_chrom_df = pd.DataFrame(v, index = b, columns = ['TN_diff'])
            current_chrom_df.dropna(inplace = True)
            current_mean = np.mean(current_chrom_df.TN_diff.values)
            current_std = np.std(current_chrom_df.TN_diff.values)
            mean_std_df.loc[chrom][cohort] = (current_mean - current_std).copy()
            mean_std_mask = current_chrom_df.TN_diff.values < (current_mean - current_std)
            globals()[cohort+'_DMR'][chrom+'_mean_std_bins'] = []
            mean_std_bins = (current_chrom_df.index.values * mean_std_mask)
            mean_std_bins_cnt = 0
            for x in mean_std_bins:
                if x != '':
                    globals()[cohort+'_DMR'][chrom+'_mean_std_bins'].append(x)
                    mean_std_bins_cnt += 1
            assert mean_std_bins_cnt == len(globals()[cohort+'_DMR'][chrom+'_mean_std_bins'])
            print("proportion_mean_std_bins: {}".format(mean_std_bins_cnt / current_chrom_df.shape[0]))
        print("---\nSave DMR bins of current cohort (per each chrom).")
        result_fname = os.path.join(cohort_dir, 'DMR_binsize_'+str(args.binsize))
        print("DMR result file: {}".format(result_fname+'.npz'))
        np.savez(result_fname, **globals()[cohort+'_DMR'])

    print("===\nSave result of all cohorts (shape: num_chrom, num_cohorts).")
    print(os.path.join(result_dir, 'chrom_cohort_mean_std_TN_diff_binsize_'+str(args.binsize)+'.csv'))
    mean_std_df.to_csv(os.path.join(result_dir, 'chrom_cohort_mean_std_TN_diff_binsize_'+str(args.binsize)+'.csv'), index = True)
