import pandas as pd
import numpy as np
import os
import argparse

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='shelf_shore'

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working_dir', type = str, required = True)
    #default: /data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}

    args.add_argument('--score_dir', help = 'directory where scores are saved in', type = str, required = True)
    #default: /data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result#/{cohort}

    args.add_argument('--result_dir', help = 'directory to save resulting files', type = str, required = True)
    #default: /data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result#/{cohort}

    args.add_argument('--usage_option', help = 'use all samples or randomly picked samples when computing reference.', type = str, default = 'part', required = True) # all, part
    #default: part

    args.add_argument('--manifest_fname', help = 'manifest file name', type = str, required = True)
    # default: /data/project/3dith/data/manifest_normal7.csv

    
    return args.parse_args()

def compute_min_max_score(result_dir, score_fname, save_full_fname, all_cohorts):
    score_df = pd.DataFrame(np.zeros((len(all_cohorts), 2), dtype = float), index = all_cohorts, columns = ['min', 'max'])
    for cohort in all_cohorts:
        cohort_dir = os.path.join(result_dir, cohort)
        score = pd.read_csv(os.path.join(cohort_dir, score_fname), index_col = 0)
        score_df.loc[cohort]['min'] = np.min(score.values)
        score_df.loc[cohort]['max'] = np.max(score.values)
    score_df.to_csv(save_full_fname)
    print("result: {}".format(save_full_fname))

if __name__=='__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    NORMAL7_COHORT = pd.read_csv(args.manifest_fname).cohort.values.tolist()
    NORMAL7_COHORT_PCBC = NORMAL7_COHORT.copy()
    NORMAL7_COHORT_PCBC.append('PCBC')
    
    distances = ['normal-distance', 'stem-distance']
    metrics = ['euclidean', 'cosine-sim']
    #score_types = ['pc1-avg', 'pc1-fluctuation']
    distance_avgs = ['simple-avg', 'weighted-avg']
    matrix_types = ['iebdm', 'bdm']
    standardize_options = ['_standardized', '']
    num_chroms = ['_to-chrom-'+str(i) for i in np.arange(1, 22)]
    num_chroms.append('')
    assert len(num_chroms)==22

    score_type = 'pc1-avg'

    for distance in distances:
        for metric in metrics:
            for distance_avg in distance_avgs:
                for matrix_type in matrix_types:
                    for standardize_option in standardize_options:
                        for num_chrom in num_chroms:
                            score_fname = distance+'_'+metric+'_'+matrix_type+'_'+score_type+'_'+distance_avg+'_'+args.usage_option+standardize_option+num_chrom+'.csv'
                            save_fname = 'minmax_'+score_fname
                            save_full_fname = os.path.join(args.result_dir, save_fname)
                            if distance=='normal-distance':
                                compute_min_max_score(args.result_dir, score_fname, save_full_fname, NORMAL7_COHORT)    
                            else: #stem-distance
                                compute_min_max_score(args.result_dir, score_fname, save_full_fname, NORMAL7_COHORT_PCBC)
