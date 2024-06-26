import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import random
import math
import os
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
random.seed(2022)
np.random.seed(2022)

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='opensea'
P_THRESHOLD = 5e-2

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working_dir', type = str, required = True)
    args.add_argument('--result_dir', help = 'directory to save resulting files', type = str, required = True)
    args.add_argument('--usage_option', help = 'use all samples or randomly picked samples when computing reference.', type = str, default = 'half', required = False) 
    args.add_argument('--score_type', help = 'pc1-avg or pc1-fluctuation', default = 'pc1-avg', type = str,  required = False)
    args.add_argument('--normalize', help = 'whether you wnat to normalize score2 and score4 or not', type = str, default = 'N', required = False)#Y or N
    args.add_argument('--minmax', help = 'use minmax scaling', default = 'N', type = str, required = False) #Y or N
    args.add_argument('--matrix_type', help = 'bdm or iebdm', type = str, default = 'iebdm', required = True)
    args.add_argument('--distance', help = 'distance metric. euclidean or cosine-sim', type = str, required = True)
    args.add_argument('--use_weighted_avg', help = 'whether to use weighted average in computing distance. Y or N (simple_avg)', type = str, required = True)
    args.add_argument('--minmax_file_dir', help = 'directory where file containing minmax scores of each cohort is in.', type = str, required = True)
    args.add_argument('--cohort', type = str, required = True)
    args.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
    args.add_argument('--num_chrom', help = 'chr1 to which chromosome you will use.', type = int, default = 22, required = False)
    
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

def compute_angle(x, y, cohort):
    if len(x) != len(y):
        raise ValueError
    
    for i in range(len(x)):
        if x[i] == 0:
            x[i] += 1e-15 #pseudocount to prevent division by zero. 
    
    slope = np.array([float(y[i]/x[i]) for i in range(len(x))])

    negative_slope_count = np.array([slope[i] < 0 for i in range(len(slope))]).sum()

    if negative_slope_count > 0:
        raise Exception("Negtive slope n {}".format(cohort))

    radian_theta = np.array([math.atan(s) for s in slope])
    degree_theta = np.array([math.degrees(r) for r in radian_theta])
    cosine_radian = np.array([math.cos(r) for r in radian_theta])

    return radian_theta, degree_theta, cosine_radian

if __name__ == '__main__':
    args = parse_arguments()

    T, N, S = get_sample_list(args.cohort) 

    print("cohort: {}".format(args.cohort))
    
    cohort_dir = os.path.join(args.result_dir, args.cohort)

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)

    weighted_avg_flag = 'weighted-avg' if args.use_weighted_avg=='Y' else 'simple-avg'
    standardize_flag = '_standardized' if args.standardize=='Y' else ''
    num_chrom_flag = '_to-chrom-'+str(args.num_chrom) if args.num_chrom < 22 else ''

    print("0. import min and max distances")
    print(args.minmax_file_dir)

    normal_minmax_fname = os.path.join(args.minmax_file_dir, 'minmax_normal-distance_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option+standardize_flag+num_chrom_flag+'.csv')
    stem_minmax_fname = os.path.join(args.minmax_file_dir, 'minmax_stem-distance_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option+standardize_flag+num_chrom_flag+'.csv')

    normal_distance_minmax = pd.read_csv(os.path.join(args.minmax_file_dir, normal_minmax_fname), index_col = 0)
    stem_distance_minmax = pd.read_csv(os.path.join(args.minmax_file_dir, stem_minmax_fname), index_col = 0)

    normal_distance_min = float(normal_distance_minmax.loc[args.cohort]['min'])
    normal_distance_max = float(normal_distance_minmax.loc[args.cohort]['max'])
    stem_distance_min = float(stem_distance_minmax.loc[args.cohort]['min'])
    stem_distance_max = float(stem_distance_minmax.loc[args.cohort]['max'])
 
    print("1. normal_distance and stem_cell_distance")    
    normal_distance_fname = os.path.join(cohort_dir, 'normal-distance_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option+standardize_flag+num_chrom_flag+'.csv')
    stem_distance_fname = os.path.join(cohort_dir, 'stem-distance_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option+standardize_flag+num_chrom_flag+'.csv')
    normal_distance = pd.read_csv(normal_distance_fname, index_col = 0).loc[S].copy().values.flatten()
    stem_distance = pd.read_csv(stem_distance_fname, index_col = 0).loc[S].copy().values.flatten()

    if args.normalize=='Y' and args.minmax=='N':
        if args.distance=='euclidean':
            print("Normalize normal- and stem- distances. X /= max(X).")
            normal_distance /= normal_distance_max
            stem_distance /= stem_distance_max
        else:
            pass

    if args.minmax=='Y' and args.normalize=='N': 
        print("Use minmax scaling. (X-min(X)) / (max(X)-min(X))")
        normal_distance = (normal_distance - normal_distance_min) / (normal_distance_max - normal_distance_min)
        stem_distance = (stem_distance - stem_distance_min) / (stem_distance_max - stem_distance_min)    
        for i in range(len(normal_distance)):
            if abs(normal_distance[i]) <= 1e-14:
                normal_distance[i] = 0
            assert normal_distance[i] >= 0
        for i in range(len(stem_distance)):
            if abs(stem_distance[i]) <= 1e-14:
                stem_distance[i] = 0
            assert stem_distance[i] >= 0

    if args.minmax=='N' and args.normalize=='N':
        print("Use raw score. No preprocessing needed.")
     
    print("2. compute angle")
    radian_theta, degree_theta, cosine_radian = compute_angle(normal_distance, stem_distance, args.cohort) 
    final_df = pd.DataFrame(zip(normal_distance, stem_distance, radian_theta, degree_theta, cosine_radian), 
                            index = S, columns = ['normal_distance', 'stem_distance', 'angle_radian', 'angle_degree','cos_radian'])
                            

    print("3. save result")
    normalize_yn = '_normalized' if args.normalize=='Y' else ''
    minmax_yn = '_minmax' if args.minmax=='Y' else ''
    
    result_fname = os.path.join(cohort_dir, 
        'stem-closeness_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option + normalize_yn+minmax_yn+standardize_flag+num_chrom_flag+'.csv') 
    final_df.to_csv(result_fname)
    print('result file: {}'.format(result_fname))
    print("----")

    plt.clf()
