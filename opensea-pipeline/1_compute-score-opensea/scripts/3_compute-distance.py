import numpy as np
import pandas as pd
import os
import sys
from scipy.spatial.distance import euclidean
import argparse
from scipy.stats import ttest_ind
from numpy import dot
from numpy.linalg import norm
from scipy.integrate import simps
from numpy import trapz

SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='opensea'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--score_type', help = 'pc1-avg or pc1-fluctuation.', type = str, default = 'pc1-avg', required = True)
    parser.add_argument('--usage_option', help = 'use all samples or randomly picked samples when computing reference.', type = str, default = 'half', required = True) # all, part
    parser.add_argument('-w_dir','--working_dir', default = 'working directory', type = str, required = True)
    parser.add_argument('--pc1_upper_dir', help = 'directory where pc1 directories are in', type = str, required = True)
    parser.add_argument('--result_dir', help = 'result directory', type = str, required = True)
    parser.add_argument('--reference_dir', help = 'normal/stem reference vector directory', type = str, required = True)#reference vector directory
    parser.add_argument('--hg19_chr_len', help = 'hg19 chromosome length file', type = str, required = True)
    parser.add_argument('--bdm_bin_dir', help = 'bdm bin directory', type = str, required = True)
    parser.add_argument('--use_weighted_avg', help = 'whether to use weighted average in computing distance. Y or N (simple_avg)', type = str, required = True)
    parser.add_argument('--distance', help = 'distance metric. cosine-sim or euclidean', type = str, required = True)
    parser.add_argument('--pseudocount', help = 'pseudocount to be used when args.distance==cosine-sim to prevent division by 0', type = float, default = 1e-15, required = False)
    parser.add_argument('-m_type', '--matrix_type', help = 'from which matrix you will compute PC1. bdm, iebdm, or all.', type = str, required = True)
    parser.add_argument('--cpg_type', help = 'cpg type. island, shelf, shore, shelf_shore, opensea', type = str, required = True)
    parser.add_argument('--cohort', help = 'cohort. TCGA-{} or PCBC', type = str, required = True)
    parser.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
    parser.add_argument('--num_chrom', help = 'chr1 to which chromosome you will use.', type = int, default = 22, required = False)
    
    return parser.parse_args()

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

def cosine_sim(v1, v2):
    # referred to https://wikidocs.net/24603
    assert len(v1) == len(v2)
    dot_product_ = dot(v1, v2)
    cosine_sim_ = dot_product_ / (norm(v1) * norm(v2))
    return cosine_sim_

def standardize(v):
    return (v - v.mean()) / v.std()

def import_pc1(pc1_dir, sample, chrom, standardize_flag, flag='inv'):
    if flag=='raw':
        fname = os.path.join(pc1_dir, sample+'.npz')
    elif flag=='inv':
        fname = os.path.join(pc1_dir, sample+'_inv_exp.npz')
    else:
        pass
    key_ = chrom+'_pc1'
    pc1 = np.load(fname)[key_]
    if standardize_flag=='Y':
        pc1 = standardize(pc1)
    return pc1


def compute_distance_pc1_avg(S, pc1_dir, reference, distance, pseudocount, use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, reference_bdm_bins, pc1_flag, standardize_flag, CHR_LIST):
    if cohort_ref_diff_flag=='Y':
        cohort_bins = {}
        reference_bins = {}
        intersecting_bins = {}
        cohort_bin_mask = {}
        reference_bin_mask = {}

        for chrom in CHR_LIST:
            cohort_bins[chrom] = bdm_bins[chrom+'_bins']
            reference_bins[chrom] = reference_bdm_bins[chrom+'_bins']
            intersecting_bins[chrom] = np.intersect1d(cohort_bins[chrom], reference_bins[chrom])
            cohort_bin_mask[chrom] = [True if x in intersecting_bins[chrom] else False for x in cohort_bins[chrom]]
            reference_bin_mask[chrom] = [True if x in intersecting_bins[chrom] else False for x in reference_bins[chrom]]
    for s in S: 
        dist = {}
        for chrom in CHR_LIST:
            cohort_pc1 = import_pc1(pc1_dir, s, chrom, standardize_flag, pc1_flag)
            reference_pc1 = reference[chrom]

            if cohort_ref_diff_flag=='Y':
                assert len(cohort_pc1) == len(cohort_bin_mask[chrom])
                assert len(reference_pc1) == len(reference_bin_mask[chrom])
                cohort_pc1 = pd.DataFrame(cohort_pc1).iloc[cohort_bin_mask[chrom],:].values.flatten()
                reference_pc1 = pd.DataFrame(reference_pc1).iloc[reference_bin_mask[chrom],:].values.flatten()
            if distance=='euclidean':
                dist[chrom] = euclidean(cohort_pc1, reference_pc1)
            else: 
                similarity = 1/cosine_sim(cohort_pc1, reference_pc1)
                dist[chrom] = 1/(similarity + pseudocount)
              
        total_dist = 0
        for chrom in CHR_LIST:
            if use_weighted_avg=='Y':
                total_dist += dist[chrom] * chrom_weight[chrom]
            else:
                total_dist += ((dist[chrom]) * (1/(len(CHR_LIST))))
        distance_df.loc[s] = total_dist

    return distance_df

if __name__ == '__main__':
    args = parse_arguments()

    weighted_avg_flag = '_weighted-avg' if args.use_weighted_avg=='Y' else '_simple-avg'
    distance_metric_flag = '_cosine-sim' if args.distance=='cosine-sim' else '_euclidean'
    matrix_type_flag = '_iebdm' if args.matrix_type=='iebdm' else '_bdm'
    pcbc_dir = os.path.join(args.result_dir, 'PCBC')
    pc1_flag = 'inv' if args.matrix_type=='iebdm' else 'raw'
    standardize_flag_yn = '_standardized' if args.standardize=='Y' else ''
    num_chrom_flag = '_to-chrom-'+str(args.num_chrom) if args.num_chrom < 22 else ''
    
    CHR_LIST = ['chr'+str(i) for i in np.arange(1, args.num_chrom+1)]

    if not os.path.exists(args.bdm_bin_dir):
        bdm_bin_dir_parent = os.path.abspath(os.path.join(args.bdm_bin_dir, os.pardir))
        if not os.path.exists(bdm_bin_dir_parent):
            os.makedirs(bdm_bin_dir_parent)
        os.makedirs(args.bdm_bin_dir)
        command_ = 'python3 /data/project/3dith/pipelines/utils/scripts/1_find-bdm-bins.py \
            -w_dir /data/project/3dith/pipelines/utils --tcga_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-v2-{}/result \
            --pcbc_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-pcbc/result \
            --hg19_chr_len /data/project/3dith/data/hg19.fa.sizes --binsize 1000000 --save_dir /data/project/3dith/data/bdm_bins \
            --cpg_type {}'.format(args.cpg_type, args.cpg_type)
        print("{} does not exists.".format(args.bdm_bin_dir))
        print("command: {}".format(command_))
        os.system(command_)

    chrom_weight = {}
    if args.use_weighted_avg=='Y':
        chr_len = pd.read_csv(args.hg19_chr_len, header = None, index_col = 0, sep = '\t')
        for chrom in CHR_LIST:
            chrom_weight[chrom] = chr_len.loc[chrom].values[0]/chr_len.values.sum()

    reference_type = ['PCBC']
 
    if 'TCGA' in args.cohort:
        reference_type.append(args.cohort)

    T, N, S = get_sample_list(args.cohort)

    print("===\ncohort: {}".format(args.cohort))
    cohort_dir = os.path.join(args.result_dir, args.cohort)
    pc1_dir = os.path.join(args.pc1_upper_dir, args.cohort, 'pc1')

    if not os.path.exists(args.result_dir):
        os.path.makedirs(args.result_dir)
    if not os.path.exists(cohort_dir):
        os.path.makedirs(cohort_dir)

    for i in range(len(reference_type)):
        current_reference = reference_type[i]
        print("---\ncurrent reference: ", current_reference)

       
        if args.cohort == current_reference:
            distance_type = 'normal-distance' if 'TCGA' in args.cohort else 'stem-distance'
            cohort_ref_diff_flag = 'N'
        else:
            distance_type = 'stem-distance'
            cohort_ref_diff_flag = 'Y'

        distance_df_fname = distance_type+distance_metric_flag+matrix_type_flag+'_'+args.score_type+weighted_avg_flag+'_'+args.usage_option+standardize_flag_yn+num_chrom_flag+'.csv'
        distance_df_fullname = os.path.join(cohort_dir, distance_df_fname)
        distance_df = pd.DataFrame(np.zeros((len(S), 1), dtype = float), index = S, columns = [distance_type])
        reference_type_flag = 'normal-reference' if distance_type=='normal-distance' else 'stem-reference'

        if current_reference != args.cohort: 
            print("---\ncohort {} != reference {}".format(args.cohort, current_reference))
            assert current_reference=='PCBC'
            reference_fname = os.path.join(pcbc_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            reference = np.load(reference_fname)
            pcbc_bdm_bins = np.load(os.path.join(args.bdm_bin_dir,'PCBC_diffmat_bins.npz'))
            bdm_bins = np.load(os.path.join(args.bdm_bin_dir, args.cohort+'_diffmat_bins.npz'))
            if args.score_type == 'pc1-avg':
                
                distance_df = compute_distance_pc1_avg(S, pc1_dir, reference, args.distance, args.pseudocount, 
                    args.use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, pcbc_bdm_bins, pc1_flag, args.standardize, CHR_LIST)
            else: 
                pass
            distance_df.to_csv(distance_df_fullname)
            print("{} filename: {}".format(distance_type, distance_df_fullname))
                
        else:  
            print("---\ncohort {} == reference {}".format(args.cohort, current_reference))
            if 'TCGA' in args.cohort:
                reference_fname = os.path.join(cohort_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            elif args.cohort=='PCBC':
                assert cohort_dir == pcbc_dir
                reference_fname = os.path.join(cohort_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            else:
                raise Exception("wrong cohort. neither TCGA nor PCBC")

            reference = np.load(reference_fname)    
            if args.score_type == 'pc1-avg': 
                distance_df = compute_distance_pc1_avg(S, pc1_dir, reference, args.distance, args.pseudocount, args.use_weighted_avg, distance_df, 
                    chrom_weight, cohort_ref_diff_flag, [], [], pc1_flag, args.standardize, CHR_LIST)
            else: 
                pass
            distance_df.to_csv(distance_df_fullname)
            print("{} filename: {}".format(distance_type, distance_df_fullname))
