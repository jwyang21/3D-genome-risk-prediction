#compute-distance
#!/usr/bin/env python
# coding: utf-8
#compute-distance
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

# global variables
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}
CPG_TYPE='shelf_shore'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--score_type', help = 'pc1-avg or pc1-fluctuation.', type = str, default = 'pc1-avg', required = True)
    #default: pc1-avg
    
    parser.add_argument('--usage_option', help = 'use all samples or randomly picked samples when computing reference.', type = str, default = 'part', required = True) # all, part
    #default: part
    
    parser.add_argument('-w_dir','--working_dir', default = 'working directory', type = str, required = True)
    #default: /data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}

    parser.add_argument('--pc1_upper_dir', help = 'directory where pc1 directories are in', type = str, required = True)
    #default: '/data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result'
    
    parser.add_argument('--result_dir', help = 'result directory', type = str, required = True)
    #default: '/data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result' 

    parser.add_argument('--reference_dir', help = 'normal/stem reference vector directory', type = str, required = True)#reference vector directory
    #default: /data/project/3dith/pipelines/{CPG_TYPE}-pipeline/1_compute-score-{CPG_TYPE}/result #/{cohort}

    parser.add_argument('--hg19_chr_len', help = 'hg19 chromosome length file', type = str, required = True)
    #default: /data/project/3dith/data/hg19.fa.sizes'

    parser.add_argument('--bdm_bin_dir', help = 'bdm bin directory', type = str, required = True)
    #/data/project/3dith/data/bdm_bins/{CPG_TYPE}/ #TCGA-COAD_diffmat_bins.npz

    parser.add_argument('--use_weighted_avg', help = 'whether to use weighted average in computing distance. Y or N (simple_avg)', type = str, required = True)
    # try both Y and N

    parser.add_argument('--distance', help = 'distance metric. cosine-sim or euclidean', type = str, required = True)
    # try both 'euclidean' and 'cosine-sim'

    parser.add_argument('--pseudocount', help = 'pseudocount to be used when args.distance==cosine-sim to prevent division by 0', type = float, default = 1e-15, required = False)
    #default: 1e-15

    parser.add_argument('-m_type', '--matrix_type', help = 'from which matrix you will compute PC1. bdm, iebdm, or all.', type = str, required = True)
    #default: 'iebdm' #try both iebdm and bdm.

    parser.add_argument('--cpg_type', help = 'cpg type. island, shelf, shore, shelf_shore, opensea', type = str, required = True)

    parser.add_argument('--cohort', help = 'cohort. TCGA-{} or PCBC', type = str, required = True)

    parser.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
    
    parser.add_argument('--num_chrom', help = 'chr1 to which chromosome you will use.', type = int, default = 22, required = False)
    
    return parser.parse_args()

def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = [] 
    else: #TCGA cohort
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S

def cosine_sim(v1, v2):
    # v1 and v2 should be individual vectors of same length. 
    # referred to https://wikidocs.net/24603
    assert len(v1) == len(v2)
    dot_product_ = dot(v1, v2)
    cosine_sim_ = dot_product_ / (norm(v1) * norm(v2))
    return cosine_sim_

def standardize(v):
    return (v - v.mean()) / v.std()

def import_pc1(pc1_dir, sample, chrom, standardize_flag, flag='inv'):
    # import pre-computed PC1 of sample-of-interest
    # flag: 'raw' or 'inv'
    if flag=='raw':
        #pass #do not consider this case #pc1 from BDM. 
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

'''
def compute_pc1_fluctuation(pc1): #pc1 절댓값 적분
    pc1_abs = np.array([abs(x) for x in pc1])
    abs_area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return abs_area 
'''
def compute_distance_pc1_avg(S, pc1_dir, reference, distance, pseudocount, use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, reference_bdm_bins, pc1_flag, standardize_flag, CHR_LIST):
    if cohort_ref_diff_flag=='Y':
        cohort_bins = {}
        reference_bins = {}#reference: pcbc
        intersecting_bins = {}
        cohort_bin_mask = {}
        reference_bin_mask = {}

        for chrom in CHR_LIST:
            cohort_bins[chrom] = bdm_bins[chrom+'_bins']
            reference_bins[chrom] = reference_bdm_bins[chrom+'_bins']
            intersecting_bins[chrom] = np.intersect1d(cohort_bins[chrom], reference_bins[chrom])
            cohort_bin_mask[chrom] = [True if x in intersecting_bins[chrom] else False for x in cohort_bins[chrom]]
            reference_bin_mask[chrom] = [True if x in intersecting_bins[chrom] else False for x in reference_bins[chrom]]
    for s in S: #iterate for all samples
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
            else: #cosine similarity
                similarity = 1/cosine_sim(cohort_pc1, reference_pc1)
                dist[chrom] = 1/(similarity + pseudocount)
                # higher similarity means smaller distance between two vectors. 
        total_dist = 0
        for chrom in CHR_LIST:
            if use_weighted_avg=='Y':
                total_dist += dist[chrom] * chrom_weight[chrom]
            else:
                total_dist += ((dist[chrom]) * (1/(len(CHR_LIST))))
        distance_df.loc[s] = total_dist

    return distance_df
'''
def compute_distance_pc1_fluctuation(S, pc1_dir, reference, distance, pseudocount, use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, reference_bdm_bins,pc1_flag, standardize_flag):
    if cohort_ref_diff_flag=='Y':#cohort는 tcga, reference는 pcbc
        cohort_bins = {}
        reference_bins = {}
        cohort_weight = {}
        reference_weight = {}
        for chrom in CHR_LIST:
            #initialize weights as 1
            cohort_weight[chrom] = 1
            reference_weight[chrom] = 1
            cohort_bins[chrom] = bdm_bins[chrom+'_bins']
            reference_bins[chrom] = reference_bdm_bins[chrom+'_bins']

            if len(cohort_bins[chrom]) != len(reference_bins[chrom]):
                if len(reference_bins[chrom]) < len(cohort_bins[chrom]):
                    #길이가 더 작은 쪽에 맞춰서 normalize. 
                    cohort_weight[chrom] = len(reference_bins[chrom]) / len(cohort_bins[chrom])
                else: # len(reference_bins[chrom]) > len(cohort_bins[chrom])
                    reference_weight[chrom] = len(cohort_bins[chrom]) / len(reference_bins[chrom])

    for s in S: #iterate for all samples
        cohort_pc1_fluctuation_v = {}
        reference_pc1_fluctuation_v = {}
        
        for chrom in CHR_LIST:
            cohort_pc1 = import_pc1(pc1_dir, s, chrom, pc1_flag)
            
            #reference_pc1_fluctuation = reference[chrom]
            reference_pc1 = reference[chrom]
            
            if standardize_flag=='Y':
                cohort_pc1 = standardize(cohort_pc1)
                reference_pc1 = standardize(reference_pc1)

            cohort_pc1_fluctuation = compute_pc1_fluctuation(cohort_pc1)
            reference_pc1_fluctuation = compute_pc1_fluctuation(reference_pc1)
            if cohort_ref_diff_flag=='Y':
                cohort_pc1_fluctuation *= cohort_weight[chrom]
                reference_pc1_fluctuation *= reference_weight[chrom]
            cohort_pc1_fluctuation_v[chrom] = cohort_pc1_fluctuation
            reference_pc1_fluctuation_v[chrom]=reference_pc1_fluctuation

        if distance=='euclidean':
            dist = {}
            total_dist = 0
            for chrom in CHR_LIST:
                if use_weighted_avg=='Y':
                    dist[chrom] = pow((cohort_pc1_fluctuation_v[chrom] - reference_pc1_fluctuation_v[chrom]),2) * chrom_weight[chrom]
                else:
                    dist[chrom] = pow((cohort_pc1_fluctuation_v[chrom] - reference_pc1_fluctuation_v[chrom]),2) * (1/len(CHR_LIST))
            for chrom in CHR_LIST:
                total_dist += dist[chrom] 
            total_dist = np.sqrt(total_dist)

        else: #disatnce=='cosine similarity'
            cohort_pc1_fluctuation_vector = []
            reference_pc1_fluctuation_vector = []
            for chrom in CHR_LIST:
                cohort_pc1_fluctuation_vector.append(cohort_pc1_fluctuation_v[chrom])
                reference_pc1_fluctuation_vector.append(reference_pc1_fluctuation_v[chrom])
            similarity = 1/cosine_sim(cohort_pc1_fluctuation_vector, reference_pc1_fluctuation_vector)
            total_dist = 1/(similarity + pseudocount)
            # higher similarity means smaller distance between two vectors. 
        
        distance_df.loc[s] = total_dist

    return distance_df
'''
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

        # define output filename
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

        if current_reference != args.cohort: #stem_distance # 이런 경우는 cohort가 TCGA, current_reference가 PCBC인 경우밖에 없음 
            print("---\ncohort {} != reference {}".format(args.cohort, current_reference))
            assert current_reference=='PCBC'
            reference_fname = os.path.join(pcbc_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            reference = np.load(reference_fname)
            pcbc_bdm_bins = np.load(os.path.join(args.bdm_bin_dir,'PCBC_diffmat_bins.npz'))
            bdm_bins = np.load(os.path.join(args.bdm_bin_dir, args.cohort+'_diffmat_bins.npz'))
            if args.score_type == 'pc1-avg':
                # 서로 다른 두 코호트 간에는 binned diffmat에 쓰인 bin들의 목록이 차이 나므로, intersecting bins를 먼저 찾아야 함.
                distance_df = compute_distance_pc1_avg(S, pc1_dir, reference, args.distance, args.pseudocount, 
                    args.use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, pcbc_bdm_bins, pc1_flag, args.standardize, CHR_LIST)
            else: #args.score_type == 'pc1-fluctuation'
                pass
                '''
                distance_df = compute_distance_pc1_fluctuation(S, pc1_dir, reference, args.distance, args.pseudocount, 
                    args.use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, bdm_bins, pcbc_bdm_bins, pc1_flag)
                '''
            distance_df.to_csv(distance_df_fullname)
            print("{} filename: {}".format(distance_type, distance_df_fullname))
                
        else: #current_reference == cohort #TCGA cohort의 normal distance 또는 PCBc 내에서의 stem distance. 
            print("---\ncohort {} == reference {}".format(args.cohort, current_reference))
            if 'TCGA' in args.cohort:
                reference_fname = os.path.join(cohort_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            elif args.cohort=='PCBC':
                assert cohort_dir == pcbc_dir
                reference_fname = os.path.join(cohort_dir, reference_type_flag+'_'+args.matrix_type+'_'+args.score_type+'_'+args.usage_option+standardize_flag_yn+'.npz')
            else:
                raise Exception("wrong cohort. neither TCGA nor PCBC")

            reference = np.load(reference_fname)    
            if args.score_type == 'pc1-avg': #TCGA cohort 내에서 normal distance를 계산하거나, pcbc 내에서 stem distance를 계산하거나. 
                distance_df = compute_distance_pc1_avg(S, pc1_dir, reference, args.distance, args.pseudocount, args.use_weighted_avg, distance_df, 
                    chrom_weight, cohort_ref_diff_flag, [], [], pc1_flag, args.standardize, CHR_LIST)
            else: #args.score_type == 'pc1-fluctuation'
                pass
                '''
                distance_df = compute_distance_pc1_fluctuation(S, pc1_dir, reference, args.distance, args.pseudocount, 
                    args.use_weighted_avg, distance_df, chrom_weight, cohort_ref_diff_flag, [],[],pc1_flag)
                '''
            distance_df.to_csv(distance_df_fullname)
            print("{} filename: {}".format(distance_type, distance_df_fullname))
