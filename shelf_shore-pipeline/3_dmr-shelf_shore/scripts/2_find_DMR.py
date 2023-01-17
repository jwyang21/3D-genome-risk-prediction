#!/usr/bin/env python
# coding: utf-8

# find DMR based on the specified binsize.(default: 1Mbp)

# In[1]:


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
    
    #initialize
    #mean_df = pd.DataFrame(np.zeros((len(CHR_LIST), len(NORMAL7_COHORT)), dtype = float), index = CHR_LIST, columns = NORMAL7_COHORT)
    mean_std_df = pd.DataFrame(np.zeros((len(CHR_LIST), len(NORMAL7_COHORT)), dtype = float), index = CHR_LIST, columns = NORMAL7_COHORT)
    #mean_2std_df = pd.DataFrame(np.zeros((len(CHR_LIST), len(NORMAL7_COHORT)), dtype = float), index = CHR_LIST, columns = NORMAL7_COHORT)

    for cohort in NORMAL7_COHORT:#real
    #for cohort in NORMAL7_COHORT[:1]:#test#TCGA-BLCA
        print("===\ncohort: ", cohort)
        cohort_dir = os.path.join(result_dir, cohort)
        print("cohort_dir: {}".format(cohort_dir))

        #initialize    
        globals()[cohort+'_DMR'] = {} #save mean, mean-std, mean-2std for each chromosome. 

        #if saving directory does not exist, make directory
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        if not os.path.exists(cohort_dir):
            os.makedirs(cohort_dir)

        # load result 
        fname = os.path.join(cohort_dir, 'binned_avg_'+args.cpg_type+'_TN_diff_binsize_'+str(args.binsize)+'.npz')
        f = np.load(fname, allow_pickle = True)

        for chrom in CHR_LIST:#real
        #for chrom in CHR_LIST[:1]:#test
        #for chrom in CHR_LIST[13:14]:#test#chr14
            print("---\n"+chrom)
            v = f[chrom].copy()
            b = f[chrom+'_bins'].copy()
            current_chrom_df = pd.DataFrame(v, index = b, columns = ['TN_diff'])
            current_chrom_df.dropna(inplace = True)

            current_mean = np.mean(current_chrom_df.TN_diff.values)
            current_std = np.std(current_chrom_df.TN_diff.values)

            #save thresholds
            #mean_df.loc[chrom][cohort] = current_mean.copy()
            mean_std_df.loc[chrom][cohort] = (current_mean - current_std).copy()
            #mean_2std_df.loc[chrom][cohort] = (current_mean - (2*current_std)).copy()

            #save thresholded result (DMR bins)
            #mean_mask = current_chrom_df.TN_diff.values < current_mean
            mean_std_mask = current_chrom_df.TN_diff.values < (current_mean - current_std)
            #mean_2std_mask = current_chrom_df.TN_diff.values < (current_mean - (2*current_std))

            #save DMR bins
            #globals()[cohort+'_DMR'][chrom+'_mean_bins'] = []#[x if x != '' for x in ]#current_chrom_df.index.values * mean_mask
            globals()[cohort+'_DMR'][chrom+'_mean_std_bins'] = []#[x if x != '' for x in ]#current_chrom_df.index.values * mean_std_mask
            #globals()[cohort+'_DMR'][chrom+'_mean_2std_bins'] = []#[x if x != '' for x in ]#current_chrom_df.index.values * mean_2std_mask

            #mean_bins = (current_chrom_df.index.values * mean_mask)
            mean_std_bins = (current_chrom_df.index.values * mean_std_mask)
            #mean_2std_bins = (current_chrom_df.index.values * mean_2std_mask)
    
            '''
            mean_bins_cnt = 0
            for x in mean_bins:
                if x != '':
                    globals()[cohort+'_DMR'][chrom+'_mean_bins'].append(x)
                    mean_bins_cnt += 1
            '''
            mean_std_bins_cnt = 0
            for x in mean_std_bins:
                if x != '':
                    globals()[cohort+'_DMR'][chrom+'_mean_std_bins'].append(x)
                    mean_std_bins_cnt += 1
            '''
            mean_2std_bins_cnt = 0
            for x in mean_2std_bins:
                if x != '':
                    globals()[cohort+'_DMR'][chrom+'_mean_2std_bins'].append(x)
                    mean_2std_bins_cnt += 1
            '''
            #assert mean_bins_cnt == len(globals()[cohort+'_DMR'][chrom+'_mean_bins'])
            assert mean_std_bins_cnt == len(globals()[cohort+'_DMR'][chrom+'_mean_std_bins'])
            #assert mean_2std_bins_cnt == len(globals()[cohort+'_DMR'][chrom+'_mean_2std_bins'])

            #print("proportion_mean_bins: {}".format(mean_bins_cnt / current_chrom_df.shape[0]))
            print("proportion_mean_std_bins: {}".format(mean_std_bins_cnt / current_chrom_df.shape[0]))
            #print("proportion_mean_2std_bins: {}".format(mean_2std_bins_cnt / current_chrom_df.shape[0]))

        # save results of current cohort #DMR bins
        print("---\nSave DMR bins of current cohort (per each chrom).")
        result_fname = os.path.join(cohort_dir, 'DMR_binsize_'+str(args.binsize))
        print("DMR result file: {}".format(result_fname+'.npz'))
        np.savez(result_fname, **globals()[cohort+'_DMR'])
        
    #save result of all cohorts #mean, mean-std, mean-2std values per chromosome, per cohort.
    print("===\nSave result of all cohorts (shape: num_chrom, num_cohorts).")
    #print(os.path.join(result_dir, 'chrom_cohort_mean_TN_diff_binsize_'+str(args.binsize)+'.csv'))
    print(os.path.join(result_dir, 'chrom_cohort_mean_std_TN_diff_binsize_'+str(args.binsize)+'.csv'))
    #print(os.path.join(result_dir, 'chrom_cohort_mean_2std_TN_diff_binsize_'+str(args.binsize)+'.csv'))
    #mean_df.to_csv(os.path.join(result_dir, 'chrom_cohort_mean_TN_diff_binsize_'+str(args.binsize)+'.csv'), index = True)
    mean_std_df.to_csv(os.path.join(result_dir, 'chrom_cohort_mean_std_TN_diff_binsize_'+str(args.binsize)+'.csv'), index = True)
    #mean_2std_df.to_csv(os.path.join(result_dir, 'chrom_cohort_mean_2std_TN_diff_binsize_'+str(args.binsize)+'.csv'), index = True)


