#!/usr/bin/env python
# coding: utf-8
'''
- Parse DMR feature results (for all 22 autosomes)
- Collect DMR features per cohort and per threshold, and write all of them into an output file. 
- DMR feature fname: /data/project/3dith/pipelines/island-pipeline/3_dmr-island/result/{cohort}/DMR_{big_category}-ANNOT_{small_category}_threshold_{threshold}_binsize_{binsize}.pickle
'''
# In[1]:


import os
import pandas as pd
import numpy as np
import pickle
import sys
import argparse


# global variables
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
THRESHOLD = ['mean', 'mean_std', 'mean_2std']
SMALL_CATEGORY_FNAME = '/data/project/3dith/data/etc/dmr-feature-small-category.npz'
if os.path.exists(SMALL_CATEGORY_FNAME):
    SMALL_CATEGORY = np.load(SMALL_CATEGORY_FNAME, allow_pickle = True)
else:
    SMALL_CATEGORY = {}
    SMALL_CATEGORY['GENE'] = ['gene','transcript']
    SMALL_CATEGORY['REG'] = ['open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', 'enhancer', 'promoter',                             'promoter_flanking_region']
    SMALL_CATEGORY['EPI'] = ['18_Quies', '13_Het', '17_ReprPCWk', '16_ReprPC', '14_TssBiv', '2_TssFlnk',     '12_ZNF/Rpts', '11_EnhWk', '1_TssA', '6_TxWk', '5_Tx', '9_EnhA1', '7_EnhG1',     '4_TssFlnkD' ,'15_EnhBiv', '10_EnhA2', '3_TssFlnkU', '8_EnhG2']
    np.savez(SMALL_CATEGORY_FNAME,**SMALL_CATEGORY)
    
cohort2eid = pd.read_csv('/data/project/3dith/data/etc/cohort2eid.txt', sep = '\t', header = None)
cohort2eid.columns = ['cohort', 'eid']
eid_cohorts = cohort2eid.cohort.values


# In[16]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    parser.add_argument('-b', '--binsize', help = 'binsize. default 1e6', type = int, default = int(1e6), required = False)
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    #default w_dir: /data/project/3dith/pipelines/island-pipeline/3_dmr-island
    return parser.parse_args()


if __name__=='__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    result_dir = os.path.join(os.getcwd(), 'result')
    cohort_dir = os.path.join(os.getcwd(), 'result', args.cohort)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)

    print("===\ncohort: {}".format(args.cohort))
    print("cohort_dir: {}".format(cohort_dir))

    for big_c in BIG_CATEGORY:
        if big_c=='EPI' and args.cohort not in eid_cohorts:
            print("{} has no matching EID.".format(args.cohort))
            continue
        #print("---\nbig_category: ", big_c)
        small_category = SMALL_CATEGORY[big_c]
        for small_c in small_category:
            if '/' in small_c:
                small_c = small_c.replace('/', '-')
            #print("small_cateogry: ", small_c)
            for t in THRESHOLD:
                if t == 'mean_std': # opensea 실험결과, mean_std일 때가 DMR에 있는 gene 개수가 제일 적절했으므로 이때만 run.
                    #print(t)
                    #print("threshold: ", t)
                    # 읽어야 할 정보가 들어있는 파일
                    if big_c != 'EPI':
                        dmr_feature_fname = os.path.join(cohort_dir, 'DMR_'+big_c+'-ANNOT_'+small_c+'_threshold_'+t+'_binsize_'+str(args.binsize)+'.pickle')
                    else:
                        dmr_feature_fname = os.path.join(cohort_dir, big_c+'-ANNOT_'+small_c+'_threshold_'+t+'_binsize_'+str(args.binsize)+'.pickle')
                    if os.path.exists(dmr_feature_fname):
                        with open(dmr_feature_fname, 'rb') as f:
                            dmr_feature = pickle.load(f)
                        f.close()
                        #print("Parsing {}".format(dmr_feature_fname))

                        # 읽어들인 정보를 쓸 파일
                        dmr_feature_parse_fname = os.path.join(cohort_dir, 'DMR_'+big_c+'-ANNOT_'+small_c+'_threshold_'+t+'_binsize_'+str(args.binsize)+'_ALL')#.npy')

                        dmr_feature_parse_all = []
                        cnt = 0
                        for k in list(dmr_feature.keys()):
                            current_value = dmr_feature[k].copy()
                            cnt += len(current_value)
                            for v in current_value:
                                dmr_feature_parse_all.append(v)
                        assert cnt == len(dmr_feature_parse_all)
                        if len(dmr_feature_parse_all) > 0:
                            dmr_feature_parse_all_array = np.array(dmr_feature_parse_all)
                            np.save(dmr_feature_parse_fname, dmr_feature_parse_all_array)
                            print(dmr_feature_parse_fname+'.npy')
                        else:
                            print("no DMR feature parsed for {}/{}/{}.".format(big_c, small_c, t))
                        #print("Parsed data is saved in {}".format(dmr_feature_parse_fname))

