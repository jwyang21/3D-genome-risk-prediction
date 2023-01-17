#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os
import sys
import argparse

# - DMR_GENE-ANNOT{}ALL.npy
# - DMR_REG-ANNOT{}_ALL.npy
# - DMR_EPI-ANNOT{}_ALL.npy

# In[3]:


# global variables
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
#THRESHOLD = ['mean', 'mean_std', 'mean_2std'] #use mean_std only. 
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
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    args.add_argument('-b', '--binsize', help = 'binsize', type = int, required = True, default = int(1e6))
    args.add_argument('-t', '--threshold', help = 'DMR threshold. mean, mean_std, mean_2std', type = str, required = True, default = 'mean_std')
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    return args.parse_args()


# To Do
# - GENE features (gene, transcript): 
#     - cohort 별로 
#         - chromosome 별로 
#             - 각 feature의 길이 summation
#     - cohort 별로, 모든 feature 모아서 저장 (ENSG ID or ENST ID)
# - REG features:
#     - cohort 별로
#         - chromosome 별로
#             - 각 feature의 길이 summation
#     - cohort 별로, 모든 feature 모아서 저장 (ENSR ID)
# - EPI features:
#     - cohort 별로, 각 feature의 길이 summation. 

# In[5]:


def parse_GENE(cohort_dir, gene_annot, threshold_key):
    cohort_ensg = []
    cohort_enst = []
    GENE = {}
    
    for chrom in CHR_LIST:
        globals()[chrom+'_ENSG_len'] = 0
        globals()[chrom+'_ENST_len'] = 0
    del(chrom)

    for f in gene_annot:

        if 'DMR_GENE-ANNOT_gene' in f and threshold_key in f:
            all_annot = np.load(os.path.join(cohort_dir, f))
            for annot in all_annot:
                current_id = annot.split('_')[-1].strip()
                if current_id not in cohort_ensg:
                    cohort_ensg.append(current_id)
                current_feature_start = int(annot.split(':')[1].split('-')[0])
                current_feature_end = int(annot.split('-')[1].split('_')[0])
                assert current_feature_end >= current_feature_start
                current_feature_length = (current_feature_end - current_feature_start)
                #current_feature_length = (current_feature_end - current_feature_start) + 1
                current_chrom = annot.split(':')[0].strip()
                globals()[current_chrom+'_ENSG_len'] += current_feature_length

        elif 'DMR_GENE-ANNOT_transcript' in f and threshold_key in f:
            all_annot = np.load(os.path.join(cohort_dir, f))
            for annot in all_annot:
                current_id = annot.split('_')[-1].strip()
                if current_id not in cohort_enst:
                    cohort_enst.append(current_id)
                current_feature_start = int(annot.split(':')[1].split('-')[0])
                current_feature_end = int(annot.split('-')[1].split('_')[0])
                assert current_feature_end >= current_feature_start
                current_feature_length = (current_feature_end - current_feature_start)
                #current_feature_length = (current_feature_end - current_feature_start) + 1

                current_chrom = annot.split(':')[0].strip()
                globals()[current_chrom+'_ENST_len'] += current_feature_length   
        else:
            pass

    for chrom in CHR_LIST:
        GENE[chrom+'_ENSG_len'] = globals()[chrom+'_ENSG_len']
        GENE[chrom+'_ENST_len'] = globals()[chrom+'_ENST_len']
    GENE['ENSG_all'] = cohort_ensg
    GENE['ENST_all'] = cohort_enst
    
    return GENE

def parse_REG(cohort_dir, reg_annot, threshold_key):
    cohort_ensr = []
    REG = {}
    for chrom in CHR_LIST:
        for small_c in SMALL_CATEGORY['REG']:
            globals()[chrom+'_'+small_c+'_len'] = 0
    del(chrom)
            
    for f in reg_annot:
        for small_c in SMALL_CATEGORY['REG']:
            if small_c in f and threshold_key in f:
                all_annot = np.load(os.path.join(cohort_dir, f))
                for annot in all_annot:
                    current_id = annot.split('_')[-1].strip()
                    if current_id not in cohort_ensr:
                        cohort_ensr.append(current_id)
                    current_feature_start = int(annot.split(':')[1].split('-')[0])
                    current_feature_end = int(annot.split('-')[1].split('_')[0])
                    assert current_feature_end >= current_feature_start
                    current_feature_length = (current_feature_end - current_feature_start)
                    #current_feature_length = (current_feature_end - current_feature_start) + 1
                    current_chrom = annot.split(':')[0].strip()
                    globals()[current_chrom+'_'+small_c+'_len'] += current_feature_length
    
        
    for chrom in CHR_LIST:
        for small_c in SMALL_CATEGORY['REG']:
            REG[chrom+'_'+small_c+'_len'] = globals()[chrom+'_'+small_c+'_len']
    REG['ENSR_all'] = cohort_ensr
    return REG

def parse_EPI(cohort_dir, epi_annot, threshold_key, ):
    EPI = {}
    
    for chrom in CHR_LIST:
        for small_c in SMALL_CATEGORY['EPI']:
            if '/' in small_c:
                small_c = small_c.replace('/', '-')
            globals()[chrom+'_'+small_c+'_len'] = 0
    del(chrom)    
    
    for f in epi_annot:
        for small_c in SMALL_CATEGORY['EPI']:
            if '/' in small_c:
                small_c = small_c.replace('/', '-')
            if small_c in f and threshold_key in f:
                all_annot = np.load(os.path.join(cohort_dir, f))
                for annot in all_annot:
                    current_feature_start = int(annot.split(':')[1].split('-')[0])
                    current_feature_end = int(annot.split('-')[1].split('_')[0])
                    assert current_feature_end >= current_feature_start
                    current_feature_length = (current_feature_end - current_feature_start)
                    #current_feature_length = (current_feature_end - current_feature_start) + 1
                    current_chrom = annot.split(':')[0].strip()
                    globals()[current_chrom+'_'+small_c+'_len'] += current_feature_length
       
    for chrom in CHR_LIST:
        for small_c in SMALL_CATEGORY['EPI']:
            if '/' in small_c:
                small_c = small_c.replace('/', '-')
            EPI[chrom+'_'+small_c+'_len'] = globals()[chrom+'_'+small_c+'_len']
    return EPI


if __name__ == '__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    cohort_dir = os.path.join(os.getcwd(), 'result', args.cohort)#이 cohrot에 대한 분석 결과 파일들이 있는 디렉토리
    threshold_key = 'threshold_'+args.threshold+'_binsize'
    print("cohort_directory: {}".format(cohort_dir))
    
    print("---\n0.list up annotation files per big category (GENE, REG, EPI)")
    gene_annot = []
    reg_annot = []
    epi_annot = []
    
    for f in os.listdir(cohort_dir):
        if f.startswith('DMR_GENE-ANNOT') and f.endswith('ALL.npy'):
            gene_annot.append(f)
        elif f.startswith('DMR_REG-ANNOT') and f.endswith('ALL.npy'):
            reg_annot.append(f)
        elif f.startswith('DMR_EPI-ANNOT') and f.endswith('ALL.npy'):
            epi_annot.append(f)
        else:
            pass
    
    print("---\n1. for GENE, collect all ENSG or ENST IDs")
    GENE = parse_GENE(cohort_dir, gene_annot, threshold_key)   
    np.savez(os.path.join(cohort_dir, 'DMR_GENE_features_threshold_'+args.threshold+'_id_len'), **GENE)    
    print("result file: {}".format(os.path.join(cohort_dir, 'DMR_GENE_features_threshold_'+args.threshold+'_id_len')+'.npz'))

    print("---\n2. for REG, collect items belonging to each category")
    REG = parse_REG(cohort_dir, reg_annot, threshold_key)
    np.savez(os.path.join(cohort_dir, 'DMR_REG_features_threshold_'+args.threshold+'_id_len'), **REG)   
    print("result file: {}".format(os.path.join(cohort_dir, 'DMR_REG_features_threshold_'+args.threshold+'_id_len')+'.npz'))
 
    print("---\n3. for EPI, compute length of each feature.")
    EPI = parse_EPI(cohort_dir, epi_annot, threshold_key)
    np.savez(os.path.join(cohort_dir, 'DMR_EPI_features_threshold_'+args.threshold+'_len'), **EPI)  
    print("result file: {}".format(os.path.join(cohort_dir, 'DMR_EPI_features_threshold_'+args.threshold+'_len')+'.npz'))
 
