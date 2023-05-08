#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
import pickle
import argparse


# In[ ]:


# # pipeline 
#  - 각 cohort마다
#      - 각 chromosome마다: globals()[chrom] = {}, globals()[chrom+'gene/reg/chromatin'] = []
#          - 각 DMR bin마다
#              - 이 bin 위에 있는 gene 찾기 -> append to globals()[chrom+'gene']
#              - 이 bin 위에 있는 regulatory feature 찾기. -> append to globals()[chrom+'reg']
#              - 이 bin 위에 있는 chromatin state 찾기.-> append to globals()[chrom+'chromatin']

# In[3]:


#CHROM_ANNOT = pd.read_csv('')
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
GENE_ANNOT = pd.read_csv('/data/project/3dith/data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf', sep = '\t', header = None, skiprows = [0, 1, 2, 3, 4])
# GENE_ANNOT column 1, 2, 3, 4: chrom, category, start, end #chrom format: 'chr{i}'
REG_ANNOT = pd.read_csv('/data/project/3dith/data/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff', sep = '\t', header = None)
# REG_ANNOT columns 1, 2, 3, 4: chrom, category, start, end #chrom format: 'str(i) or int(i)'
#NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
SCORE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

cohort2eid_fname = '/data/project/3dith/data/etc/cohort2eid.txt'
eid_dir = '/data/project/3dith/data/roadmap_epigenomics/'#{eid}_18_core_K27ac_mnemonics.bed'

cohort2eid = {}#tcga cohort to epigenome id (matched by tissue type)
f = open(cohort2eid_fname, 'r')
lines = f.readlines()
for current_line in lines:
#with open(cohort2eid_fname, 'r') as f:
    current_line = current_line.strip()
    current_key = current_line.split('\t')[0].strip()
    current_value = current_line.split('\t')[1].strip()
    cohort2eid[current_key] = current_value
f.close()

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    args.add_argument('-th', '--threshold', help = 'threshold for dmr; mean, mean_std, mean_2std', type = str, default = 'mean_std', required = True)
    args.add_argument('-b', '--binsize', help = 'binsize', type = int, default = int(1e6), required = False)
    args.add_argument('-df', '--dmr_feature', help = 'which features in DMR you will focus on. gene_reg or epigenome', type = str, default = 'gene_reg', required = True)
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    
    args.add_argument('--dmr_type', help = 'DMR type. TN (tumor-normal), HL (high score - low score) risk_HL', type = str, required = True)
    args.add_argument('--event', help = 'survival event', type = str, required = True)
    args.add_argument('--fold', help = 'fold number', type = str, required = True, default = 'fold1')
    args.add_argument('--version', help = 'feature version', type = str, required = True, default = 'v7.1')
    args.add_argument('--lr', help = 'learning rate. 0.001, 0.0005, 0.0001', type = float, required = True, default = 0.001)
    return args.parse_args()

def preprocess_gene_annot(category_, chrom):#chrom: chr{i}
    #change index of GENE_ANNOT dataframe to ENSG ID of each gene. 
    globals()['GENE_ANNOT_'+category_] = GENE_ANNOT[GENE_ANNOT.iloc[:,2] == category_].copy()
    if category_ == 'gene':
        GENE_ANNOT_indices = [x.split('gene_id')[1].split(';')[0].split('"')[1].split('.')[0] for x in globals()['GENE_ANNOT_'+category_].iloc[:,8].values]
    elif category_ == 'transcript':
        GENE_ANNOT_indices = [x.split('transcript_id')[1].split(';')[0].split('"')[1].split('.')[0] for x in globals()['GENE_ANNOT_'+category_].iloc[:,8].values]
    else:
        pass
    globals()['GENE_ANNOT_'+category_].index = GENE_ANNOT_indices
    
    assert GENE_ANNOT[GENE_ANNOT.iloc[:,2]==category_].shape[0] == globals()['GENE_ANNOT_'+category_].shape[0]
    
    current_chrom_annot = globals()['GENE_ANNOT_'+category_][globals()['GENE_ANNOT_'+category_].iloc[:,0]==chrom].copy()
    
    return current_chrom_annot

def preprocess_reg_annot(category_, chrom): #chrom: chr{i}
    #category in REG_ANNOT: 'open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', enhancer', 'promoter', 'promoter_flanking_region'
    #change index of REG_ANNOT to the id of its annotated item.  
    globals()['REG_ANNOT_'+category_] = REG_ANNOT[REG_ANNOT.iloc[:,2] == category_].copy()
    REG_ANNOT_indices = [x.split(':')[1].split(';')[0].strip() for x in globals()['REG_ANNOT_'+category_].iloc[:,8].values]
    globals()['REG_ANNOT_'+category_].index = REG_ANNOT_indices
    
    assert REG_ANNOT[REG_ANNOT.iloc[:,2]==category_].shape[0] == globals()['REG_ANNOT_'+category_].shape[0]
    
    chrom_ = str(chrom.split('chr')[1].strip())
    #print(chrom_)
    chrom_mask = np.array([str(x)==chrom_ for x in globals()['REG_ANNOT_'+category_].iloc[:,0].values])
    
    current_reg_annot = globals()['REG_ANNOT_'+category_].iloc[chrom_mask,:]
    
    return current_reg_annot


# In[18]:


def preprocess_epi_annot(category_, chrom, epi_annot):
    globals()['EPI_ANNOT_'+category_] = epi_annot[epi_annot.iloc[:,-1] == category_].copy()
    #chrom_mask = globals()['EPI_ANNOT_'+category_].iloc[:,0] == chrom
    chrom_mask = np.array([str(x)==chrom for x in globals()['EPI_ANNOT_'+category_].iloc[:,0].values])
    current_epi_annot = globals()['EPI_ANNOT_'+category_].iloc[chrom_mask,:]
    #display(current_epi_annot)
    assert current_epi_annot.iloc[:,0].unique() == chrom
    
    return current_epi_annot


if __name__ == '__main__':
    
    args = parse_arguments()
    
    os.chdir(args.working_dir)
    
    resultdir =os.path.join(os.getcwd(), 'result') #전체 result directory
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)

    #cohort_dir = os.path.join(resultdir, args.cohort) # 개별 cohort의 result directory
    #if not os.path.exists(cohort_dir):
    #    os.makedirs(cohort_dir)
    
    # 개별 cohort의 result를 저장하는 디렉토리
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
        
    print("cohort: {}".format(args.cohort))
    #print("cohort_dir: {}".format(cohort_dir))
    print(f'SAVEDIR: {SAVEDIR}')
    
    dmr_fname = os.path.join(SAVEDIR, 'DMR_binsize_'+str(args.binsize)+'.npz')
    print("DMR fname: {}".format(dmr_fname))
    dmr = np.load(dmr_fname, allow_pickle = True)
    
    if args.dmr_feature == 'gene_reg':
        # find DMR features (gene annotations, regulatory feature annotations)
        gene_annot_category = ['gene', 'transcript']
        reg_annot_category = REG_ANNOT.iloc[:,2].unique() 
        all_annot_category = np.concatenate((gene_annot_category, reg_annot_category), axis=None)

        for c in all_annot_category:
            globals()[c+'_ALL'] = {}
            #result_pickle_fname1 = os.path.join(cohort_dir, 'DMR_GENE-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
            #result_pickle_fname2 = os.path.join(cohort_dir, 'DMR_REG-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')

        for chrom in CHR_LIST:
            print("---\n"+chrom)
            processed_annot = {}
            for c in gene_annot_category:
                processed_annot['GENE_'+c] = preprocess_gene_annot(c, chrom)
            for c in reg_annot_category:
                processed_annot['REG_'+c] = preprocess_reg_annot(c, chrom)

            chrom_dmr = dmr[chrom+'_'+args.threshold+'_'+'bins'].copy()
            for i in range(len(chrom_dmr)):#iterate for each differentially-methylated genomic bin in this chromosome
                current_binname = chrom_dmr[i].copy()            
                start_ = int(current_binname.split('-')[0].split(':')[-1].strip())            
                end_ = int(current_binname.split('-')[-1].strip())            
                current_dmr_start = (start_ // int(args.binsize)) * int(args.binsize)            
                current_dmr_end = end_
                #find genes in this region
                for c in gene_annot_category:
                    #print("---\n"+c)
                    for k in range(processed_annot['GENE_'+c].shape[0]):
                        if int(processed_annot['GENE_'+c].iloc[k,3])>= int(current_dmr_start) and int(processed_annot['GENE_'+c].iloc[k,4]) <= int(current_dmr_end):
                            if current_binname not in list(globals()[c+'_ALL'].keys()):
                                globals()[c+'_ALL'][current_binname] = []
                            #globals()[c+'_ALL'][current_binname].append(processed_annot['GENE_'+c].index.values[k])
                            feature_to_append = str(processed_annot['GENE_'+c].iloc[k,0])+':'+str(processed_annot['GENE_'+c].iloc[k,3])+'-'                            +str(processed_annot['GENE_'+c].iloc[k,4]) + '_'+ processed_annot['GENE_'+c].index.values[k]
                            globals()[c+'_ALL'][current_binname].append(feature_to_append)

                #find regulatory features in this region
                for c in reg_annot_category:
                    #print("---\n"+c)
                    for k in range(processed_annot['REG_'+c].shape[0]):
                        if int(processed_annot['REG_'+c].iloc[k,3])>= int(current_dmr_start) and int(processed_annot['REG_'+c].iloc[k,4]) <= int(current_dmr_end):
                            if current_binname not in list(globals()[c+'_ALL'].keys()):
                                globals()[c+'_ALL'][current_binname] = []
                            #globals()[c+'_ALL'][current_binname].append(processed_annot['REG_'+c].index.values[k])
                            feature_to_append = 'chr'+str(processed_annot['REG_'+c].iloc[k,0])+':'+str(processed_annot['REG_'+c].iloc[k,3])                            +'-'+str(processed_annot['REG_'+c].iloc[k,4])+'_'+str(processed_annot['REG_'+c].index.values[k])
                            #print("feature to append (reg):", feature_to_append)#debug
                            globals()[c+'_ALL'][current_binname].append(feature_to_append)

                #find chromatin state in this region.
                #should specify epigenome id of each tissue type (i.e., cancer type.)

        for c in gene_annot_category:
            pickle_fname = os.path.join(SAVEDIR, 'DMR_GENE-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
            print(pickle_fname)
            with open(pickle_fname, 'wb') as f:
                pickle.dump(globals()[c+'_ALL'], f)
            f.close()

        for c in reg_annot_category:
            pickle_fname = os.path.join(SAVEDIR, 'DMR_REG-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
            print(pickle_fname)
            with open(pickle_fname, 'wb') as f:
                pickle.dump(globals()[c+'_ALL'], f)
            f.close()
        
    elif args.dmr_feature == 'epigenome':
        cohort = args.cohort
        if cohort in list(cohort2eid.keys()): # if this cohort has matching EID (epigenome id)
            cohort_eid = cohort2eid[cohort]
            epi_annot_fname = '/data/project/3dith/data/roadmap_epigenomics/'+cohort_eid+'_18_core_K27ac_mnemonics.bed'
            epi_annot_df = pd.read_csv(epi_annot_fname, sep = '\t', header = None)
            epi_annot_category = epi_annot_df.iloc[:,-1].unique()
            
            for c in epi_annot_category:
                globals()[c+'_ALL'] = {}#all DMR bin names which belong to this category. 
                #result_pickle_fname = os.path.join(cohort_dir, 'DMR_EPI-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')

            for chrom in CHR_LIST:
                print("---\n"+chrom)
                processed_annot = {}
                
                # In current chromosome, iterate for all epigenome feature categories.
                
                for c in epi_annot_category:
                    processed_annot['EPI_'+c] = preprocess_epi_annot(c, chrom, epi_annot_df)

                chrom_dmr = dmr[chrom+'_'+args.threshold+'_'+'bins'].copy()
                for i in range(len(chrom_dmr)):#iterate for each differentially-methylated genomic bin in this chromosome
                    current_binname = chrom_dmr[i].copy()            
                    start_ = int(current_binname.split('-')[0].split(':')[-1].strip())            
                    end_ = int(current_binname.split('-')[-1].strip())            
                    current_dmr_start = (start_ // int(args.binsize)) * int(args.binsize)            
                    current_dmr_end = end_

                    #find epigenome features in this region
                    for c in epi_annot_category:
                        #print("---\n"+c)
                        for k in range(processed_annot['EPI_'+c].shape[0]):
                            if int(processed_annot['EPI_'+c].iloc[k,1])>= int(current_dmr_start) and int(processed_annot['EPI_'+c].iloc[k,2]) <= int(current_dmr_end):
                                if current_binname not in list(globals()[c+'_ALL'].keys()):
                                    globals()[c+'_ALL'][current_binname] = []
                                feature_to_append = chrom+':'+str(processed_annot['EPI_'+c].iloc[k,1])+'-'+str(processed_annot['EPI_'+c].iloc[k,2]) #이 epigenome feature의 위치 정보. 
                                globals()[c+'_ALL'][current_binname].append(feature_to_append)


            for c in epi_annot_category:
                pickle_fname = os.path.join(SAVEDIR, 'EPI-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
                # remove existing EPI ANNOT pickle file 
                if os.path.exists(pickle_fname):
                    os.system("rm -rf %s" % pickle_fname)
                    print("Existing {} is removed. Writing {} from scratch again for clearance.".format(pickle_fname, pickle_fname))
                
                # write new pickle file (current result)
                if '/' in c:
                    c_remove_slash = c.replace('/', '-')
                    pickle_fname = os.path.join(SAVEDIR, 'EPI-ANNOT_'+c_remove_slash+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')

                print(pickle_fname)
                with open(pickle_fname, 'wb') as f:
                    pickle.dump(globals()[c+'_ALL'], f)
                f.close()
        else:
            print("This cohort has no matching Epigenome ID (EID)")
    else:
        pass

