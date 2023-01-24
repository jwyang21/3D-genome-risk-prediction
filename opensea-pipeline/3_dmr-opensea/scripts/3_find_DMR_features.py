import pandas as pd
import numpy as np
import os
import pickle
import argparse

CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
GENE_ANNOT = pd.read_csv('/data/project/3dith/data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf', sep = '\t', header = None, skiprows = [0, 1, 2, 3, 4])
REG_ANNOT = pd.read_csv('/data/project/3dith/data/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff', sep = '\t', header = None)
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
cohort2eid_fname = '/data/project/3dith/data/etc/cohort2eid.txt'
eid_dir = '/data/project/3dith/data/roadmap_epigenomics/'
cohort2eid = {}
f = open(cohort2eid_fname, 'r')
lines = f.readlines()
for current_line in lines:
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
    return args.parse_args()

def preprocess_gene_annot(category_, chrom):
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

def preprocess_reg_annot(category_, chrom):
    globals()['REG_ANNOT_'+category_] = REG_ANNOT[REG_ANNOT.iloc[:,2] == category_].copy()
    REG_ANNOT_indices = [x.split(':')[1].split(';')[0].strip() for x in globals()['REG_ANNOT_'+category_].iloc[:,8].values]
    globals()['REG_ANNOT_'+category_].index = REG_ANNOT_indices
    
    assert REG_ANNOT[REG_ANNOT.iloc[:,2]==category_].shape[0] == globals()['REG_ANNOT_'+category_].shape[0]
    
    chrom_ = str(chrom.split('chr')[1].strip())
    chrom_mask = np.array([str(x)==chrom_ for x in globals()['REG_ANNOT_'+category_].iloc[:,0].values])
    
    current_reg_annot = globals()['REG_ANNOT_'+category_].iloc[chrom_mask,:]
    
    return current_reg_annot

def preprocess_epi_annot(category_, chrom, epi_annot):
    globals()['EPI_ANNOT_'+category_] = epi_annot[epi_annot.iloc[:,-1] == category_].copy()
    chrom_mask = np.array([str(x)==chrom for x in globals()['EPI_ANNOT_'+category_].iloc[:,0].values])
    current_epi_annot = globals()['EPI_ANNOT_'+category_].iloc[chrom_mask,:]
    assert current_epi_annot.iloc[:,0].unique() == chrom
    
    return current_epi_annot

if __name__ == '__main__':
    
    args = parse_arguments()
    
    os.chdir(args.working_dir)
    
    resultdir =os.path.join(os.getcwd(), 'result') 
    if not os.path.exists(resultdir):
        os.makedirs(resultdir)

    cohort_dir = os.path.join(resultdir, args.cohort) 
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
        
    print("cohort: {}".format(args.cohort))
    print("cohort_dir: {}".format(cohort_dir))
    
    dmr_fname = os.path.join(cohort_dir, 'DMR_binsize_'+str(args.binsize)+'.npz')
    print("DMR fname: {}".format(dmr_fname))
    dmr = np.load(dmr_fname, allow_pickle = True)
    
    if args.dmr_feature == 'gene_reg':
        gene_annot_category = ['gene', 'transcript']
        reg_annot_category = REG_ANNOT.iloc[:,2].unique() 
        all_annot_category = np.concatenate((gene_annot_category, reg_annot_category), axis=None)

        for c in all_annot_category:
            globals()[c+'_ALL'] = {}

        for chrom in CHR_LIST:
            print("---\n"+chrom)
            processed_annot = {}
            for c in gene_annot_category:
                processed_annot['GENE_'+c] = preprocess_gene_annot(c, chrom)
            for c in reg_annot_category:
                processed_annot['REG_'+c] = preprocess_reg_annot(c, chrom)

            chrom_dmr = dmr[chrom+'_'+args.threshold+'_'+'bins'].copy()
            for i in range(len(chrom_dmr)):
                current_binname = chrom_dmr[i].copy()            
                start_ = int(current_binname.split('-')[0].split(':')[-1].strip())            
                end_ = int(current_binname.split('-')[-1].strip())            
                current_dmr_start = (start_ // int(args.binsize)) * int(args.binsize)            
                current_dmr_end = end_
                for c in gene_annot_category:
                    for k in range(processed_annot['GENE_'+c].shape[0]):
                        if int(processed_annot['GENE_'+c].iloc[k,3])>= int(current_dmr_start) and int(processed_annot['GENE_'+c].iloc[k,4]) <= int(current_dmr_end):
                            if current_binname not in list(globals()[c+'_ALL'].keys()):
                                globals()[c+'_ALL'][current_binname] = []
                            feature_to_append = str(processed_annot['GENE_'+c].iloc[k,0])+':'+str(processed_annot['GENE_'+c].iloc[k,3])+'-'  +str(processed_annot['GENE_'+c].iloc[k,4]) + '_'+ processed_annot['GENE_'+c].index.values[k]
                            globals()[c+'_ALL'][current_binname].append(feature_to_append)

                for c in reg_annot_category:
                    for k in range(processed_annot['REG_'+c].shape[0]):
                        if int(processed_annot['REG_'+c].iloc[k,3])>= int(current_dmr_start) and int(processed_annot['REG_'+c].iloc[k,4]) <= int(current_dmr_end):
                            if current_binname not in list(globals()[c+'_ALL'].keys()):
                                globals()[c+'_ALL'][current_binname] = []
                            feature_to_append = 'chr'+str(processed_annot['REG_'+c].iloc[k,0])+':'+str(processed_annot['REG_'+c].iloc[k,3]) +'-'+str(processed_annot['REG_'+c].iloc[k,4])+'_'+str(processed_annot['REG_'+c].index.values[k])
                            globals()[c+'_ALL'][current_binname].append(feature_to_append)

        for c in gene_annot_category:
            pickle_fname = os.path.join(cohort_dir, 'DMR_GENE-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
            print(pickle_fname)
            with open(pickle_fname, 'wb') as f:
                pickle.dump(globals()[c+'_ALL'], f)
            f.close()

        for c in reg_annot_category:
            pickle_fname = os.path.join(cohort_dir, 'DMR_REG-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
            print(pickle_fname)
            with open(pickle_fname, 'wb') as f:
                pickle.dump(globals()[c+'_ALL'], f)
            f.close()
        
    elif args.dmr_feature == 'epigenome':
        cohort = args.cohort
        if cohort in list(cohort2eid.keys()): 
            cohort_eid = cohort2eid[cohort]
            epi_annot_fname = '/data/project/3dith/data/roadmap_epigenomics/'+cohort_eid+'_18_core_K27ac_mnemonics.bed'
            epi_annot_df = pd.read_csv(epi_annot_fname, sep = '\t', header = None)
            epi_annot_category = epi_annot_df.iloc[:,-1].unique()
            
            for c in epi_annot_category:
                globals()[c+'_ALL'] = {}

            for chrom in CHR_LIST:
                print("---\n"+chrom)
                processed_annot = {}
                
                for c in epi_annot_category:
                    processed_annot['EPI_'+c] = preprocess_epi_annot(c, chrom, epi_annot_df)

                chrom_dmr = dmr[chrom+'_'+args.threshold+'_'+'bins'].copy()
                for i in range(len(chrom_dmr)):
                    current_binname = chrom_dmr[i].copy()            
                    start_ = int(current_binname.split('-')[0].split(':')[-1].strip())            
                    end_ = int(current_binname.split('-')[-1].strip())            
                    current_dmr_start = (start_ // int(args.binsize)) * int(args.binsize)            
                    current_dmr_end = end_

                    for c in epi_annot_category:
                        for k in range(processed_annot['EPI_'+c].shape[0]):
                            if int(processed_annot['EPI_'+c].iloc[k,1])>= int(current_dmr_start) and int(processed_annot['EPI_'+c].iloc[k,2]) <= int(current_dmr_end):
                                if current_binname not in list(globals()[c+'_ALL'].keys()):
                                    globals()[c+'_ALL'][current_binname] = []
                                feature_to_append = chrom+':'+str(processed_annot['EPI_'+c].iloc[k,1])+'-'+str(processed_annot['EPI_'+c].iloc[k,2]) 
                                globals()[c+'_ALL'][current_binname].append(feature_to_append)


            for c in epi_annot_category:
                pickle_fname = os.path.join(cohort_dir, 'EPI-ANNOT_'+c+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')
                if os.path.exists(pickle_fname):
                    os.system("rm -rf %s" % pickle_fname)
                    print("Existing {} is removed. Writing {} from scratch again for clearance.".format(pickle_fname, pickle_fname))
                
                if '/' in c:
                    c_remove_slash = c.replace('/', '-')
                    pickle_fname = os.path.join(cohort_dir, 'EPI-ANNOT_'+c_remove_slash+'_threshold_'+str(args.threshold)+'_binsize_'+str(args.binsize)+'.pickle')

                print(pickle_fname)
                with open(pickle_fname, 'wb') as f:
                    pickle.dump(globals()[c+'_ALL'], f)
                f.close()
        else:
            print("This cohort has no matching Epigenome ID (EID)")
    else:
        pass
