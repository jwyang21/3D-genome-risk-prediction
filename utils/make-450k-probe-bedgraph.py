
import pandas as pd
import numpy as np
import os
import sys
from tqdm.auto import tqdm
import argparse

cohorts = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'CpG type. island or shelf or shore of shelf_shore')
    return args.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    meta = pd.read_csv('/data/project/3dith/data/450k_metadata.'+args.cpg_type+'.sorted.bed', sep = '\t', names = ['chrom', 'start', 'end', 'name'])
    for cohort in cohorts:
        print("===\ncohort: {}".format(cohort))
        beta = pd.read_csv(f'/data/project/3dith/data/450k_xena/{cohort}.HumanMethylation450.tsv', sep = '\t', index_col = 0).dropna()
        beta_merged = beta.merge(meta, left_index = True, right_on = 'name')

        result_dir = os.path.join('/data/project/jeewon/research/3D-ITH/data/450k_bdgs_v2', cohort+'-'+args.cpg_type) #permission issue
        #result_dir = os.path.join('/data/project/3dith/data/450k_bdgs_v2/', cohort+'-'+args.cpg_type) #original result_dir
        print("result_dir: {}".format(result_dir))

        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        barcodes = [barcode for barcode in beta_merged.columns if barcode.startswith('TCGA-')]
        for barcode in tqdm(barcodes, desc=cohort):
            tmp = beta_merged[['chrom', 'start', 'end', barcode]].rename({barcode: 'score'}, axis=1).dropna()
            result_fname = os.path.join(result_dir, barcode+'.bedGraph')

            if barcodes.index(barcode)%100==0:
                print("result_fname: {}".format(result_fname)) #print result filename for every 100-th sample.

            tmp.to_csv(result_fname, sep = '\t', index = False, header = False)

            sorted_result_fname = os.path.join(result_dir, barcode+'.sorted.bedGraph')

            if barcodes.index(barcode)%100==0:
                print("sorted_result_fname: {}".format(sorted_result_fname)) #print result filename for every 100-th sample.

            command_ = 'bedtools sort -i {} > {} && rm {}'.format(result_fname, sorted_result_fname, result_fname)
            print(command_)

            os.system(command_)
