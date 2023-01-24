import pandas as pd
import numpy as np
import os
import sys
import glob
import argparse

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--version', type = str, required = True)
    args.add_argument('--log_dir', type = str, required = True)
    args.add_argument('--log_fname', type = str, required = True)
    return args.parse_args()

cohorts = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
cohorts.remove('TCGA-ESCA')
cohorts.remove('TCGA-HNSC')

p_threshold = 5e-2

if __name__=='__main__':
    args = parse_arguments()
    
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)
    print("result filename: {}".format(os.path.join(args.log_dir, args.log_fname)))
    f_ = open(os.path.join(args.log_dir, args.log_fname), 'w')
    
    for cohort in cohorts:
        f_.write("\n===\ncohort: {}\n".format(cohort))
        files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/{cohort}/cox_v{args.version}/*.csv')

        for f in files:
            current_event = f.split('.csv')[0].split('-')[-1]
            df = pd.read_csv(f, index_col = 0)
            for i in range(df.shape[0]):
                current_p = float(df.p.values[i])
                if current_p <= p_threshold:
                    f_.write("---\n")
                    f_.write("{}\n".format(current_event))
                    f_.write('covariate: {}, coef: {}, p-value: {}\n'.format(df.index.values[i], df.coef.values[i], df.p.values[i]))
    f_.close()
