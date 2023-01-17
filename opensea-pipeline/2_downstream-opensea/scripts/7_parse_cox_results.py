#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[6]:


import os
import sys
import glob
import argparse


# In[ ]:

# command: python3 7_parse_cox_results.py --version i --log_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/log --log_fname 7_parse_cox_results-v{i}.txt
def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--version', type = str, required = True)
    args.add_argument('--log_dir', type = str, required = True)
    args.add_argument('--log_fname', type = str, required = True)
    return args.parse_args()


# In[3]:


cohorts = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
cohorts.remove('TCGA-ESCA')
cohorts.remove('TCGA-HNSC')


# In[9]:


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
            #print("---")
            current_event = f.split('.csv')[0].split('-')[-1]
            df = pd.read_csv(f, index_col = 0)
            for i in range(df.shape[0]):
                current_p = float(df.p.values[i])
                if current_p <= p_threshold:
                    f_.write("---\n")
                    f_.write("{}\n".format(current_event))
                    f_.write('covariate: {}, coef: {}, p-value: {}\n'.format(df.index.values[i], df.coef.values[i], df.p.values[i]))
    f_.close()



