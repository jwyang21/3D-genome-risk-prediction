import pandas as pd
import numpy as np
import os
import argparse

COHORTS = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ') 

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'CpG type. island, opensea, shelf_shore', type = str, required = True)
    return args.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    log = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'log', 'pcc-avg_beta-stem_closeness-ALL.log')
    working_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type)
    result_dir = os.path.join('/data/project/3dith/pipelines/', args.cpg_type+'-pipeline', '2_downstream-'+args.cpg_type, 'result')
    
    all_df = pd.DataFrame(np.zeros((len(COHORTS), 2), dtype = float), index = COHORTS, columns = ['pcc_all', 'pcc_T'])
    
    all_lines = []
    
    f = open(log, 'r') 
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        all_lines.append(line)
    f.close()
    
    for i in range(len(all_lines)):
        line = all_lines[i]
        if line.startswith(f'PCC(avg_beta({args.cpg_type}), stem-closeness) in'):
            cohort = line.split('in ')[1].split(',')[0]
            if line.endswith('(Tumor + Normal)'):
                pcc = float(all_lines[i+1].split('pcc: ')[1].split(',')[0].replace(' ',''))
                all_df.loc[cohort]['pcc_all'] = pcc
            elif line.endswith('tumor samples'):
                pcc = float(all_lines[i+1].split('pcc: ')[1].split(',')[0].replace(' ',''))
                all_df.loc[cohort]['pcc_T'] = pcc
                
    all_df.to_csv(os.path.join(result_dir, 'pcc_sc_avg-beta.csv'))
    print("result file: ", os.path.join(result_dir, 'pcc_sc_avg-beta.csv'))
