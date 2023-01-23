import pandas as pd
import numpy as np
import os
import sys
import argparse

survival_types = ['OS', 'DSS', 'DFI', 'PFI']
clinical_file = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'

# command: python3 5_write-cox-commands.py --cpg_type {cpg_type} --cox_version $i --command_fname 6_cox-commands --score_fname /data/project/3dith/data/cohort-1-best-score-km.csv
def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', type = str, required = True)
    args.add_argument('--command_fname', type = str, required = False, default = '6_cox-commands') 
    args.add_argument('--cox_version', type = str, required = True)
    args.add_argument('--score_fname', help = 'stem-closeness score file name, based on opensea K-M result', type = str, default = '/data/project/3dith/data/cohort-1-best-score-km.csv', required = True)
    return args.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    #scores_fname = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/cohort-1-best-score-km.csv'
    scores = pd.read_csv(args.score_fname, index_col = 0)
        
    command_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/scripts'
    
    command_full_name = f'{command_dir}/{args.command_fname}_v{args.cox_version}.sh'
    if not os.path.exists(command_dir):
        os.makedirs(command_dir)
    
    f = open(command_full_name, 'w')
        
    if args.cox_version == '5':
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]            
            score_group_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/sc-group.csv'
            avg_beta_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/{args.cpg_type}_tumors_avg_beta.csv'
            avg_beta_group_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/avg-beta-group.csv'
            result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/cox'

            for survival_type in survival_types:
                result_file = 'cox-regression-'+survival_type
                command_ = f'python3 cox_v{args.cox_version}.py --cohort {cohort} --survival_type {survival_type} ' \
                + f'--result_dir {result_dir}_v{args.cox_version} --result_file {result_file} ' \
                + f'--use_score N --use_score_group Y --use_avg_beta N --use_beta_group Y ' \
                + f'--score_group_file {score_group_file} --avg_beta_file {avg_beta_file} --beta_group_file {avg_beta_group_file} --clinical_file {clinical_file}'
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))
