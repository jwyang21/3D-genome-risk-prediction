#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import sys
import argparse

# input: 각 TCGA cohort 및 그 cohort에서 사용해야 할 stem closeness filename (absolute path)가 포함되어 있는 csv file
# output: 입력받은 csv file에 있는 cohort들에 대해 cox regression을 돌리기 위한 python command들이 쓰여져 있는 file


survival_types = ['OS', 'DSS', 'DFI', 'PFI']
clinical_file = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'

# command: python3 5_write-cox-commands.py --cpg_type {cpg_type} --cox_version $i --command_fname 6_cox-commands --score_fname /data/project/3dith/data/cohort-1-best-score-km.csv
def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', type = str, required = True)
    args.add_argument('--command_fname', type = str, required = False, default = '6_cox-commands') #command 파일을 어떤 파일명으로 저장할 건지 #확장자 없이 입력. 
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
    if args.cox_version == '1': #covariate: age, gender, score, avg_beta
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]
            #score_file = scores.filename.values[i]
            score_file = scores.loc[cohort][f'filename_{args.cpg_type}']

            avg_beta_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/{args.cpg_type}_tumors_avg_beta.csv'
            result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/cox'

            for survival_type in survival_types:
                result_file = 'cox-regression-'+survival_type
                command_ = f'python3 cox_v{args.cox_version}.py --cohort {cohort} --survival_type {survival_type} ' \
                + f'--result_dir {result_dir}_v{args.cox_version} --result_file {result_file} --score_file {score_file} '   \
                + f'--avg_beta_file {avg_beta_file} --clinical_file {clinical_file}'
                #print(command_)
                #os.system(command_)
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))
    if args.cox_version == '2': #covariate: age, gender, score_group, avg_beta
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]
            #score_file = scores.filename.values[i]
            score_file = scores.loc[cohort][f'filename_{args.cpg_type}']
            
            score_group_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/sc-group.csv'
            avg_beta_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/{args.cpg_type}_tumors_avg_beta.csv'
            result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/cox'

            for survival_type in survival_types:
                result_file = 'cox-regression-'+survival_type
                command_ = f'python3 cox_v{args.cox_version}.py --cohort {cohort} --survival_type {survival_type} ' \
                + f'--result_dir {result_dir}_v{args.cox_version} --result_file {result_file} ' \
                + f'--use_score N --use_score_group Y --use_avg_beta Y --score_file {score_file} ' \
                + f'--score_group_file {score_group_file} --avg_beta_file {avg_beta_file} --clinical_file {clinical_file}'
                #print(command_)
                #os.system(command_)
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))
    
    if args.cox_version == '3': #covariate: age, gender, score
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]
            #score_file = scores.filename.values[i]
            score_file = scores.loc[cohort][f'filename_{args.cpg_type}']
            
            #score_group_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/sc-group.csv'
            avg_beta_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/{args.cpg_type}_tumors_avg_beta.csv'
            result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/cox'
            

            for survival_type in survival_types:
                result_file = 'cox-regression-'+survival_type
                command_ = f'python3 cox_v{args.cox_version}.py --cohort {cohort} --survival_type {survival_type} ' \
                + f'--result_dir {result_dir}_v{args.cox_version} --result_file {result_file} ' \
                + f'--use_score Y --use_score_group N --use_avg_beta N --score_file {score_file} ' \
                + f'--avg_beta_file {avg_beta_file} --clinical_file {clinical_file}'
                #print(command_)
                #os.system(command_)
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))
        
    if args.cox_version == '4': #covariate: age, gender, score_group
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]
            #score_file = scores.filename.values[i]
            #score_file = scores.loc[cohort][f'filename_{args.cpg_type}']
            
            score_group_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/sc-group.csv'
            avg_beta_file = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/{args.cpg_type}_tumors_avg_beta.csv'
            result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result/{cohort}/cox'

            for survival_type in survival_types:
                result_file = 'cox-regression-'+survival_type
                command_ = f'python3 cox_v{args.cox_version}.py --cohort {cohort} --survival_type {survival_type} ' \
                + f'--result_dir {result_dir}_v{args.cox_version} --result_file {result_file} ' \
                + f'--use_score N --use_score_group Y --use_avg_beta N ' \
                + f'--score_group_file {score_group_file} --avg_beta_file {avg_beta_file} --clinical_file {clinical_file}'
                #print(command_)
                #os.system(command_)
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))
        
    if args.cox_version == '5': #covariate: age, gender, score_group, beta_group
        for i in range(scores.shape[0]):
        #for i in range(1):
            cohort = scores.index.values[i]
            #score_file = scores.filename.values[i]
            #score_file = scores.loc[cohort][f'filename_{args.cpg_type}']
            
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
                #print(command_)
                #os.system(command_)
                f.write(command_)
                f.write('\n\n')
        f.close()
        print("command file: {}".format(command_full_name))

