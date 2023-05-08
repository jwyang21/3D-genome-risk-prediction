#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import csv

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('--dmr_type', help = 'dmr type', type = str, default = 'HL', required = True)    
    parser.add_argument('--cpg_type', help = 'island, opensea, shelf_shore', type = str, required = True)    
    parser.add_argument('--cohort2eid', help = 'cohort2eid fname', type = str, default = '/data/project/3dith/data/etc/cohort2eid.txt', required = False)
    parser.add_argument('--chromatin_states', help = 'chromatin state fname', type = str, default = '/data/project/3dith/data/chromatin_states.npy', required = False)
    parser.add_argument('--event', help = 'survival event', type = str, required = True)
    parser.add_argument('--fold', help = 'fold number', type = str, required = True, default = 'fold1')
    parser.add_argument('--version', help = 'feature version', type = str, required = True, default = 'v7.1')
    parser.add_argument('--lr', help = 'learning rate. 0.001, 0.0005, 0.0001', type = float, required = True, default = 0.001)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    parser.add_argument('--input_fname', help = 'input filename', type = str, default = 'DMR_EPI_features_threshold_mean_std_len.npz', required = False)
    return parser.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    working_dir = args.working_dir    
    save_dir = os.path.join(working_dir, 'result') 
    cohort2eid = pd.read_csv(args.cohort2eid, sep = '\t', header = None, names = ['cohort', 'eid'])
    chromatin_states = np.load(args.chromatin_states, allow_pickle = True)
    
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
    SAVEDIR_basename = os.path.basename(SAVEDIR)
    
    
    save_fname_l = f'EPI-category-len-stacked-bar-chart-{args.cpg_type}-{args.dmr_type}.csv' 
    save_full_fname_l = os.path.join(save_dir, save_fname_l) 
    
    if os.path.exists(save_full_fname_l):
        f_ = open(save_full_fname_l, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        f_ = open(save_full_fname_l, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        colnames = ['cohort', 'category']
        for i in range(len(chromatin_states)):
            colnames.append(chromatin_states[i])
        f_writer.writerow(colnames)
        
    if args.cohort in cohort2eid.cohort.values:
        print(args.cohort)
        if args.dmr_type == 'HL':
            input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/4_HL-DMR-{args.cpg_type}/result/{cohort}'
        elif args.dmr_type == 'TN':
            input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/result/{cohort}'
        elif args.dmr_type == 'risk_HL':
            input_dir = os.path.join(args.working_dir, 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
        else:
            pass
        
        data = np.load(os.path.join(input_dir, args.input_fname))
        
        cnt = {} 
                
        for k in list(data.keys()):
            for s in chromatin_states:
                if s not in list(cnt.keys()): 
                    cnt[s] = 0
                if s in k:
                    cnt[s] += data[k]
        
        content = []
        content.append(args.cohort)
        content.append(SAVEDIR_basename)
        for k in list(cnt.keys()):
            content.append(cnt[k])
        f_writer.writerow(content)
                    
    print(os.path.join(save_dir, save_fname_l))
    f_.close()
    
    all_df = pd.read_csv(os.path.join(save_dir, save_fname_l))
    save_fname_p = f'EPI-category-proportion-stacked-bar-chart-{args.cpg_type}-{args.dmr_type}.csv'
    save_full_fname_p = os.path.join(save_dir, save_fname_p) 
    
    if os.path.exists(save_full_fname_p):
        f_ = open(save_full_fname_p, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        f_ = open(save_full_fname_p, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        colnames = []
        colnames.append('cohort')
        colnames.append('category')
        for i in range(len(chromatin_states)):
            colnames.append(chromatin_states[i])
        f_writer.writerow(colnames)
    
    if args.cohort in cohort2eid.cohort.values:
        content = []
        content.append(args.cohort)
        content.append(SAVEDIR_basename)
        print('all_df')
        print(all_df)
        all_df_data = all_df[all_df['cohort']==args.cohort].copy()
        all_df_data = all_df_data[all_df_data['category']==SAVEDIR_basename].copy()
        
        for s in chromatin_states:
            current_state_proportion = all_df_data[s].values[0]/np.sum(all_df_data[chromatin_states].values)
            content.append(current_state_proportion)
            print(s, current_state_proportion)
        f_writer.writerow(content)
    print(os.path.join(save_dir, save_fname_p))    
    f_.close()
