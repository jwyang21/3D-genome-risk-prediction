#!/usr/bin/env python
# coding: utf-8

# In[27]:

# 전체 DMR에서 각 chromatin state가 차지하는 길이 및 비율 구하기 -> 나중에 stacked bar chart 그릴 때 데이터로 활용.

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
    
    #+
    parser.add_argument('--cpg_type', help = 'island, opensea, shelf_shore', type = str, required = True)    
    #parser.add_argument('-s_dir', '--save_dir', help = 'saving directory', type = str, required = True)
    
    #0
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
    #working_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/'
    #save_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/result/'
    #working_dir, save_dir = args.working_dir, args.save_dir
    working_dir = args.working_dir
    
    save_dir = os.path.join(working_dir, 'result') #전체 cohort들의 결과를 한꺼번에 저장.
    cohort2eid = pd.read_csv(args.cohort2eid, sep = '\t', header = None, names = ['cohort', 'eid'])
    chromatin_states = np.load(args.chromatin_states, allow_pickle = True)
    
    
    # SAVEDIR: 현재 cohort의 결과를 저장하는 디렉토리.
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
    SAVEDIR_basename = os.path.basename(SAVEDIR)
    
    
    save_fname_l = f'EPI-category-len-stacked-bar-chart-{args.cpg_type}-{args.dmr_type}.csv' #l: length (of each chromatin state)
    
    #all_df = pd.DataFrame(np.zeros((len(cohort2eid.cohort.values), len(chromatin_states)), dtype = float), index = cohort2eid.cohort.values, columns = chromatin_states)
    save_full_fname_l = os.path.join(save_dir, save_fname_l) #length
    
    if os.path.exists(save_full_fname_l):
        #f_ = open(save_full_fname_l, 'a')
        f_ = open(save_full_fname_l, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        #f_ = open(save_full_fname_l, 'w')
        f_ = open(save_full_fname_l, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        #f_writer.writerow(['cohort', 'category', 'mean_std'])
        colnames = ['cohort', 'category']
        for i in range(len(chromatin_states)):
            colnames.append(chromatin_states[i])
        f_writer.writerow(colnames)
    
    '''
    if cohort in cohort2eid.cohort.values:
        if os.path.exists(save_full_fname_l):
            # write column names
            f_.write('cohort')
            f_.write(',')
            f_.write('category')
            f_.write(SAVEDIR_basename)
            f_.write(',')
            for i in range(len(chromatin_states)-1):
                f_.write(chromatin_states[i])
                f_.write(',')
            f_.write(chromatin_states[-1])
            f_.write('\n')    
    '''
    
    #for cohort in cohort2eid.cohort.values:
    if args.cohort in cohort2eid.cohort.values:
        print(args.cohort)
        if args.dmr_type == 'HL':
            input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/4_HL-DMR-{args.cpg_type}/result/{cohort}'
        elif args.dmr_type == 'TN':
            input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/result/{cohort}'
        elif args.dmr_type == 'risk_HL':
            #input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/5_risk-HL-DMR-{args.cpg_type}/result/{cohort}'
            input_dir = os.path.join(args.working_dir, 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
        else:
            pass
        
        data = np.load(os.path.join(input_dir, args.input_fname))
        
        cnt = {} #각 chromatin state의 총 길이를 셈 (sum over 22 autosomes)
                
        for k in list(data.keys()):
            for s in chromatin_states:
                if s not in list(cnt.keys()): 
                    cnt[s] = 0
                if s in k:
                    #print(data[k])
                    #all_df.loc[cohort][s] += data[k]
                    cnt[s] += data[k]
        
        # write data
        content = []
        content.append(args.cohort)
        content.append(SAVEDIR_basename)
        for k in list(cnt.keys()):
            content.append(cnt[k])
        f_writer.writerow(content)
        '''
        f_.write(cohort)
        f_.write(',')
        for k in list(cnt.keys()):
            f_.write(cnt[s])
            if k != list(cnt.keys())[1]:
                f_.write(',')
            else:
                f_.write('\n')    
        '''
                    
    #all_df.to_csv(os.path.join(save_dir, save_fname_l))
    print(os.path.join(save_dir, save_fname_l))
    f_.close()
    
    #-------------------- all_df (proportion)
    #all_df_proportion = pd.DataFrame(np.zeros((len(cohort2eid.cohort.values), len(chromatin_states)), dtype = float), index = cohort2eid.cohort.values, columns = chromatin_states)
    #all_df = pd.read_csv(os.path.join(save_dir, save_fname_l), index_col = 0) #length data
    all_df = pd.read_csv(os.path.join(save_dir, save_fname_l))
    save_fname_p = f'EPI-category-proportion-stacked-bar-chart-{args.cpg_type}-{args.dmr_type}.csv'
    save_full_fname_p = os.path.join(save_dir, save_fname_p) #p: proportion (of each chromatin state)
    
    if os.path.exists(save_full_fname_p):
        f_ = open(save_full_fname_p, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        f_ = open(save_full_fname_p, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        # 새로 파일을 만드는 경우엔 column name을 맨 처음에 써줘야 함.
        colnames = []
        colnames.append('cohort')
        colnames.append('category')
        for i in range(len(chromatin_states)):
            colnames.append(chromatin_states[i])
        f_writer.writerow(colnames)
    
    #for cohort in cohort2eid.cohort.values:
    if args.cohort in cohort2eid.cohort.values:
        '''
        if not os.path.exists(save_full_fnape_p):
            # 새로 파일을 만드는 경우엔 column name을 맨 처음에 써줘야 함. 
            # write column names
            f_.write('cohort')
            f_.write(',')
            f_.write('category')
            f_.write(SAVEDIR_basename)
            f_.write(',')
            for i in range(len(chromatin_states)-1):
                f_.write(chromatin_states[i])
                f_.write(',')
            f_.write(chromatin_states[-1])
            f_.write('\n')  
        '''
        #cnt = {}
        content = []
        content.append(args.cohort)
        content.append(SAVEDIR_basename)
        print('all_df')
        print(all_df)
        all_df_data = all_df[all_df['cohort']==args.cohort].copy()
        all_df_data = all_df_data[all_df_data['category']==SAVEDIR_basename].copy()
        
        for s in chromatin_states:
            #if s not in list(cnt.keys()):
                #cnt[s] = 0
            #all_df_proportion.loc[cohort][s] = all_df.loc[cohort][s]/all_df.loc[cohort].sum()
            #all_df_proportion.loc[args.cohort][s] = all_df.loc[args.cohort][s]/all_df.loc[args.cohort].sum()
            #current_state_proportion = all_df.loc[args.cohort][s]/all_df.loc[args.cohort].sum()
            current_state_proportion = all_df_data[s].values[0]/np.sum(all_df_data[chromatin_states].values)
            content.append(current_state_proportion)
            print(s, current_state_proportion)
        f_writer.writerow(content)
        
        
        
        #for k in list(cnt.keys()):
        #    content.append(cnt[s])
        #f_writer.writerow(content)
        '''
        f_.write(cohort)
        f_.write(',')
        for k in list(cnt.keys()):
            f_.write(cnt[s])
            if k != list(cnt.keys())[1]:
                f_.write(',')
            else:
                f_.write('\n') 
        ''' 
    #all_df_proportion.to_csv(os.path.join(save_dir, save_fname_p))
    print(os.path.join(save_dir, save_fname_p))    
    f_.close()