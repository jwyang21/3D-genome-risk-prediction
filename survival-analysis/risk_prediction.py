#!/usr/bin/env python
# coding: utf-8

# Modified the codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/notebooks/DeepSurv%20Example.ipynb

import sys
sys.path.append('/data/project/jeewon/research/DeepSurv/deepsurv')
import deep_surv
from deepsurv_logger import DeepSurvLogger, TensorboardLogger
import utils
import viz
import numpy as np
import pandas as pd
import os
import lasagne
import matplotlib
import matplotlib.pyplot as plt
import sys
import lifelines
import argparse
import pickle
import time
import math

# set default params for matplotlib.pyplot
matplotlib.rcParams['figure.dpi'] = 300 #150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

P_THRESHOLD = 5e-2
all_survival_targets = ['OS', 'DSS', 'DFI', 'PFI']

def parse_arguments():
    # necessary
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_type', type = str, help = 'opensea, island, shelf_shore', required = True)
    parser.add_argument('--cohort', type = str, help = 'TCGA cohort', required = True)
    parser.add_argument('--version', type = str, help = 'code version', required = True, default = None)
    parser.add_argument('--result_dir', type = str, help = 'result directory', required = True, default = '/data/project/3dith/result')
    parser.add_argument('--data_dir', type = str, help = 'data directory. train,valid,test sets are stored here', required = True, default = '/data/project/3dith/data/')
    parser.add_argument('--num_folds', type = int, help = 'number of K-fold CV', required = True, default = 5)
    parser.add_argument('--lr', type = float, help = 'initial learning rate', required = True, default = 1e-05)
    
    # optional
    parser.add_argument('--test_df_with_risk_load', type = str, help = 'whether to load test_df_with_risk file. Y/N', required = False, default = 'N')
    parser.add_argument('--model_load', type = str, help = 'whether to load model. Y/N', required = False, default = 'N')
    parser.add_argument('--load_model_name', type = str, help = 'fname of model to be loaded', required =False, default = '')
    parser.add_argument('--load_model_weight_name', type = str, help = 'fname of model weight to be loaded', required = False, default = '')
    parser.add_argument('--load_metrics_name', type = str, help = 'fname of metrics, train log, to be loaded', required = False, default = '')
    
    return parser.parse_args()

# used codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/notebooks/DeepSurv%20Example.ipynb#L77
def dataframe_to_deepsurv_ds(df, event_col = 'Event', time_col = 'Time'):
    # input df should include columns 'Event' and 'Time' 
    # Extract the event and time columns as numpy arrays
    e = df[event_col].values.astype(np.int32) #event indicator
    t = df[time_col].values.astype(np.float32) # failure time. 

    x_df = df.drop([event_col, time_col], axis = 1)
    x = x_df.values.astype(np.float32)
    
    return {
        'x' : x,
        'e' : e,
        't' : t
    }

def survival_analysis_single_target(df, t, directory, fig_width, figname, cohort):#t = target
    fig = plt.figure(figsize = (fig_width, fig_width)) 

    valid_t = []
    valid_t_pvalue = []
    sig_t = []
    sig_t_pvalue = []

    d = df[df[f'Time'].notnull() & df[f'Event'].notnull()].copy()
    if d.shape[0] !=0:
        valid_t.append(t)

        print("target = {}, num_samples = {}".format(t, d.shape[0]))

        groups = df.risk_group.unique()

        ax = fig.add_subplot(111)

        for group in groups:#group 0 (low), group 1(High)
            T = d[d.risk_group==group][f'Time'].values
            E = d[d.risk_group==group][f'Event'].astype(bool).values #0, 1을 False, True로 바꿈.

            kmf = lifelines.KaplanMeierFitter()
            kmf.fit(T, E, label = group)

            kmf.plot_survival_function(ax = ax, ci_show = False, linewidth = 3, xticks = [0, 2000], yticks = [0.2, 0.4, 0.6, 0.8, 1.0])
        ax.get_xaxis().set_visible(True)
        ax.grid(False)
        ax.set_facecolor('white')
        ax.spines['top'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['right'].set_color('black')

        res = lifelines.statistics.logrank_test( #input: T_group0, T_group1, E_group0, E_group1
            d[d.risk_group=='Low'][f'Time'].values,#group0
            d[d.risk_group=='High'][f'Time'].values,#group1
            d[d.risk_group=='Low'][f'Event'].values,#group0
            d[d.risk_group=='High'][f'Event'].values#group1
        )
        ax.legend(frameon = False)
        valid_t_pvalue.append(res.p_value)
        if res.p_value < P_THRESHOLD:
            sig_t.append(t)
            sig_t_pvalue.append(res.p_value)
        cohort_type = cohort.split('-')[1]
        ax.set_title(f'{cohort_type} ({t}, p = {res.p_value:.2g})', pad = 5)
    fig.subplots_adjust(wspace = 0.3)
    fig.tight_layout()
    print("figure file: {}".format(os.path.join(directory, figname)))
    plt.savefig(os.path.join(directory, figname))

    return valid_t, valid_t_pvalue, sig_t, sig_t_pvalue

# used codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/deepsurv/viz.py#L18
def extract_value_list(arr):
    return list(np.array(arr)[:,1])
  
# used codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/deepsurv/viz.py#L21
def plot_log_nll(log, savedir, fig_name):
    TRAIN_LOSS = 'loss'
    TRAIN_CI = 'c-index'
    VALID_LOSS = 'valid_loss'
    VALID_CI = 'valid_c-index'
    num_epochs = len(log[TRAIN_LOSS])
    
    fig = plt.figure(figsize = (3, 3))
    ax1 = fig.add_subplot(111)
    handles = []
    if TRAIN_LOSS in log:
        epochs = range(num_epochs)
        values = extract_value_list(log[TRAIN_LOSS])
        train, = ax1.plot(epochs, values, 'b', label = 'Training')
        handles.append(train)
    if VALID_LOSS in log:
        epochs = np.linspace(0,num_epochs-1,num=len(log[VALID_LOSS]))
        values = extract_value_list(log[VALID_LOSS])
        valid, = ax1.plot(epochs,values, 'r', label = 'Validation')
        #ax2.tick_params('y', colors='r')
        handles.append(valid)
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Negative Log Likelihood')
    plt.legend(handles = handles)
    fig.tight_layout()
    plt.savefig(os.path.join(savedir, fig_name))
    print(os.path.join(savedir, fig_name))

# used codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/deepsurv/viz.py#L21
def plot_log_c_index(log, savedir, fig_name):
    TRAIN_LOSS = 'loss'
    TRAIN_CI = 'c-index'
    VALID_LOSS = 'valid_loss'
    VALID_CI = 'valid_c-index'
    num_epochs = len(log[TRAIN_LOSS])
    
    fig = plt.figure(figsize = (3, 3))
    ax1 = fig.add_subplot(111)

    handles = []
    if TRAIN_CI in log:
        epochs = np.linspace(0,num_epochs-1,num=len(log[TRAIN_CI]))
        train, = ax1.plot(epochs, extract_value_list(log[TRAIN_CI]), label = 'Training')
        handles.append(train)
    if VALID_CI in log:
        epochs = np.linspace(0,num_epochs-1,num=len(log[VALID_CI]))
        valid, = ax1.plot(epochs, extract_value_list(log[VALID_CI]), label = 'Validation')
        handles.append(valid)
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Concordance Index')
    plt.legend(handles = handles)    
    fig.tight_layout()
    plt.savefig(os.path.join(savedir, fig_name))
    print(os.path.join(savedir, fig_name))

if __name__ == '__main__':
    # used the codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/notebooks/DeepSurv%20Example.ipynb#L77, after proper modification
    args = parse_arguments()
    cohort_type = args.cohort.split('-')[1]
    
    colnames = ['target', 'fold', 'train_c', 'valid_c', 'test_c', 'LogRank_p', 'time', 'model', 'weight', 'metrics', 'train_valid_nll', 'train_valid_c_plot', 'LogRank_plot']
    rownames = []
    for i in range(args.num_folds):
        for target in all_survival_targets:
            rownames.append(f'fold{i+1}_{target}')

    result_df = pd.DataFrame(columns = colnames, index = rownames)
    
    cohort_result_dir = os.path.join(args.result_dir, args.cpg_type, args.cohort, f'v{args.version}')
    if not os.path.exists(cohort_result_dir):
        os.makedirs(cohort_result_dir)
    print(f'cohort result dir: {cohort_result_dir}')
    
    cohort_data_dir = os.path.join(args.data_dir, args.cpg_type, args.cohort, f'v{args.version}')
    
    for target in all_survival_targets:
        if args.cohort=='TCGA-GBM' and target=='DFI': continue # since most of the samples have null value for both the event indicator and time, only 3 samples are available in (GBM, DFI)
        
        for i in range(args.num_folds):
            current_item = f'fold{i+1}_{target}'
            print(f"===\n{current_item}")
            result_df.loc[current_item]['target'] = target
            result_df.loc[current_item]['fold'] = i+1
           
            train_df_fname = os.path.join(cohort_data_dir, f'{target}_train_fold{i+1}.csv')
            valid_df_fname = os.path.join(cohort_data_dir, f'{target}_valid_fold{i+1}.csv')
            test_df_fname = os.path.join(cohort_data_dir, f'{target}_test_fold{i+1}.csv')
            
            train_df = pd.read_csv(train_df_fname, index_col = 0)
            valid_df = pd.read_csv(valid_df_fname, index_col = 0)
            test_df = pd.read_csv(test_df_fname, index_col = 0)
            
            train_data = dataframe_to_deepsurv_ds(train_df)
            valid_data = dataframe_to_deepsurv_ds(valid_df)
            test_data = dataframe_to_deepsurv_ds(test_df)
            
            hyperparams = {
                'L2_reg': 10.0,
                'batch_norm': True,
                'dropout': 0.4,
                'hidden_layers_sizes': [128, 128], #number of nodes in each hidden layer #length of this array equals the total number of hidden layers inside the DeepSurv model.
                'learning_rate': args.lr,
                'lr_decay': 0.001,
                'momentum': 0.9,
                'n_in': train_data['x'].shape[1],  # input covariate의 개수. 
                'standardize': False
            }
           
            if args.model_load=='Y':
                load_model_name = os.path.join(cohort_result_dir, args.load_model_name)
                load_weight_name = os.path.join(cohort_result_dir, args.load_model_weight_name)
                metrics_basename = os.path.join(cohort_result_dir, args.load_metrics_name)
                model = deep_surv.load_model_from_json(load_model_name, load_weight_name)
                with open(metrics_basename, 'rb') as f:
                    metrics = pickle.load(f)
                f.close()
                
            else:        
                model = deep_surv.DeepSurv(**hyperparams)
                
            experiment_name = f'{args.cpg_type}-{cohort_type}-{current_item}-{args.version}-lr_{args.lr}'
            logdir = './logs/tensorboard/'
            logger = TensorboardLogger(experiment_name, logdir=logdir)

            update_fn=lasagne.updates.nesterov_momentum
            n_epochs = 500 
            start_time = time.time()
            metrics = model.train(train_data, valid_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
            end_time = time.time()

            train_tmp = []      
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            train_c_index = np.max(train_tmp)

            valid_tmp = []
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            valid_c_index = np.max(valid_tmp)
            
            result_df.loc[current_item]['train_c'] = train_c_index
            result_df.loc[current_item]['valid_c'] = valid_c_index
            result_df.loc[current_item]['time'] = end_time - start_time

            print(f'train c-index of the last epoch: {train_c_index}')
            print(f'validation c-index of the last epoch: {valid_c_index}')
            print(f'spent time: {end_time-start_time}')
        
            metrics_basename = f'{current_item}_lr_{args.lr}_metrics'
            metrics_full_name = os.path.join(cohort_result_dir, metrics_basename)
            with open(metrics_full_name, 'wb') as f:
                pickle.dump(metrics, f)
            f.close()
            print(f'metrics filename: {metrics_full_name}')
            result_df.loc[current_item]['metrics'] = metrics_full_name
                
            model_name = os.path.join(cohort_result_dir, f'{current_item}_lr_{args.lr}_model')
            weight_name = os.path.join(cohort_result_dir, f'{current_item}_lr_{args.lr}_model_weights')
            print(f'model: {model_name}')
            print(f'weight: {weight_name}')
            model.save_model(model_name, weight_name)
            result_df.loc[current_item]['model'] = model_name
            result_df.loc[current_item]['weight'] = weight_name            

            plot_log_nll(metrics, cohort_result_dir, f'{current_item}_lr_{args.lr}_train_log_nll.png') #plots negative log likelihood and concordance index
            plot_log_c_index(metrics, cohort_result_dir, f'{current_item}_lr_{args.lr}_train_log_c-index.png') #plots negative log likelihood and concordance index
            result_df.loc[current_item]['train_valid_nll'] = os.path.join(cohort_result_dir, f'{current_item}_lr_{args.lr}_train_log_nll.png')
            result_df.loc[current_item]['train_valid_c_plot'] = os.path.join(cohort_result_dir, f'{current_item}_lr_{args.lr}_train_log_c-index.png')

            test_set_predicted_risk = model.predict_risk(test_data['x'])

            test_df_with_risk_fname = f'{current_item}_lr_{args.lr}_test_df_with_risk.csv'
            test_df_with_risk = test_df.copy()
            test_df_with_risk['risk'] = test_set_predicted_risk
            test_df_with_risk['risk_group'] = pd.qcut(test_df_with_risk['risk'], q = 2, labels = ['Low', 'High'])
            test_df_with_risk.to_csv(os.path.join(cohort_result_dir, test_df_with_risk_fname))
            
            # (Function get_concordance_index) used codes from https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/deepsurv/deep_surv.py#L275
            test_set_c_index = model.get_concordance_index(test_data['x'], test_data['t'], test_data['e'])
            print(f'test set c-index: {test_set_c_index}')
            result_df.loc[current_item]['test_c'] = test_set_c_index

            valid_t, valid_t_pvalue, sig_t, sig_t_pvalue = survival_analysis_single_target(test_df_with_risk, target, cohort_result_dir, 2, f'LogRank_{current_item}_lr_{args.lr}.png', args.cohort)
            result_df.loc[current_item]['LogRank_p'] = valid_t_pvalue[0]
            result_df.loc[current_item]['LogRank_plot'] = os.path.join(cohort_result_dir, f'LogRank_{current_item}_lr_{args.lr}.png')
            del(model)
            result_pvalue = valid_t_pvalue[0]
            print(f'fold {i+1} {target} (lr={args.lr}) LogRank p-value: {result_pvalue}')
            plt.clf()
    result_df.to_csv(os.path.join(cohort_result_dir, f'all_versions_deepsurv_logrank_result_lr_{args.lr}.csv'))
    print(os.path.join(cohort_result_dir, f'all_versions_deepsurv_logrank_result_lr_{args.lr}.csv'))
    
