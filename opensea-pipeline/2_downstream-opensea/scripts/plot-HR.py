import pandas as pd
import numpy as np
import os
import glob

import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

mpl.rcParams['figure.dpi'] = 300 
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

def parse_argument():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'opensea, island, shelf_shore', type = str, required = True)
    args.add_argument('--cox_version', help = 'which type of cox version you want to plot', type = int, required = True)
    return args.parse_args()

def draw_cox_hr(tmp, save_dir, fig_name, cox_version, fig_right=3):
    fig = plt.figure(figsize = (fig_right, 2))
    ax = fig.add_subplot(232)
    yticks = []
    all_covariates = []
    all_heights = []
    ax.set_xlim([0, fig_right])
    x_mins = []
    x_maxs = []
    for i in range(tmp.shape[0]):
        covariate = tmp.index.values[i]
        all_covariates.append(covariate)
        height = 0.1 * (i+1)
        all_heights.append(height)
        yticks.append(tmp.index.values[i])
        x_min = float(tmp['exp(coef) lower 95%'][covariate])
        x_max = float(tmp['exp(coef) upper 95%'][covariate])
        x_mins.append(x_min)
        x_maxs.append(x_max)
    assert len(all_heights) == len(x_mins) and len(x_mins) == len(x_maxs)
    for i in range(len(x_mins)):
        x_min = x_mins[i]
        x_max = x_maxs[i]
        height = all_heights[i]
        ax.hlines(height, xmin = x_min, xmax = x_max, color='lightgray', linestyle='-', linewidth=2)
    ax.axvline(x=1, ymin=0, ymax=1, color = 'black', linestyle = '--')
    ax.set_xlabel('Hazard Ratio (HR)', fontsize = 7)
    ax.set_yticks(all_heights)
    ax.set_yticklabels(all_covariates)
    for i in range(tmp.shape[0]):    
        covariate = tmp.index.values[i]
        height = 0.1 * (i+1)
        hr = float(tmp['exp(coef)'][covariate])
        x_min = float(tmp['exp(coef) lower 95%'][covariate])
        x_max = float(tmp['exp(coef) upper 95%'][covariate])
        p = float(tmp['p'][covariate])
        ax.plot(hr, height, 'rD', ms = 2)
        str_ = f'HR = {round(hr, 2)} ({round(x_min, 2)}, {round(x_max, 2)}), p = {round(p, 2)}'
        ax.text(fig_right+0.2, height-0.01, str_, fontsize = 7)
    mb_length = 0.25
    ax.plot([0.57, 0.57+mb_length-0.035], [1.2, 1.2], linestyle = '--', lw=1.5, c='k', clip_on=False, transform=ax.transAxes)
    ax.text(0.57 + mb_length , 1.2, 'HR = 1', transform=ax.transAxes, ha='left', va='center', fontsize = 5)
    ax.plot([0.57, 0.57+mb_length-0.035], [1.1, 1.1], linestyle = '-', lw=1.5, color = 'lightgray', clip_on=False, transform=ax.transAxes)
    ax.text(0.57 + mb_length , 1.1, '95% CI', transform=ax.transAxes, ha='left', va='center', fontsize = 5)
    plt.savefig(os.path.join(save_dir, fig_name))

def draw_log_cox_hr(tmp, save_dir, fig_name, cox_version, fig_right=3):
    fig = plt.figure(figsize = (fig_right, 2))
    ax = fig.add_subplot(232)
    yticks = []
    all_covariates = []
    all_heights = []
    ax.set_xlim([-1, fig_right])
    x_mins = []
    x_maxs = []
    for i in range(tmp.shape[0]):
        covariate = tmp.index.values[i]
        all_covariates.append(covariate)
        height = 0.1 * (i+1)
        all_heights.append(height)
        yticks.append(tmp.index.values[i])
        x_min = float(tmp['coef lower 95%'][covariate])
        x_max = float(tmp['coef upper 95%'][covariate])
        x_mins.append(x_min)
        x_maxs.append(x_max)
    assert len(all_heights) == len(x_mins) and len(x_mins) == len(x_maxs)
    for i in range(len(x_mins)):
        x_min = x_mins[i]
        x_max = x_maxs[i]
        height = all_heights[i]
        ax.hlines(height, xmin = x_min, xmax = x_max, color='lightgray', linestyle='-', linewidth=2)
    ax.axvline(x=0, ymin=0, ymax=1, color = 'black', linestyle = '--')
    ax.set_xlabel('log(HR)', fontsize = 7)
    ax.set_yticks(all_heights)
    ax.set_yticklabels(all_covariates)
    for i in range(tmp.shape[0]):    
        covariate = tmp.index.values[i]
        height = 0.1 * (i+1)
        hr = float(tmp['coef'][covariate])
        x_min = float(tmp['coef lower 95%'][covariate])
        x_max = float(tmp['coef upper 95%'][covariate])
        p = float(tmp['p'][covariate])
        #ax.hlines(height, xmin = x_min, xmax = x_max, color='lightgray', linestyle='-', linewidth=2)
        ax.plot(hr, height, 'rD', ms = 2)
        str_ = f'HR = {round(hr, 2)} ({round(x_min, 2)}, {round(x_max, 2)}), p = {round(p, 2)}'
        ax.text(fig_right+0.2, height-0.01, str_, fontsize = 5)
    mb_length = 0.25
    ax.plot([0.57, 0.57+mb_length-0.035], [1.2, 1.2], linestyle = '--', lw=1.5, c='k', clip_on=False, transform=ax.transAxes)
    ax.text(0.57 + mb_length , 1.2, 'log(HR) = 0', transform=ax.transAxes, ha='left', va='center', fontsize = 5)
    ax.plot([0.57, 0.57+mb_length-0.035], [1.1, 1.1], linestyle = '-', lw=1.5, color = 'lightgray', clip_on=False, transform=ax.transAxes)
    ax.text(0.57 + mb_length , 1.1, '95% CI', transform=ax.transAxes, ha='left', va='center', fontsize = 5)
    plt.savefig(os.path.join(save_dir, fig_name))
    print(f'figure: {os.path.join(save_dir, fig_name)}')

if __name__=='__main__':
    args = parse_arguments()
    
    result_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/2_downstream-{args.cpg_type}/result'

    cohorts = []
    for f in os.listdir(result_dir):
        if 'TCGA' in f:
            cohorts.append(os.path.basename(f))
            
    if args.cox_version in [4, 5]:
        fig_right = 6
    else:
        fig_right = 3
    '''
    for cohort in cohorts:
        print(f"---\n{cohort}")
        cohort_dir = os.path.join(result_dir, cohort, f'cox_v{args.cox_version}')
        for f in glob.glob(os.path.join(cohort_dir, f'*.csv')):
            print(f'input: {f}')
            f_basename = os.path.basename(f)
            fig_name = f_basename.split('.csv')[0]+'-v2-log.png'
            tmp = pd.read_csv(f, index_col = 0)
            draw_log_cox_hr(tmp, cohort_dir, fig_name, args.cox_version)
     '''       
    for cohort in cohorts:
        print(f"---\n{cohort}")
        cohort_dir = os.path.join(result_dir, cohort, f'cox_v{args.cox_version}')
        for f in glob.glob(os.path.join(cohort_dir, f'*.csv')):
            print(f)
            f_basename = os.path.basename(f)
            fig_name = f_basename.split('.csv')[0]+'-v2.png'
            tmp = pd.read_csv(f, index_col = 0)
            draw_cox_hr(tmp, cohort_dir, fig_name, args.cox_version, fig_right)
