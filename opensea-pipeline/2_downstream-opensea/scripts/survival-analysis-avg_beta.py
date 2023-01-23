#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from collections import defaultdict
import lifelines
import matplotlib as mpl
from statannot import add_stat_annotation
import argparse
import pickle
from sklearn.preprocessing import MinMaxScaler


# In[2]:


# lifelines 참고: https://www.notion.so/dohlee/Lifelines-survival-analysis-e929aae590d94037a585b8b1f42bbc2e
# statannot 참고: https://partrita.github.io/posts/statannot/

# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))


# In[3]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


# In[4]:


def parse_arguments(): #default usage: score=stem-closeness, cohort=each TCGA cohort, s_type=avg-pc1, use_option=part, normalize = N, minmax = Y, w_dir = your working directory.
    args = argparse.ArgumentParser()
    args.add_argument('--score_name', help = 'name of score. It should be exactly same with the name of column where the score-of-interest is saved in.', type = str, required = True)
    args.add_argument('--cohort', help = 'TCGA cohort', type = str, required = True)
    
    # below are all arguments for calculating stem-closeness.
    
    ## score 관련 옵션들 (avg_beta 및 brca subtype으로 돌릴 땐 빼도 됨)
    #args.add_argument('--score_type', help = 'which one you will use among PC1 fluctuation or averaged PC1 vector. avg_pc1 or pc1_fluctuation', default = 'avg_pc1', required = False)
    #args.add_argument('--usage_option', help = 'use all samples or randomly picked samples. all or part', default = 'part', required = False) # all, part
    #args.add_argument('--normalize', help = 'whether you wnat to normalize score2 and score4 or not', default = 'N', required = False)#Y or N
    #args.add_argument('--minmax', help = 'use minmax scaling', default = 'N', required = False) #Y or N
    #args.add_argument('--matrix_type', help = 'bdm or iebdm', type = str, default = 'iebdm', required = True)
    #args.add_argument('--distance', help = 'distance metric. euclidean or cosine-sim', type = str, required = True)
    args.add_argument('--score_result_dir', help = 'directory where the scores of samples are located in.', type = str, required = True)
    #example: /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result
    
    #args.add_argument('--cpg_type', type = str, required = True)
    
    args.add_argument('--score_fname', help = 'score filename', type = str, required = True)
    #opensea_tumors_avg_beta.csv
    
    #args.add_argument('--use_weighted_avg', help = 'whether to use weighted average in computing distance. Y or N (simple_avg)', type = str, required = True)
    #args.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
    #args.add_argument('--num_chrom', help = 'chr1 to which chromosome you will use.5, 10, 15, 22.', type = int, default = 22, required = False)
    
    args.add_argument('-w_dir', '--working_dir', help = 'working_directory', type = str, required = True)
    args.add_argument('--result_dir', help = 'result directory to save files', type = str, required = True)
    args.add_argument('--lowest_save_dir', help = 'name of the lowest subdirectory where the result file is saved in', type = str, required = True)

    # try both Y and N
    #default: /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/#TCGA-BRCA/

    #score fname example: stem-closeness_euclidean_bdm_pc1-fluctuation_simple-avg_part_normalized.csv
    return args.parse_args()


# In[5]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = []
    else: #TCGA cohort
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S


# In[6]:


def merge(clinical_cohort, score_df, T):
    # function 'merge': merge score and clinical data    
    clinical_cohort.index = clinical_cohort['bcr_patient_barcode'].values.flatten()
    #clinical_cohort.drop(['bcr_patient_barcode'], axis = 1, inplace = True)

    print("score shape: {}".format(score_df.shape))
    score_df['barcode'] = [x[:12] for x in score_df.index.values]
    score_df['full_barcode'] = score_df.index.values
    score_df['sample_type'] = ['Tumor' if score_df.index.values[i] in T else 'Normal' for i in range(score_df.shape[0])]
    score_df2 = score_df.drop_duplicates(subset = ['barcode']).copy() 
    # barcode가 중복될 경우 keep = 'first'. 나머진 drop. # 동일 샘플에서 여러번 데이터가 수집된 경우가 여기에 해당됨
    
    score_df2.index = score_df2.barcode.values
    score_df2.drop(['barcode'], axis = 1, inplace=True)
    
    intersection = np.intersect1d(clinical_cohort.index.values, score_df2.index.values)
    intersection_tumor_only = np.intersect1d(clinical_cohort.index.values, score_df2[score_df2['sample_type']=='Tumor'].index.values)
    score_df2.drop(['sample_type'], axis = 1, inplace=True)
    
    # clinical과 score dataframe (현재 코호트의 전체 샘플들) 에서 겹치는 애들만 남기기
    clinical_cohort2 = clinical_cohort.loc[intersection].copy().sort_index() # 전체 코호트 샘플들과 겹치는 애들만 남김.
    score_df3 = score_df2.loc[intersection].copy().sort_index() # clinical과 겹치는 애들만 남김
    
    merged = pd.merge(clinical_cohort2, score_df3, left_index=True, right_index=True)

    # clinical 과 score dataframe 에서 겹치는 'tumor' 샘플들만 남기기
    clinical_cohort3 = clinical_cohort.loc[intersection_tumor_only].copy().sort_index()
    score_df4 = score_df2.loc[intersection_tumor_only].copy().sort_index()
    print("final clinical data shape (tumor only): {}".format(clinical_cohort3.shape))
    print("final score data shape (tumor only): {}".format(score_df4.shape))
    
    merged_tumor = pd.merge(clinical_cohort3, score_df4, left_index=True, right_index=True)
    
    if 'sample_type' in merged.columns:
        print('Error! type is still in merged.columns')
    
    if 'sample_type' in merged_tumor.columns:
        print('Error! type is still in merged_tumor.columns')
    
    return merged_tumor #make sure that 'group' column does not exist in 'merged'


# In[39]:

'''
def survival_analysis(df, target, q, directory, fig_width, figname, cohort):
    # function 'survival_analysis': plot survival analysis results and save resulting figure.
    
    #fig_width = 4
    fig = plt.figure(figsize = (4 * fig_width, fig_width)) #각 subplot은 정사각형 모양.
    
    valid_t = []
    valid_t_pvalue = []
    sig_t = []
    sig_t_pvalue = []
    
    for i, t in enumerate(target): # Iterate for all targets. 
        d = df[df[f'{t}.time'].notnull() & df[f'{t}'].notnull()].copy()
        if d.shape[0] !=0:
            valid_t.append(t)
            
            d['group'] = pd.qcut(d.score, q = q, labels = list(range(q)))

            print("target = {}, num_samples = {}".format(t, d.shape[0]))

            groups = list(range(q)) # 총 q개만큼의 그룹이 있음.

            ax = fig.add_subplot(1, 4, i+1)

            for group in groups:#group 0 (low), group 1(High)
                # T: event가 일어나기까지의 시간.
                # E: event의 발생 여부.
                T = d[d.group==group][f'{t}.time'].values
                E = d[d.group==group][f'{t}'].astype(bool).values #0, 1을 False, True로 바꿈.

                kmf = lifelines.KaplanMeierFitter()
                kmf.fit(T, E, label = ['Low', 'High'][group])

                kmf.plot_survival_function(ax = ax, ci_show = False, linewidth = 3, xticks = [0, 2000], yticks = [0.2, 0.4, 0.6, 0.8, 1.0]) #ci_show: show confidence intervals if True.
                #kmf.plot_survival_function(ax = ax, ci_show = False) #ci: confidence interval
            ax.get_xaxis().set_visible(True)
            ax.grid(False)
            ax.set_facecolor('white')
            ax.spines['top'].set_color('black')
            ax.spines['bottom'].set_color('black')
            ax.spines['left'].set_color('black')
            ax.spines['right'].set_color('black')

            res = lifelines.statistics.logrank_test( #input: T_group0, T_group1, E_group0, E_group1
                d[d.group==0][f'{t}.time'].values,
                d[d.group==q-1][f'{t}.time'].values,
                d[d.group==0][f'{t}'].values,
                d[d.group==q-1][f'{t}'].values
            )
            ax.legend(frameon = False)
            valid_t_pvalue.append(res.p_value)
            if res.p_value < P_THRESHOLD:
                sig_t.append(t)
                sig_t_pvalue.append(res.p_value)
            ax.set_title(f'{t}, p = {res.p_value:.2g}', fontsize = 15, pad = 5) #pad: The offset of the title from the top of the axes, in points. Default is None to use rcParams['axes.titlepad'].
    fig.suptitle('Survival analysis ('+cohort+')', fontsize = 15)
    fig.subplots_adjust(wspace = 0.3)
    fig.tight_layout()
    print("figure file: {}".format(os.path.join(directory, figname)))
    plt.savefig(os.path.join(directory, figname))  
    
    return valid_t, valid_t_pvalue, sig_t, sig_t_pvalue
'''
def survival_analysis(df, target, q, directory, fig_width, figname, cohort):
    # function 'survival_analysis': plot survival analysis results and save resulting figure.
    
    #fig_width = 4
    fig = plt.figure(figsize = (4 * fig_width, fig_width)) #각 subplot은 정사각형 모양.
    
    valid_t = []
    valid_t_pvalue = []
    sig_t = []
    sig_t_pvalue = []
    
    for i, t in enumerate(target): # Iterate for all targets. 
        d = df[df[f'{t}.time'].notnull() & df[f'{t}'].notnull()].copy()
        if d.shape[0] !=0:
            valid_t.append(t)
            
            d['group'] = pd.qcut(d.score, q = q, labels = list(range(q)))

            print("target = {}, num_samples = {}".format(t, d.shape[0]))

            groups = list(range(q)) # 총 q개만큼의 그룹이 있음.

            ax = fig.add_subplot(1, 4, i+1)

            for group in groups:#group 0 (low), group 1(High)
                # T: event가 일어나기까지의 시간.
                # E: event의 발생 여부.
                T = d[d.group==group][f'{t}.time'].values
                E = d[d.group==group][f'{t}'].astype(bool).values #0, 1을 False, True로 바꿈.

                kmf = lifelines.KaplanMeierFitter()
                kmf.fit(T, E, label = ['Low', 'High'][group])

                kmf.plot_survival_function(ax = ax, ci_show = False, linewidth = 3, xticks = [0, 2000], yticks = [0.2, 0.4, 0.6, 0.8, 1.0]) #ci_show: show confidence intervals if True.
                #kmf.plot_survival_function(ax = ax, ci_show = False) #ci: confidence interval
            ax.get_xaxis().set_visible(True)
            ax.grid(False)
            ax.set_facecolor('white')
            ax.spines['top'].set_color('black')
            ax.spines['bottom'].set_color('black')
            ax.spines['left'].set_color('black')
            ax.spines['right'].set_color('black')

            res = lifelines.statistics.logrank_test( #input: T_group0, T_group1, E_group0, E_group1
                d[d.group==0][f'{t}.time'].values,
                d[d.group==q-1][f'{t}.time'].values,
                d[d.group==0][f'{t}'].values,
                d[d.group==q-1][f'{t}'].values
            )
            
            num_two_thousands = int((d[f'{t}.time'].dropna().values.max() // 2000) + 1)
            print('num_thousands: ', num_two_thousands)
            xticks_ = [i * 2000 for i in range(num_two_thousands)]
            ax.set_xticks(xticks_)
            
            
            ax.legend(frameon = False)
            valid_t_pvalue.append(res.p_value)
            if res.p_value < P_THRESHOLD:
                sig_t.append(t)
                sig_t_pvalue.append(res.p_value)
            ax.set_title(f'{t}, p = {res.p_value:.2g}', fontsize = 13, pad = 5) #pad: The offset of the title from the top of the axes, in points. Default is None to use rcParams['axes.titlepad'].
            ax.set_xlabel('Days', fontsize = 7)
            ax.set_ylabel('Survival Probability', fontsize = 7)
            #ax.set_xlim([0, 10000])
    fig.suptitle('Survival analysis ('+cohort+')', fontsize = 15)
    fig.subplots_adjust(wspace = 0.3)
    fig.tight_layout()
    print("figure file: {}".format(os.path.join(directory, figname)))
    plt.savefig(os.path.join(directory, figname))  
    
    return valid_t, valid_t_pvalue, sig_t, sig_t_pvalue



if __name__=='__main__':
    
    args = parse_arguments()
    #working_dir = f''#working directory
    #result_dir = f''# 결과가 저장될 최종 디렉토리의 상위 디렉토리
    #score_result_dir = f''# score file들이 있는 디렉토리
    os.chdir(args.working_dir)
    
    clinical = pd.read_csv(CLINICAL_FNAME)

    survival_target = ['OS', 'DSS', 'DFI', 'PFI'] # survival analysis #these variables are included in 'clinical_cohort'

    #score_fname = os.path.join(args.score_dir, args.cohort, args.score_fname)
    #score_fname: opensea_tumors_avg_beta.csv
    
    cohort_dir = os.path.join(args.result_dir, args.cohort, args.lowest_save_dir) #결과 저장할 때 씀. 

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
    if not os.path.exists(os.path.join(args.result_dir, args.cohort)):
        os.makedirs(os.path.join(args.result_dir, args.cohort))
    if not os.path.exists(os.path.join(args.result_dir, args.cohort, args.lowest_save_dir)):
        os.makedirs(os.path.join(args.result_dir, args.cohort, args.lowest_save_dir))
    
    print("===\ncohort: {}".format(args.cohort))
    
    cohort = args.cohort
    clinical_cohort = clinical[clinical['type'] == cohort.split('-')[-1]].copy() # clinical data of current cohort

    T, N, S = get_sample_list(args.cohort)
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))
 
    full_score_fname = os.path.join(args.score_result_dir, args.cohort, args.score_fname)

    score_df = pd.read_csv(full_score_fname, index_col=0)
    original_index = score_df.index.values
    score = score_df[f'{args.score_name}'].values.flatten()
    
    # score_df 에는 target score 한 개만 있어야 함 (1 column). 
    # score_df의 이름은 'score'로 할 것.
    score_df = pd.DataFrame(score, index = original_index, columns = ['score']).loc[T].copy()

    merged_tumor = merge(clinical_cohort, score_df, T) 
    merged_tumor_fname = f'clinical_{args.score_name}_merged.csv'
    merged_tumor.to_csv(os.path.join(cohort_dir, merged_tumor_fname))
    print("merged_tumor fname: ", os.path.join(cohort_dir, merged_tumor_fname))

    valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, cohort_dir, 3, 
            f'survival_analysis_{args.score_name}.png', args.cohort) 
    
    #Tumor samples only #생존분석 변수 column 전체가 NaN인 경우도 있음. 
    print("survival analysis: {}".format(os.path.join(cohort_dir, f'survival_analysis_pvalues_{args.score_name}.npz')))
    np.savez(os.path.join(cohort_dir, f'survival_analysis_pvalues_{args.score_name}'), 
            valid_target = np.array(valid_t), valid_target_pvals = np.array(valid_t_pvals), significant_target = np.array(sig_t), significant_target_pvals = np.array(sig_t_pvals))

