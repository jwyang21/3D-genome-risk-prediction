#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# Add score5 to 'survival-analysis_v2'

# statannot 참고: https://partrita.github.io/posts/statannot/

# In[2]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/')


# In[3]:


# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))


# In[4]:


# GLOBAL VARIABLES
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
#CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] + ['chrX','chrY'] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/' 
#/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/TCGA-BLCA/score2345.npz
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2


# In[5]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #for score2
    #s2r, s2d, s4r, s4d: normal, euclidean, all, euclidean
    return parser.parse_args()


# In[6]:


def get_sample_list(directory):
    print("Get_sample_list")
    t = []
    n = []
    for f in os.listdir(directory):
        if 'TCGA' in f and '.npz' in f:
            if f[13:15] in TUMOR_BARCODE:
                t.append(f[:15])
            elif f[13:15] in NORMAL_BARCODE:
                n.append(f[:15])
            else:
                print("Neither tumor nor normal")
    T = list(set(t))
    N = list(set(n))
    S = T + N
    return T, N, S


# In[7]:


def ttest_plot(score_df, T, N, cohort, directory, score, figname):
    # function 'ttest_plot': compare tumor and normal samples' 3dith scores + print mean and std of tumor and normal samples.

    t = score_df.loc[T].values.flatten()
    n = score_df.loc[N].values.flatten()
    #print("Independet t-test between score of tumors and normals: ", end = '')
    #print(ttest_ind(t, n))
    fig = plt.figure(figsize = (14, 7))
    ax1 = fig.add_subplot(121)
    ax1.hist(t)
    ax1.set_xlabel(score)
    ax1.set_ylabel('Frequency')
    ax1.set_title(score+' of tumor samples ({}, n = {})'.format(cohort, len(t)))
    ax2 = fig.add_subplot(122)
    ax2.hist(n)
    ax2.set_xlabel(score)
    ax2.set_ylabel('Frequency')
    ax2.set_title(score+' of normal samples ({}, n = {})'.format(cohort, len(n)))
    plt.tight_layout()

    print("figure file: {}".format(os.path.join(directory, figname)))
    
    #print("mean and std of tumor samples: {} and {}, respectively.".format(np.mean(t), np.std(t)))
    #print("mean and std of normal samples: {} and {}, respectively.".format(np.mean(n), np.std(n)))
    
    plt.savefig(os.path.join(directory, figname))
    
    return np.mean(t), np.std(t), np.mean(n), np.std(n), ttest_ind(t, n)[1]


# In[8]:


def merge(clinical_cohort, score_df, T):
    # function 'merge': merge fire 3dith and clinical data
    
    clinical_cohort.index = clinical_cohort['bcr_patient_barcode'].values.flatten()
    #clinical_cohort.drop(['bcr_patient_barcode'], axis = 1, inplace = True)
    
    #fire = pd.read_pickle(fire_fname)# (22, 1) shape의 df.
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
    #print("final clinical data shape: {}".format(clinical_cohort2.shape))
    #print("final score data shape: {}".format(score_df3.shape))
    
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
    
    return merged, merged_tumor #make sure that 'group' column does not exist in 'merged'


# In[9]:


def save_pickle(data, directory, fname):
    print("fname: {}".format(os.path.join(directory, fname)))
    data.to_pickle(os.path.join(directory, fname))


# In[10]:


def boxplots_clinical_variable(df, directory, figname, cohort): #merged_tumor, SAVEDIR, 'clinical_variables.png', respectively
    
    clinical_variables = ['treatment_outcome_first_course', 'new_tumor_event_site', 'new_tumor_event_type', 'tumor_status',
                          'vital_status', 'histological_grade', 'ajcc_pathologic_tumor_stage', 'race', 'gender']
    
    fig = plt.figure(figsize = (21, 21))
    
    subtype_dict = {}
    for i, c in enumerate(clinical_variables):
        ax_ = fig.add_subplot(3, 3, i+1)
        categories = df[c].unique() # 이 clinical variable의 subtype들
        if '[Not Evaluated]' in categories:
            categories = np.delete(categories, np.where(categories=='[Not Evaluated]')[0][0])
        if '[Unknown]' in categories:
            categories = np.delete(categories, np.where(categories=='[Unknown]')[0][0])
        if '[Not Available]' in categories:
            categories = np.delete(categories, np.where(categories=='[Not Available]')[0][0])
        if '[Not Applicable]' in categories:
            categories = np.delete(categories, np.where(categories=='[Not Applicable]')[0][0])
        if 'nan' in [str(x) for x in categories]:
            categories = np.delete(categories, [str(x) for x in categories].index('nan'))
        
        if len(categories)>0:
            subtype_dict[c] = categories
        
        # 현재 clinical variable의 all-subtype-pairwise independent t-test 진행. 
        all_pairs = [(a, b) for idx, a in enumerate(categories) for b in categories[idx + 1:]]
        
        if len(categories) > 1 and len(categories)<=4: #subgroup 개수가 1 초과 4 이하앤 clinical variable들만 plot
            ax = sns.boxplot(data = df, x = c, y = df.columns[-2], order = categories, ax = ax_)
            ax, test_results = add_stat_annotation(ax, data = df, x = c, y = df.columns[-2], order = categories, 
                                                   box_pairs = all_pairs, test = 't-test_ind', text_format = 'star', loc = 'outside', verbose = 2)
            sns.despine()
    fig.suptitle('Boxplots from each clinical variable ('+cohort+')', fontsize = 15)
    fig.tight_layout()
    #fig.show()
    print("figure file: {}".format(os.path.join(directory, figname)))
    plt.savefig(os.path.join(directory, figname))
    return subtype_dict


# In[31]:


def survival_analysis(df, target, q, directory, fig_width, figname, cohort, score):
    #survival_analysis(merged_tumor, survival_target, 2, SAVEDIR, 3, 'survival_analysis_'+args.score+'.png')
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
            
            if score != 'score5':
                d['group'] = pd.qcut(d.simple_avg, q = q, labels = list(range(q))) # label이 0이면 Low group, 1이면 High group #score2,3,4
            else:
                d['group'] = pd.qcut(d.score5, q = q, labels = list(range(q)))#score5
            
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
            #p = res.p_value
            #target_p.append(p)
            ax.set_title(f'{t}, p = {res.p_value:.2g}', fontsize = 15, pad = 5) #pad: The offset of the title from the top of the axes, in points. Default is None to use rcParams['axes.titlepad'].
            #ax.set_title(f'{target}, p = {p:.2g}', fontsize = 7, pad = 5)
    fig.suptitle('Survival analysis ('+cohort+')', fontsize = 15)
    fig.subplots_adjust(wspace = 0.3)
    fig.tight_layout()
    print("figure file: {}".format(os.path.join(directory, figname)))
    plt.savefig(os.path.join(directory, figname))  
    
    return valid_t, valid_t_pvalue, sig_t, sig_t_pvalue


if __name__=='__main__':
    
    args = parse_arguments()
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    print("cohort: {}".format(args.cohort))
    
    clinical = pd.read_csv(CLINICAL_FNAME)
    
    cohort = args.cohort
    clinical_cohort = clinical[clinical['type'] == cohort.split('-')[-1]].copy() # clinical data of current cohort
    #del(cohort)

    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))


    if args.score == 'score2':
        SCORE_FNAME = os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/', 
                                   args.cohort, args.score+'_'+args.reference+'_'+args.distance+'_simple_avg.pickle')
        score_df = pd.read_pickle(SCORE_FNAME)
        tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval = ttest_plot(score_df, T, N, args.cohort, SAVEDIR, args.score,
                                                                                            'tumor_normal_'+args.score+'_'+args.reference+'_'+args.distance+'_histogram.png') 
    
        merged, merged_tumor = merge(clinical_cohort, score_df, T) 
    
        save_pickle(merged, SAVEDIR, 'clinical_'+args.score+'_'+args.reference+'_'+args.distance+'_merged.pickle')

        save_pickle(merged_tumor, SAVEDIR, 'clinical_'+args.score+'_'+args.reference+'_'+args.distance+'_merged_TumorOnly.pickle')

        subtype_dict = boxplots_clinical_variable(merged_tumor, SAVEDIR, 
                                                  'clinical_variables_'+args.score+'_'+args.reference+'_'+args.distance+'.png', args.cohort) #Tumor samples only

        survival_target = ['OS', 'DSS', 'DFI', 'PFI'] # survival analysis #these variables are included in 'clinical_cohort'

        valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, SAVEDIR, 3, 
                                                                       'survival_analysis_'+args.score+'_'+args.reference+'_'+args.distance+'.png', 
                                                                       args.cohort) #Tumor samples only #생존분석 변수 column 전체가 NaN인 경우도 있음. 

        print("survival analysis: {}".format(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_survival_analysis_pvalues.npz')))
        np.savez(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_survival_analysis_pvalues'), 
                 valid_target = np.array(valid_t), valid_target_pvals = np.array(valid_t_pvals), # 생존분석 타겟들 중에서, 전체 값이 NaN이지 않은 항목들.
                 significant_target = np.array(sig_t), significant_target_pvals = np.array(sig_t_pvals))

        print("tumor v.s. normal ttest: {}".format(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_tumor_normal_ttest.npz')))
        np.savez(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_tumor_normal_ttest'), 
                 values = np.array([tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval]), 
                 items = np.array(['tumor_mean', 'tumor_std', 'normal_mean', 'normal_std', 'tumor_normal_ttest_pval']))

        print("subtype dictionary: {}".format(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_subtype_dictionary.pickle')))
        with open(os.path.join(SAVEDIR, args.score+'_'+args.reference+'_'+args.distance+'_subtype_dictionary.pickle'), 'wb') as f:
            pickle.dump(subtype_dict, f) 
            
    elif args.score == 'score3':
        SCORE_FNAME = os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/', args.cohort, args.score+'_simple_avg.pickle')
        score_df = pd.read_pickle(SCORE_FNAME)
        
        tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval = ttest_plot(score_df, T, N, args.cohort, SAVEDIR, args.score,
                                                                                            'tumor_normal_'+args.score+'_histogram.png') 
    
        merged, merged_tumor = merge(clinical_cohort, score_df, T) 
        
        save_pickle(merged, SAVEDIR, 'clinical_'+args.score+'_merged.pickle')

        save_pickle(merged_tumor, SAVEDIR, 'clinical_'+args.score+'_merged_TumorOnly.pickle')

        subtype_dict = boxplots_clinical_variable(merged_tumor, SAVEDIR, 'clinical_variables_'+args.score+'.png', args.cohort) #Tumor samples only

        survival_target = ['OS', 'DSS', 'DFI', 'PFI'] # survival analysis #these variables are included in 'clinical_cohort'

        valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, SAVEDIR, 3, 'survival_analysis_'+args.score+'.png', args.cohort) #Tumor samples only #생존분석 변수 column 전체가 NaN인 경우도 있음. 
        
        print("survival analysis: {}".format(os.path.join(SAVEDIR, 'survival_analysis_pvalues.npz')))
        np.savez(os.path.join(SAVEDIR, 'survival_analysis_pvalues'), valid_target = np.array(valid_t), 
                 valid_target_pvals = np.array(valid_t_pvals), significant_target = np.array(sig_t), significant_target_pvals = np.array(sig_t_pvals))

        print("tumor v.s. normal ttest: {}".format(os.path.join(SAVEDIR, 'tumor_normal_ttest.npz')))
        np.savez(os.path.join(SAVEDIR, 'tumor_normal_ttest'), values = np.array([tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval]), 
                 items = np.array(['tumor_mean', 'tumor_std', 'normal_mean', 'normal_std', 'tumor_normal_ttest_pval']))

        print("subtype dictionary: {}".format(os.path.join(SAVEDIR, 'subtype_dictionary.pickle')))
        with open(os.path.join(SAVEDIR, 'subtype_dictionary.pickle'), 'wb') as f:
            pickle.dump(subtype_dict, f)
    
    elif args.score == 'score5':
        
        SCORE_FNAME = os.path.join(SCORE_DIR, args.cohort,'score2345.npz')
        samples=np.load(SCORE_FNAME, allow_pickle = True)['rownames']
        score = np.load(SCORE_FNAME, allow_pickle = True)['score5']
        score_df = pd.DataFrame(score, index = samples, columns = ['score5']).loc[S]

        tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval = ttest_plot(score_df, T, N, args.cohort, SAVEDIR, args.score,
                                                                                            'tumor_normal_'+args.score+'_histogram.png') 

        merged, merged_tumor = merge(clinical_cohort, score_df, T) 

        save_pickle(merged, SAVEDIR, 'clinical_'+args.score+'_merged.pickle')

        save_pickle(merged_tumor, SAVEDIR, 'clinical_'+args.score+'_merged_TumorOnly.pickle')

        subtype_dict = boxplots_clinical_variable(merged_tumor, SAVEDIR, 'clinical_variables_'+args.score+'.png', args.cohort) #Tumor samples only

        survival_target = ['OS', 'DSS', 'DFI', 'PFI'] # survival analysis #these variables are included in 'clinical_cohort'

        valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, 
                                                                       SAVEDIR, 3, 'survival_analysis_'+args.score+'.png', args.cohort, args.score) 
        #Tumor samples only #생존분석 변수 column 전체가 NaN인 경우도 있음. 

        print("survival analysis: {}".format(os.path.join(SAVEDIR, 'survival_analysis_pvalues_'+args.score+'.npz')))
        np.savez(os.path.join(SAVEDIR, 'survival_analysis_pvalues_'+args.score), valid_target = np.array(valid_t), 
                 valid_target_pvals = np.array(valid_t_pvals), significant_target = np.array(sig_t), significant_target_pvals = np.array(sig_t_pvals))

        print("tumor v.s. normal ttest: {}".format(os.path.join(SAVEDIR, 'tumor_normal_ttest_'+args.score+'.npz')))
        np.savez(os.path.join(SAVEDIR, 'tumor_normal_ttest_'+args.score), values = np.array([tumor_mean, tumor_std, normal_mean, normal_std, tumor_normal_ttest_pval]), 
                 items = np.array(['tumor_mean', 'tumor_std', 'normal_mean', 'normal_std', 'tumor_normal_ttest_pval']))

        print("subtype dictionary: {}".format(os.path.join(SAVEDIR, 'subtype_dictionary_'+args.score+'.pickle')))
        with open(os.path.join(SAVEDIR, 'subtype_dictionary_'+args.score+'.pickle'), 'wb') as f:
            pickle.dump(subtype_dict, f)
        
        print("----\n")
        
    else:
        pass

