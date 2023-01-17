#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# In[3]:


mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)


# In[4]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


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


# In[11]:


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


# In[82]:


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
            
            d['group'] = pd.qcut(d.stem_closeness, q = q, labels = list(range(q)))

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


# In[84]:


score = 'stem-closeness'
cohorts = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
cohorts.remove('TCGA-ESCA')
cohorts.remove('TCGA-HNSC')

cpg_type = 'opensea'

working_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}'

result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/230116_KM_plots'


# In[85]:


os.chdir(working_dir)


# In[86]:


clinical = pd.read_csv(CLINICAL_FNAME)


# In[87]:


all_score_fnames = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)


# In[88]:


survival_target = ['OS', 'DSS', 'DFI', 'PFI'] 


# In[90]:


for cohort in cohorts:
    full_score_fname = all_score_fnames.loc[f'{cohort}'][f'filename_{cpg_type}']
    cohort_dir = os.path.join(result_dir, cohort)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    print("===\ncohort: {}".format(cohort))
    clinical_cohort = clinical[clinical['type'] == cohort.split('-')[-1]].copy() # clinical data of current cohort
    T, N, S = get_sample_list(cohort)
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))
    stem_closeness_df = pd.read_csv(full_score_fname, index_col=0)
    original_index = stem_closeness_df.index.values
    stem_closeness = stem_closeness_df.cos_radian.values.flatten()
    score_df = pd.DataFrame(stem_closeness, index = original_index, columns = ['stem_closeness']).loc[S].copy()
    merged_tumor = merge(clinical_cohort, score_df, T) 
    result_filename_default = full_score_fname.split(cohort)[1].split('.csv')[0].replace('/', '').replace(' ','')
    valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, cohort_dir, 3, 
            'survival_analysis_'+result_filename_default+'.png', cohort) 


# In[ ]:




