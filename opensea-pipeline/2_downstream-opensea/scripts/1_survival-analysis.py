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

# lifelines 참고: https://www.notion.so/dohlee/Lifelines-survival-analysis-e929aae590d94037a585b8b1f42bbc2e
# statannot 참고: https://partrita.github.io/posts/statannot/

# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))

CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] 
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2
SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}

def parse_arguments(): #default usage: score=stem-closeness, cohort=each TCGA cohort, s_type=avg-pc1, use_option=part, normalize = N, minmax = Y, w_dir = your working directory.
    args = argparse.ArgumentParser()
    args.add_argument('--score', help = 'Type of score', type = str, required = True)
    args.add_argument('--cohort', help = 'TCGA cohort', type = str, required = True)
    
    # below are all arguments for calculating stem-closeness.
    
    ## score 관련 옵션들 (avg_beta 및 brca subtype으로 돌릴 땐 빼도 됨)
    args.add_argument('--score_type', help = 'which one you will use among PC1 fluctuation or averaged PC1 vector. avg_pc1 or pc1_fluctuation', default = 'avg_pc1', required = False)
    args.add_argument('--usage_option', help = 'use all samples or randomly picked samples. all or part', default = 'part', required = False) # all, part
    args.add_argument('--normalize', help = 'whether you wnat to normalize score2 and score4 or not', default = 'N', required = False)#Y or N
    args.add_argument('--minmax', help = 'use minmax scaling', default = 'N', required = False) #Y or N
    args.add_argument('--matrix_type', help = 'bdm or iebdm', type = str, default = 'iebdm', required = True)
    args.add_argument('--distance', help = 'distance metric. euclidean or cosine-sim', type = str, required = True)
    args.add_argument('--score_result_dir', help = 'directory where the scores of samples are located in.', type = str, required = True)
    args.add_argument('--use_weighted_avg', help = 'whether to use weighted average in computing distance. Y or N (simple_avg)', type = str, required = True)
    args.add_argument('--standardize', help = 'whether to standardize PC1 or not when calculating distance. Y/N', type = str, required = True)
    args.add_argument('--num_chrom', help = 'chr1 to which chromosome you will use.5, 10, 15, 22.', type = int, default = 22, required = False)
    
    args.add_argument('-w_dir', '--working_dir', help = 'working_directory', type = str, required = True)
    args.add_argument('--result_dir', help = 'result directory to save files', type = str, required = True)

    # try both Y and N
    #default: /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/#TCGA-BRCA/

    #score fname example: stem-closeness_euclidean_bdm_pc1-fluctuation_simple-avg_part_normalized.csv
    return args.parse_args()

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


# In[11]:


def save_pickle(data, directory, fname):
    print("fname: {}".format(os.path.join(directory, fname)))
    data.to_pickle(os.path.join(directory, fname))

#이건 빼도 됨
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


if __name__=='__main__':
    
    args = parse_arguments()
    os.chdir(args.working_dir)
    weighted_avg_flag = 'weighted-avg' if args.use_weighted_avg=='Y' else 'simple-avg'
    normalize_yn = '_normalized' if args.normalize=='Y' else ''
    minmax_yn = '_minmax' if args.minmax=='Y' else ''
    standardize_flag = '_standardized' if args.standardize=='Y' else ''
    num_chrom_flag = '_to-chrom-'+str(args.num_chrom) if args.num_chrom < 22 else ''
    
    clinical = pd.read_csv(CLINICAL_FNAME)
    score_fname = args.score+'_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option + normalize_yn+minmax_yn+standardize_flag+num_chrom_flag+'.csv'
    #score filename 예시: 'stem-closeness_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option + normalize_yn+minmax_yn+'.csv')
    survival_target = ['OS', 'DSS', 'DFI', 'PFI'] # survival analysis #these variables are included in 'clinical_cohort'
    result_filenme_default = args.score+'_'+args.distance+'_'+args.matrix_type+'_'+args.score_type+'_'+weighted_avg_flag+'_'+args.usage_option + normalize_yn+minmax_yn+standardize_flag+num_chrom_flag
    
    cohort_dir = os.path.join(args.result_dir, args.cohort, 'kaplan-meier') #결과 저장할 때 씀. 

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir)
    if not os.path.exists(os.path.join(args.result_dir, args.cohort)):
        os.makedirs(os.path.join(args.result_dir, args.cohort))
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    
    print("===\ncohort: {}".format(args.cohort))
    
    cohort = args.cohort
    clinical_cohort = clinical[clinical['type'] == cohort.split('-')[-1]].copy() # clinical data of current cohort

    T, N, S = get_sample_list(args.cohort)
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))
 
    full_score_fname = os.path.join(args.score_result_dir, args.cohort, score_fname)
    print("full_score_fname: {}".format(full_score_fname))

    stem_closeness_df = pd.read_csv(full_score_fname, index_col=0)
    original_index = stem_closeness_df.index.values
    stem_closeness = stem_closeness_df.cos_radian.values.flatten()
    score_df = pd.DataFrame(stem_closeness, index = original_index, columns = ['stem_closeness']).loc[S].copy()

    merged_tumor = merge(clinical_cohort, score_df, T) 
    merged_tumor_fname = 'clinical_'+result_filenme_default+'_merged.csv'
    merged_tumor.to_csv(os.path.join(cohort_dir, merged_tumor_fname))
    print("merged_tumor fname: ", os.path.join(cohort_dir, merged_tumor_fname))

    subtype_dict = boxplots_clinical_variable(merged_tumor, cohort_dir, 'clinical_variables_'+result_filenme_default+'.png', args.cohort) #Tumor samples only

    valid_t, valid_t_pvals, sig_t, sig_t_pvals = survival_analysis(merged_tumor, survival_target, 2, cohort_dir, 3, 
            'survival_analysis_'+result_filenme_default+'.png', args.cohort) 
    #Tumor samples only #생존분석 변수 column 전체가 NaN인 경우도 있음. 

    print("survival analysis: {}".format(os.path.join(cohort_dir, 'survival_analysis_pvalues_'+result_filenme_default+'.npz')))
    np.savez(os.path.join(cohort_dir, 'survival_analysis_pvalues_'+result_filenme_default), 
            valid_target = np.array(valid_t), valid_target_pvals = np.array(valid_t_pvals), significant_target = np.array(sig_t), significant_target_pvals = np.array(sig_t_pvals))

    print("subtype dictionary: {}".format(os.path.join(cohort_dir, 'subtype_dictionary_'+result_filenme_default+'.pickle')))
    with open(os.path.join(cohort_dir, 'subtype_dictionary'+result_filenme_default+'.pickle'), 'wb') as f:
        pickle.dump(subtype_dict, f) 




