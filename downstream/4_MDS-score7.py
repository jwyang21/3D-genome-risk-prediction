#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import sys
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from sklearn.manifold import MDS
from sklearn.utils.validation import check_symmetric
from collections import Counter
import matplotlib as mpl
import argparse


# - sample type에 따라 라벨링 (모든 샘플 사용)
# - 각 clinical variable의 subtype에 따라 라벨링 (tumor sample only)
# - score_group에 따라 라벨링 (tumor sample only)

# score7 all-sample-pairwise distance file: 
# - /data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/TCGA-BLCA/score7_sample_sample_euclidean.csv

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
CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
#SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/' #scores other than score7
SCORE_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta/result/'#TCGA-BLCA/score2-score4-cosine.csv'

#score3 pickle file example: /data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/{cohort}/{cohort}_score3.pickle

#jsd file example: /data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/TCGA-ACC/score3_sample_sample_jsd.csv
#FIRE_DIR = '/data/project/3dith/data/'
#FIRE_AB_FNAME = 'fire_a_b_labels.csv'
#FIRE_PC1_FNAME = 'fire_pc1.csv'
#COHORT_TCGA_FIRE_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA_FIRE_intersection.csv'
#COHORT_TCGA_FIRE = pd.read_csv(COHORT_TCGA_FIRE_FNAME)
RANDOM_STATE = 42
COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
CLINICAL_FNAME = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv'
P_THRESHOLD = 5e-2
SUBTYPE_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA_subtypes.csv'
SUBTYPE = pd.read_csv(SUBTYPE_FNAME)
DISTANCE_MATRIX_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/'
CLINICAL_VARIABLES = ['treatment_outcome_first_course', 'new_tumor_event_site', 'new_tumor_event_type', 'tumor_status', 'vital_status', 'histological_grade', 'ajcc_pathologic_tumor_stage', 'race', 'gender']
EXCLUDE_SUBTYPES = ['[Not Evaluated]', '[Unknown]', '[Not Available]', '[Not Applicable]', 'nan']


# In[5]:


def parse_arguments(): 
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True)
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    parser.add_argument('-d', '--distance', help = 'Distance metric used for computing how similar/different two samples are', type = str, default = 'euclidean', required = False)
    return parser.parse_args()


# In[26]:


def get_sample_list(directory):
    #print("Get_sample_list")
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
    print("Tumor: {}, Normal: {}, All: {}".format(len(T), len(N), len(S)))
    return T, N, S



def tumor_normal_mds_score7(cohort, score, metric, T, N, S, savedir, figure_name, fig_size):
    # import all sample-sample pairwise distance matrix
    #fname = '/data/project/jeewon/research/3D-ITH/binned_diff/score2_pc1/fire_compare_result/'+cohort+'/sample_sample_pc1_dist_'+chrom+'.pickle' #for score 2
    distance_matrix_fname = os.path.join(DISTANCE_MATRIX_DIR, cohort, score+'_sample_sample_'+metric+'.csv') # all sample MDS

    distance = pd.read_csv(distance_matrix_fname, index_col = 0) # All samples
    
    if metric == 'cosine-similarity':
        distance2 = 1/np.exp(distance)
    else:
        distance2 = distance.copy()
    
    
    mds = MDS(dissimilarity = 'precomputed', random_state = RANDOM_STATE)
    embedding = mds.fit_transform(distance2.values)

    # save embedding
    embedding_df = pd.DataFrame(embedding, index = distance2.index.values, columns = ['0', '1'])
    #embedding_df_fname = 'all_Samples_MDS_'+score+'_'+metric+'.pickle'
    #embedding_df.to_pickle(os.path.join(savedir, embedding_df_fname))
    fig = plt.figure(figsize = (fig_size, fig_size))
    ax = fig.add_subplot(111)
    ax.scatter(embedding_df.loc[T].values[:,0], embedding_df.loc[T].values[:,1], c = COLORS[0], label = 'Tumor') #tumor
    ax.scatter(embedding_df.loc[N].values[:,0], embedding_df.loc[N].values[:,1], c = COLORS[1], label = 'Normal') #normal
    ax.set_title(cohort+' (tumor: '+str(len(T))+', normal: '+str(len(N))+', '+score+')', fontsize = 15)
    ax.legend()
    fig.tight_layout()
    
    #plt.show()
    print("figure file: {}".format(os.path.join(savedir, figure_name)))
    plt.savefig(os.path.join(savedir, figure_name))
    return embedding_df




def clinical_variable_mds_score7(cohort, T, savedir, figure_name, fig_size, metric, score):
    merged_tumor_fname = os.path.join(os.getcwd(), 'result', cohort, 'clinical_'+score+'_merged_TumorOnly.pickle')
    merged_tumor = pd.read_pickle(merged_tumor_fname)
    merged_tumor.index = merged_tumor.full_barcode.values.flatten()
    
    # each subplot for each clinical variable
    fig = plt.figure(figsize = (3*fig_size, 3*fig_size)) #number of clinical variables == 9
    
    distance_matrix_fname = os.path.join(DISTANCE_MATRIX_DIR, cohort, score+'_sample_sample_'+metric+'.csv')
    distance = pd.read_csv(distance_matrix_fname, index_col = 0) # All samples
    
    if metric == 'cosine-similarity':
        distance2 = 1/np.exp(distance)
    else:
        distance2 = distance.copy()
    
    tumor_distance = distance2.loc[T, T]
    mds = MDS(dissimilarity = 'precomputed', random_state = RANDOM_STATE)
    tumor_embedding = mds.fit_transform(tumor_distance.values)
    tumor_embedding_df = pd.DataFrame(tumor_embedding, index = T, columns = ['0', '1'])
    
    # full barcode 말고 bcr_patient_barcode로 줄였을 때 중복되는 샘플이 존재하는 경우가 있음. (동일 샘플에게서 여러번 데이터를 얻었을 때.) -> drop duplicates
    bcr_barcodes = [x[:12] for x in tumor_embedding_df.index.values]
    tumor_embedding_df['full_barcode'] = tumor_embedding_df.index.values
    tumor_embedding_df['bcr_barcodes'] = bcr_barcodes
    tumor_embedding_df.index = bcr_barcodes
    tumor_embedding_df2 = tumor_embedding_df[~tumor_embedding_df.index.duplicated(keep = 'last')]
    
    #print("tumor-tumor-pairwise distance: {}".format(tumor_distance.shape))
    
    with open(os.path.join(os.getcwd(), 'result', cohort, 'subtype_dictionary.pickle'), 'rb') as f:
        subtype_dictionary = pickle.load(f)
    
    clinical_variables2 = np.intersect1d(np.array(CLINICAL_VARIABLES), np.array(list(subtype_dictionary.keys())))
    
    for i, c in enumerate(CLINICAL_VARIABLES):
        #print("----------\nclinical variable: "+c)
        if c in clinical_variables2:
            ax = fig.add_subplot(3, 3, i+1)

            #현재 clinical variable 값이 notnull인 샘플들만 남김.
            merged_tumor2 = merged_tumor[merged_tumor[c].notnull()].copy() 
            exclude_mask = np.array([True if str(merged_tumor2[c].values.flatten()[i]) in EXCLUDE_SUBTYPES else False for i in range(merged_tumor2.shape[0])])
            merged_tumor3 = merged_tumor2.iloc[~exclude_mask, :]

            # 각 clinical variable의 subtype의 개수 세기
            # 이 clinical variable은 임의로 정한 것이기 때문에 2_survival-analysis에서 구한 subtype_dictionary와 다를 수 있다.
            # 2번 코드에서 subtype_dictionary 만들 때는, 'not available' 등의 값은 빼고 subtype을 셌기 때문.


            subtypes = subtype_dictionary[c] 

            if len(subtypes)<len(COLORS):
                sample_num = 0
                for k in range(len(subtypes)):
                    current_subtype_samples = merged_tumor3[merged_tumor3[c]==subtypes[k]].bcr_patient_barcode.values.flatten()
                    current_subtype_samples_embedding = tumor_embedding_df2.loc[current_subtype_samples]
                    sample_num += current_subtype_samples_embedding.shape[0]

                    if k==0:
                        ax.scatter(current_subtype_samples_embedding.values[:,0], current_subtype_samples_embedding.values[:,1], 
                                   label = subtypes[k], color = COLORS[k])

                    else:
                        ax.scatter(current_subtype_samples_embedding.values[:,0], current_subtype_samples_embedding.values[:,1], 
                                   label = subtypes[k], color = COLORS[k])
                ax.set_title(c + ' (n = ' + str(sample_num) + ')' ) #이 clinical variable 값이 NaN이 아닌 샘플들의 전체 수
                ax.legend()
    fig.subplots_adjust(wspace = 0.3)
    fig.suptitle(cohort, fontsize = 15)
    fig.tight_layout()
    print("figure file: {}".format(os.path.join(savedir, figure_name)))
    plt.savefig(os.path.join(savedir, figure_name))     
    return tumor_embedding_df




def score_group_mds_score7(cohort, T, savedir, tumor_embedding_df, figure_name, fig_size, score):
    
    merged_tumor_fname = os.path.join(os.getcwd(), 'result', cohort, 'clinical_'+score+'_merged_TumorOnly.pickle')
    merged_tumor = pd.read_pickle(merged_tumor_fname)
    
    
    #clinical_score = pd.read_pickle(os.path.join(directory, CLINICAL_SCORE_FNAME))
    merged_tumor['score_group'] = ['Low' if float(merged_tumor.score7.values.flatten()[i]) < np.median(merged_tumor.score7.values.flatten()) 
                                   else 'High' for i in range(merged_tumor.shape[0])]
    #score_group_dictionary = {} 
    
    #distance_fname = '/data/project/jeewon/research/3D-ITH/binned_diff/score2_pc1/fire_compare_result/'+cohort+'/sample_sample_pc1_dist_'+str(chrom)+'.pickle'
    #distance = pd.read_pickle(distance_fname)
    #print("all-sample-pairwise distance: {}".format(distance.shape))
    
    #bcr_barcodes = [x[:12] for x in tumor_embedding_df.index.values]
    
    
    
    
    tumor_embedding_df.index = tumor_embedding_df['bcr_barcodes'].values.flatten()

    tumor_embedding_df2 = tumor_embedding_df[~tumor_embedding_df.index.duplicated(keep = 'last')]
    
    merged_tumor2 = merged_tumor[~merged_tumor.index.duplicated(keep = 'last')]
    
    tumor_embedding_df3 = pd.merge(tumor_embedding_df2, merged_tumor2, left_index = True, right_index = True)
    
    #merged_tumor.index = merged_tumor.full_barcode.values.flatten()
    
    fig = plt.figure(figsize = (fig_size, fig_size))
    ax = fig.add_subplot(111)
    sample_num = 0
    for k in range(tumor_embedding_df3['score_group'].nunique()): 
        #current_group_samples = merged_tumor[merged_tumor['score_group']==merged_tumor['score_group'].unique()[k]].bcr_patient_barcode.values.flatten()
        #current_group_embedding = tumor_embedding_df2.loc[current_group_samples]
        current_embedding = tumor_embedding_df3[tumor_embedding_df3['score_group'] == tumor_embedding_df3['score_group'].unique()[k]]
        ax.scatter(current_embedding['0'].values.flatten(), current_embedding['1'].values.flatten(), label = tumor_embedding_df3['score_group'].unique()[k], color = COLORS[k])
        sample_num += current_embedding.shape[0]
    ax.set_title('score_group MDS (' + cohort+ ' , n = '+str(sample_num)+')')
    ax.legend()

    fig.tight_layout()

    print("figure file: {}".format(os.path.join(savedir, figure_name)))
    plt.savefig(os.path.join(savedir, figure_name))      





def tumor_subtype_mds_score7(subtype_cohort, cohort, tumor_embedding_df, fig_size, savedir, fig_name, score):
    #change index of tumor_embedding_df to full barcode
    
    subtype_cohort.index = [x[:12] for x in subtype_cohort['pan.samplesID'].values.flatten()] 
    #pan.samplesID가 bcr_barcodes // full_barcode // full_barcode+['A' or 'B'] 등으로 들쭉날쭉해서 bcr_barcodes로 통일
    
    tumor_embedding_df.index = tumor_embedding_df['bcr_barcodes'].values.flatten()
    
    
    subtype_cohort2 = subtype_cohort[~subtype_cohort.index.duplicated(keep = 'last')]
    tumor_embedding_df2 = tumor_embedding_df[~tumor_embedding_df.index.duplicated(keep = 'last')]
    
    merged = pd.merge(tumor_embedding_df2, subtype_cohort2, left_index = True, right_index = True)
    #print('merged: ')
    #display(merged.head(3))
    #print(merged.shape)
    #print("\n")
    
    fig = plt.figure(figsize = (fig_size * 3, fig_size * 3)) 
    #subtype의 column 개수가 8개이기 때문에, 8개를 각각 다 plot할 수 있도록 9개의 subplot을 준비.
    
    targets = [] # subtypes to be plotted.
    for i in range(len(subtype_cohort2.columns)):
        if 'Subtype' in subtype_cohort2.columns[i]:
            targets.append(subtype_cohort2.columns[i])
    #print(targets)
    
    for i, t in enumerate(targets):
        merged2 = merged[merged[t].notnull()].copy()
        ax = fig.add_subplot(3, 3, i+1)
        
        if merged2[t].nunique() <= len(COLORS):
        
            for k in range(merged2[t].nunique()): # iterate for all categories belonging to current subtype
                #current_sample_full_barcodes = merged2[merged2[t] == merged2[t].unique()[k]].index.values
                #current_sample_embedding = merged2.loc[current_sample_full_barcodes].copy()
                current_embedding = merged2[merged2[t]==merged2[t].unique()[k]][['0', '1']]
                ax.scatter(current_embedding['0'].values.flatten(), current_embedding['1'].values.flatten(), label = merged2[t].unique()[k], color = COLORS[k])

            ax.legend()
            ax.set_title(t + ' (n = ' + str(merged2.shape[0]) + ')')
        
    fig.suptitle(cohort+' subtype MDS (n = ' + str(merged.shape[0]) + ', ' + score + ')', fontsize = 15)
    fig.tight_layout()
    
    print("figure file: {}".format(os.path.join(savedir, fig_name)))
    plt.savefig(os.path.join(savedir, fig_name)) 




if __name__=='__main__':
    
    args = parse_arguments()
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    print("cohort: {}".format(args.cohort))
        
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))
    print("len(tumor), len(normal), len(all): {}, {}, {}, respectively.".format(len(T), len(N), len(S)))
    
    cohort = args.cohort
    subtype_cohort = SUBTYPE[SUBTYPE['cancer.type'] == cohort.split('-')[-1]]
    
    # 함수 사용 예시
    #def tumor_normal_mds_score7(cohort, score, metric, T, N, S, savedir, figure_name, fig_size):
    #def clinical_variable_mds_score7(cohort, T, savedir, figure_name, fig_size):
    #def score_group_mds_score7(cohort, T, savedir, tumor_embedding_df, figure_name, fig_size, score):
    #def tumor_subtype_mds_score7(subtype_cohort, cohort, tumor_embedding_df, fig_size, savedir, fig_name, score):  

    
    all_sample_embedding_df = tumor_normal_mds_score7(args.cohort, args.score, args.distance, T, N, S, SAVEDIR, 'tumor_normal_MDS_'+args.score+'_'+args.distance+'.png', 5)
    plt.clf()
    
    tumor_embedding_df = clinical_variable_mds_score7(args.cohort, T, SAVEDIR, 'clinical_variable_MDS_' + args.score+'_'+args.distance+'.png', 5, args.distance, args.score) 
    plt.clf()
    
    score_group_mds_score7(args.cohort, T, SAVEDIR, tumor_embedding_df, 'score_group_MDS_' + args.score+'_'+args.distance+'.png', 5, args.score)
    plt.clf()
    
    tumor_subtype_mds_score7(subtype_cohort, args.cohort, tumor_embedding_df, 5, SAVEDIR, 'cancer_subtype_MDS_'+args.score+'_'+args.distance+'.png', args.score)
    plt.clf()
    
    # save MDS embedding for further analyses
    all_sample_embedding_fname = os.path.join(SAVEDIR, 'all_sample_MDS_embedding_'+args.score+'_'+args.distance+'.pickle')
    
    tumor_embedding_fname = os.path.join(SAVEDIR, 'tumor_sample_MDS_embedding_'+args.score+'_'+args.distance+'.pickle')
    
    print("all_sample_embedding_fname: {}".format(all_sample_embedding_fname))
    print("tumor_embedding_fname: {}".format(tumor_embedding_fname))
    
    all_sample_embedding_df.to_pickle(all_sample_embedding_fname)
    
    tumor_embedding_df.to_pickle(tumor_embedding_fname)

