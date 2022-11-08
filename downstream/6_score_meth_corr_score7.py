#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
from scipy.stats import pearsonr
import argparse


# ## 각 score의 filename
# - score2
#     - /data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/{cohort}/score2_{reference}_{metric}.pickle
#         - reference: 'fire', 'normal'
#         - metric: 'euclidean', 'cosine-similarity', 'jsd'
# - score3
#     - /data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/{cohort}/score3_simple_avg.pickle
# - score4
#     - /data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/{cohort}/score4_{ref_type}_{distance}.pickle
#         - distance: euclidean
#         - cohort: 33 TCGA cohorts
#         - ref_type: all, SC, nonSC
# - score5
#     - /data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/{cohort}/score2345.npz 
# - score7:
#    - /data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta/result/{cohort}/score2-score4-cosine.csv

# In[4]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses')


# In[5]:


mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[6]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)] #성염색체도 포함하면 chrY에서 mask시키면 아무 bin도 안 남아서 PCA할때 에러 남.
#BETA_FNAME = '/data/project/jeewon/research/3D-ITH/data/TCGA-LIHC_hg19/TCGA.LIHC.sampleMap%2FHumanMethylation450'
#BETA = pd.read_csv(BETA_FNAME, sep = '\t', index_col = 0)
SAVEDIR = os.path.join(os.getcwd(), 'result')
TUMOR_BARCODE = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
NORMAL_BARCODE = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result/'


# In[7]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--score', help = 'Type of score', type = str, required = True) #어떤 score로 구한 PC1으로 sample 간 거리를 계산하는지.
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    #parser.add_argument('-t', '--threshold', help = 'Threshold for score1', type = float, required = True) #for score 1
    #parser.add_argument('-r', '--reference', help = 'Reference for score2', type = str, required = True) #for score2
    #parser.add_argument('-d', '--distance', help = 'Distance metric for score2', type = str, required = True) #which distance metric you want to use
    return parser.parse_args()


# In[8]:


def get_sample_list(directory):
#def get_sample_list(directory):
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


# In[9]:


def get_cpg_list(chr_list):
    total_list = np.array([])
    for chrom in chr_list:
        fname = '/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+chrom+'_opensea_CpG.pickle'
        cpgs = pd.read_pickle(fname).index.values
        total_list = np.union1d(total_list, cpgs)
        
    return total_list 


# In[10]:


def get_beta_opensea(cohort, cpg_list, S):
    # 현재 cohort의 전체 beta 데이터 중 opensea cpg들의 beta 값만 반환.
    beta_fname = '/data/project/3dith/data/450k_xena/'+cohort+'.HumanMethylation450.tsv'
    beta = pd.read_csv(beta_fname, sep = '\t', index_col = 0)
    
    beta_opensea_df = beta.loc[cpg_list].dropna()
    print("num_opensea_cpg_after_dropna: {}".format(beta_opensea_df.shape[0]))
    samples_beta = beta_opensea_df.columns.values
    avg_opensea_beta = beta_opensea_df.mean().values #column mean
    
    avg_beta_opensea_df = pd.DataFrame(avg_opensea_beta, index = samples_beta, columns = ['avg_beta'])
    
    #avg_beta_opensea_df = pd.DataFrame(zip(samples_beta, avg_opensea_beta), columns = ['sample', 'avg_beta'])

    return beta_opensea_df[S], avg_beta_opensea_df.loc[S] 


# In[24]:


def get_sample_score(cohort, score_type, S):
    # 이 코호트의 score 불러오기 (score2, score3, score4)
    if score_type == 'score2':
        score2 = pd.read_pickle('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/score2_normal_euclidean.pickle')
        #display(score2.head(3))
        #print(score2.index.values[:5])
        #print(score2.mean(axis = 1).values[:5])
        #print(len(score2.index.values)==len(score2.mean(axis = 1)))
        samples = score2.index.values.flatten()
        scores = score2.mean(axis = 1).values.flatten()
    elif score_type == 'score3':
        # score3 예시
        score3 = pd.read_pickle('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/'+cohort+'/score3_simple_avg.pickle')
        #display(score3.head(3))
        #print(score3.index.values[:5])
        #print(score3.simple_avg.values.flatten()[:5])
        #print(len(score3.index.values)==len(score3.simple_avg.values.flatten()))
        samples = score3.index.values.flatten()
        scores = score3.simple_avg.values.flatten()
    elif score_type == 'score4':
        score4 = pd.read_pickle('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/score4_'+'all'+'_'+'euclidean'+'.pickle')
        #display(score4.head(3))
        #print(score4.index.values[:5])
        #print(score4.simple_avg.values.flatten()[:5])
        #print(len(score4.index.values)==len(score4.simple_avg.values.flatten()))
        samples = score4.index.values.flatten()
        scores = score4.simple_avg.values.flatten()
        
    elif score_type == 'score5':
        scores = np.load('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/score2345.npz', allow_pickle=True)['score5']
        samples = np.load('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/score2345.npz', allow_pickle=True)['rownames']
    
    elif score_type == 'score7':
        df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/score2-score4-theta/result/', cohort, 'score2-score4-cosine.csv'), index_col=0) #index should be sample name.
        scores = df['cos_radian'].values.flatten()
        samples = df.index.values
    else:
        pass
    df = pd.DataFrame(scores, index = samples, columns = ['score'])
    
    #df = pd.DataFrame(zip(samples, scores), columns = ['sample', 'score'])
    return df.loc[S]


# In[13]:


def scatter_plot(fig_size, score, beta_score_df, cohort_dir, cohort):
    fig = plt.figure(figsize = (fig_size, fig_size))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Average opensea CpG methylation level', fontsize = 13)
    ax.set_ylabel(score, fontsize = 13)
    ax.scatter(beta_score_df.avg_beta.values.flatten(), beta_score_df.score.values.flatten())
    figname = 'scatter-avg_opensea_beta-'+score+'.png'
    print('figure name: {}'.format(os.path.join(cohort_dir, figname)))
    ax.set_title(cohort, fontsize = 15)
    plt.savefig(os.path.join(cohort_dir, figname))
    #plt.clf()


# In[14]:



if __name__=='__main__':
    
    #for cohort in COHORTS:
    args = parse_arguments()
    cohort_dir = os.path.join(SAVEDIR, args.cohort)  #현재 cohort의 결과가 이 디렉토리에 저장돼야 함. 
    print("cohort: {}".format(args.cohort))
    T, N, S = get_sample_list(os.path.join(BINNED_DIFFMAT_DIR, args.cohort))

    
    #/data/project/3dith/data/450k_xena/TCGA-[XXXX].HumanMethylation450.tsv
    
    # 먼저, 22개 염색체들의 opensea CpG probe ID들의 합집합을 구하기 
    cpg_list = get_cpg_list(CHR_LIST) 
    print("num_opensea_cpg: {}".format(len(cpg_list)))
    
    # 현재 cohort의 beta 값 불러들이기
    beta_opensea_df, avg_beta_opensea_df = get_beta_opensea(args.cohort, cpg_list, S) #beta_opensea_df: column이 sample, index가 CpG -> 여기서 column mean 
                                                                            #-> avg_beta_opensea_df: index가 sample, column이 average opensea beta value.
    
    # beta_opensea의 각 column별 mean (각 sample 별 average opensea beta value)를 구하기

    '''
    # 나중에 쓸 수도 있으니까 저장
    npz_fname = os.path.join(cohort_dir, 'avg_beta_opensea')
    np.savez(npz_fname, sample = avg_beta_opensea_df.index.values.flatten(), value = avg_beta_opensea_df.avg_beta.values.flatten())
    print('npz filename: ', end = '')
    print(npz_fname+'.npz')
    '''    
    
    # 현재 cohort의 각 sample의 score 값 import
    score_df = get_sample_score(args.cohort, args.score, S) #score_df: index가 sample, column이 average opensea beta value
    
    beta_score_df = avg_beta_opensea_df.merge(score_df, left_index = True, right_index = True)
    
    # score과 avg beta value 간 PCC
    print('PCC between {} and {}'.format('avg_opensea_beta', args.score))
    pcc, pvalue = pearsonr(beta_score_df.avg_beta.values.flatten(), beta_score_df.score.values.flatten())
    print("PCC: {}, pvalue: {}".format(pcc, pvalue))
    
    # Scatter plot (X: average opensea beta value, Y: score)
    #def scatter_plot(fig_size, score, beta_score_df, cohort_dir, figname):
    scatter_plot(5, args.score, beta_score_df, cohort_dir, args.cohort)


# In[ ]:




