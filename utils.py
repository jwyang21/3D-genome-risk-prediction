# frequently used global variables and functions
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.integrate import simps
import random
random.seed(2022)
np.random.seed(2022)
import os
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import ttest_ind


CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
# TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
NORMAL_BARCODES = ['10', '11', '12', '13','14', '15', '16', '17', '18', '19']
## TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
#CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
CPG_ANNOT = pd.read_csv('/data/project/jeewon/research/3D-ITH/data/illumina/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)

CPG_METADATA = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', header=None, sep = '\t') # columns = ['chrom', 'start', 'end', 'cpg'] # use this!!


TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{cohort}/{sample}.npz or /{cohort}/{sample}_inv_exp.npz  #'chr1'
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/PCBC/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
TCGA_SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'#{cohort}/score3_simple_avg.pickle
PCBC_SCORE3_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv'
SAVEDIR = os.path.join(os.getcwd(), 'result')
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
ALL_COHORT_W_PCBC = ALL_COHORT.copy()
ALL_COHORT_W_PCBC.append('PCBC')
print("SAVEDIR: {}".format(SAVEDIR))
BINNED_DIFFMAT_BINS = '/data/project/jeewon/research/3D-ITH/pipelines/etc/binned-diffmat-bins' #pcbc_bins.npz #{TCGA_cohort}_diffmat_bins.npz
P_THRESHOLD = 5e-2
SCORE2_4_MINMAX_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/score2_score4_minmax.txt'
NORMAL_STEM_DISTANCE_MINMAX = pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/normal-stem-distance-min-max.csv', index_col = 0)


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True) #'TCGA-{}' or 'PCBC'
    #parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)
    parser.add_argument('-r_type', '--reference_type', help = 'reference type. PCBC or TCGA', required = True) #computing reference-v2
    parser.add_argument('-s_type', '--score_type', help = 'score type', required = True) #avg_pc1 #pc1_fluctuation #computing reference-v2
    parser.add_argument('-s2r', '--score2_reference', help = 'reference to be used in score2. fire or normal', default = 'NORMAL', required = False)
    parser.add_argument('-d', '--distance_metric', help = 'distance metric to be used in either score2 or score4. euclidean, jsd, or cosine-similarity', defulat = 'euclidean', required = False)
    parser.add_argument('-s4r', '--score4_reference', help = 'reference to be used in score4. all, sc, nonsc', default = 'all', required = False)
    parser.add_argument('-use_option', '--usage_option', help = 'use all samples or randomly picked samples', required = True) # all, part
    parser.add_argument('-n', '--normalize', help = 'whether you wnat to normalize score2 and score4 or not', default = 'N', required = False)#Y or N
    parser.add_argument('-m', '--minmax', help = 'use minmax scaling', default = 'N', required = False) #Y or N
    return parser.parse_args()


def preprocess_cpg_annot(CPG_ANNOT):
    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()]
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 
    return CPG_ANNOT2

def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort+'_samples']
    if cohort=='PCBC':
        return samples.tolist()
    else: #TCGA cohort
        T = []
        N = []
        S = [] #all samples
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
        S = T + N
        print("{}: tumor {}, normal {}, total {}".format(cohort, len(T), len(N), len(S)))

        return T, N, S


# default figure setting
mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family = 'FreeSans', size = 7)
plt.rc('figure', figsize = (1.5, 1.5))
    
def import_binned_diffmat(cohort, sample, chrom):
    # return a binned diffmat
    fname = os.path.join(BINNED_DIFFMAT_DIR, cohort, sample+'.npz')
    raw_diffmat = np.load(fname)[chrom]
    raw_mask = np.load(fname)['{}_mask'.format(chrom)]
    diffmat_masked = raw_diffmat[~raw_mask].T[~raw_mask].T

    return diffmat_masked, 1/np.exp(diffmat_masked)

def pc1(m):
    #print("pc1(m)")
    pca = PCA(n_components=3)
    pc = pca.fit_transform(m)
    #print(pc)
    pc1 = pc[:,0]
    #print(pc1)
    #print('-----')
    return pc1

def import_pc1(cohort, sample, chrom, flag, ref_type):
    # import pre-computed PC1 of sample-of-interest
    # ref_type: 'TCGA' or 'PCBC'. Type of reference you want to import. 
    # flag: 'raw' or 'inv'
    if ref_type=='TCGA':
        if flag=='raw':
            fname = os.path.join(TCGA_PC1_DIR, cohort, sample+'.npz')
        elif flag=='inv':
            fname = os.path.join(TCGA_PC1_DIR, cohort, sample+'_inv_exp.npz')
        else:
            pass

    elif ref_type == 'PCBC':
        if flag=='raw':
            fname = os.path.join(PCBC_PC1_DIR, sample+'.npz')
        elif flag=='inv':
            fname = os.path.join(PCBC_PC1_DIR, sample+'_inv_exp.npz')
        else:
            pass
        
    pc1 = np.load(fname)[chrom]
    return pc1

def integrate_abs_pc1(pc1_450k): 
    pc1_abs = np.array([abs(x) for x in pc1_450k])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area

def import_score(cohort, score, reference, distance):#import score2 and score4 #여기에 score3 추가. 
    if score=='score2':
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/'+score+'_'+reference+'_'+distance+'.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        mean_score = raw_score_df.mean(axis=1).values
        score_df = pd.DataFrame(mean_score, index = raw_score_df.index.values, columns = ['score2'])
         
    elif score=='score3':
        #SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/'#{cohort}/score3_simple_avg.pickle
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+'TCGA-LIHC/'+'TCGA-LIHC'+'_score3.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        
    elif score=='score4':
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/score4_'+reference_score4+'_'+distance_score4+'.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        score_df = pd.DataFrame(raw_score_df.simple_avg.values.flatten(), index = raw_score_df.index.values, columns = ['score4'])
    else: #add here if other types of scores are needed. (score1, score3, ...)
        pass
    return score_df

def compute_angle(x, y):#1사분면에서 (x,y)가 주어졌을때, (0,0)과 (x,y)를 잇는 직선이 x축과 이루는 각도 theta를 구하고 cosine(theta) 
    x = score2_norm.copy().flatten()
    y = score4_norm.copy().flatten()
    if len(x) != len(y):
        raise ValueError
    slope = np.array([float(y[i]/x[i]) for i in range(len(x))])
    radian_theta = np.array([math.atan(s) for s in slope])
    degree_theta = np.array([math.degrees(r) for r in radian_theta])
    cosine_radian = np.array([math.cos(r) for r in radian_theta])

    '''
    slope = y / x #slope of line = y / x
    radian_theta = math.atan(slope)
    degree_theta = math.degrees(radian_theta)
    cosine_radian = math.cos(radian_theta) 
    '''
    return radian_theta, degree_theta, cosine_radian

def compare_tumor_normal(score_df):
    #score_df should be a dataframe whose indices are sample names and have one column (score)
    ## here, name of the score column is 'cos_radian' (stem-closeness)
    tumor_mask = np.array([int(x[13:15])<=9 for x in score_df.index.values])
    tumor_score = score_df.iloc[tumor_mask,:].cos_radian.values.flatten()
    normal_score = score_df.iloc[~tumor_mask,:].cos_radian.values.flatten()
    print("---\nIndependent t-test between tumor and normal score")
    print(ttest_ind(tumor_score, normal_score))
    if ttest_ind(tumor_score, normal_score)[1]<=P_THRESHOLD:
        print("significant. p-value {} <= {}".format(ttest_ind(tumor_score, normal_score)[1], P_THRESHOLD))
    print("---\nTumor score mean and std")
    print("(mean, std) = ({}, {})".format(np.mean(tumor_score), np.std(tumor_score)))
    print("---\nNormal score mean and std")
    print("(mean, std) = ({}, {})".format(np.mean(normal_score), np.std(normal_score)))
    
def get_cpg_list(chr_list): #get list of opensea CpG probes
    total_list = np.array([])
    for chrom in chr_list:
        fname = '/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+chrom+'_opensea_CpG.pickle'
        cpgs = pd.read_pickle(fname).index.values
        total_list = np.union1d(total_list, cpgs)
        
    return total_list 

def get_sample_score(cohort, score_type, S, stemness_type, normalize, minmax):
    # 이 코호트의 모든 샘플들의 score 불러오기 (score2, score3, score4)
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
    
    elif score_type == 'score7': #version 1 (use all available samples to compute reference)
        df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result', cohort, 'score7.csv'), index_col=0) #index should be sample name.
        scores = df['cos_radian'].values.flatten()
        samples = df.index.values
        print("samples")
        print(samples)
    
    elif score_type == 'score7-without-minmax':
        df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result', cohort, 'score7-without-minmax.csv'), index_col = 0).loc[S]
        scores = df['cos_radian'].values.flatten()
        samples = df.index.values
        
    elif score_type == 'stem_closeness':
        if stemness_type == 'avg_pc1':
            if normalize == 'Y':
                df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', 
                                              cohort, 'stem-closeness_avg_pc1_part_normalized.csv'), index_col=0) #index should be sample name.
            else:
                if minmax=='Y':
                    pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', 
                                              cohort, 'stem-closeness_avg_pc1_part_minmax.csv'), index_col = 0)
                else:
                    df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', 
                                              cohort, 'stem-closeness_avg_pc1_part.csv'), index_col=0) #index should be sample name.
            
        
        else: #pc1_fluctuation
            if normalize == 'Y':
                df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', 
                                              cohort, 'stem-closeness_pc1_fluctuation_part_normalized.csv'), index_col=0)
            else:
                df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', 
                                              cohort, 'stem-closeness_pc1_fluctuation_part.csv'), index_col=0)
                
        #df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', cohort, 'stem-closeness_avg_pc1_part.csv'), index_col=0) #index should be sample name.
        # #index should be sample name.
        #df = pd.read_csv(os.path.join('/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/', cohort, 'stem-closeness_avg_pc1_all.csv'), index_col=0)
        scores = df['cos_radian'].values.flatten()
        samples = df.index.values     
    
    else:
        pass
    
    df = pd.DataFrame(scores, index = samples, columns = ['score'])
    
    #df = pd.DataFrame(zip(samples, scores), columns = ['sample', 'score'])
    return df.loc[S]

def get_cpg_list(chr_list):
    total_list = np.array([])
    for chrom in chr_list:
        fname = '/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+chrom+'_opensea_CpG.pickle'
        cpgs = pd.read_pickle(fname).index.values
        total_list = np.union1d(total_list, cpgs)
        
    return total_list 

def confirm_cpg_list_concordance(CHR_LIST):#binned diffmat 제작할때 쓴 opensea cpg probe들의 목록이 일치하는지 확인 
    #지금은 binned diffmat 버전을 CPG_METADATA 바탕으로 만들어진 버전으로 일괄적으로 사용하므로 필요x
    for x in CHR_LIST:
        df_tmp = df[df['chrom']==x]
        value1 = pd.read_pickle('/data/project/jeewon/research/3D-ITH/binned_diff/snake/'+x+'_opensea_CpG.pickle').shape[0]
        value2 = df_tmp.shape[0]
        assert value1==value2

'''
# starting main()
if __name__=='__main__':
    
    args = parse_arguments()
    
    SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    
    if not os.path.exists(os.path.join(os.getcwd(), 'result')):
        os.makedirs(os.path.join(os.getcwd(), 'result'))
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    
    print("cohort: {}".format(args.cohort))
    print("SAVEDIR: {}".format(SAVEDIR))
'''
