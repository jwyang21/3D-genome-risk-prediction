# frequently used global variables and functions
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.integrate import simps
import random
random.seed(2022)
np.random.seed(2022)

CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
# TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
#NORMAL_BARCODES = ['10', '11', '12', '13','14', '15', '16', '17', '18', '19']
#CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
#CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #/{cohort}/{sample}.npz or /{cohort}/{sample}_inv_exp.npz  #'chr1'
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
FIRE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD'.split(' ')
NORMAL_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-GBM TCGA-READ TCGA-KIRC TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')
ALL_COHORT = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')
TCGA_SCORE3_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/'#{cohort}/score3_simple_avg.pickle
PCBC_SCORE3_FILE = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/PCBC/integrate-pcbc-abs-pc1.csv
SAVEDIR = os.path.join(os.getcwd(), 'result')
SAMPLE_NAME_FILE = '/data/project/jeewon/research/3D-ITH/data/samplename.npz'#item: {cohort}_samples
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]
print("SAVEDIR: {}".format(SAVEDIR))


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
    return parser.parse_args()


def preprocess_cpg_annot(CPG_ANNOT):
    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()]
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 
    return CPG_ANNOT2

def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort+'_samples']
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
