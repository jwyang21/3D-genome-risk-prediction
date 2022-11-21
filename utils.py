# frequently used global variables and functions
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.integrate import simps

CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask'
TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
NORMAL_BARCODES = ['11', '12', '13','14', '15', '16', '17', '18', '19']
## TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
TCGA_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #{sample}.npz or {sample}_inv_exp.npz  #'chr1'
PCBC_PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/pcbc/'#{sample}.npz or {sample}_inv_exp.npz #'chr1'
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv #columns: ['chrom','cpg']
SAVEDIR = os.path.join(os.getcwd(), 'result')
print("SAVEDIR: {}".format(SAVEDIR))

def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    #parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True)
    parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)

    return parser.parse_args()

def preprocess_cpg_annot(CPG_ANNOT):
    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()]
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 
    return CPG_ANNOT2

def get_sample_list(cohort):
    # sample list of input TCGA cohort
    cohort_binned_diffmat_dir = os.path.join(BINNED_DIFFMAT_DIR, cohort)
    T = []
    N = []
    S = [] #all samples
    for l in os.listdir(cohort_binned_diffmat_dir):
        if l.startswith('TCGA') and l.endswith('.npz'):
            if l[13:15] in TUMOR_BARCODES:
                T.append(l.split('.')[0].strip())
            elif l[13:15] in NORMAL_BARCODES:
                N.append(l.split('.')[0].strip())
            else:
                pass
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

def import_pc1(cohort, sample, chrom, flag):
    # import pre-computed PC1 of sample-of-interest
    if flag=='raw':
        fname = os.path.join(PC1_DIR, cohort, sample+'.npz')
    elif flag=='inv':
        fname = os.path.join(PC1_DIR, cohort, sample+'_inv_exp.npz')
    pc1 = np.load(fname)[chrom]

    return pc1

def integrate_abs_pc1(pc1_450k): 
    pc1_abs = np.array([abs(x) for x in pc1_450k])
    area = simps(pc1_abs, np.arange(len(pc1_abs)))
    return area

def import_score(cohort, score, reference, distance):#import score2 and score4
    if score=='score2':
        pickle_fname = '/data/project/jeewon/research/3D-ITH/pipelines/compute-score/result/'+cohort+'/'+score+'_'+reference+'_'+distance+'.pickle'
        raw_score_df = pd.read_pickle(pickle_fname)
        mean_score = raw_score_df.mean(axis=1).values
        score_df = pd.DataFrame(mean_score, index = raw_score_df.index.values, columns = ['score2'])
         
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
