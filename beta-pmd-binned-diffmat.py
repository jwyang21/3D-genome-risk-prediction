#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# To-do: compute binned difference matrix using CpG probes which (1) belong to PMD and (2) has notnull beta value


import pandas as pd
import numpy as np
import os
import argparse


# In[ ]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/beta-pmd-binned-diffmat/')


# ## To-do: compute binned difference matrix using CpG probes in PMD only.

# In[ ]:


CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
DATA_DIR = '/data/project/jeewon/research/3D-ITH/data'
BINNED_DIFFMAT_DIR = '/data/project/3dith/pipelines/binned-difference-matrix-v2/result' #'chr1', 'chr1_mask' #used when extracting sample list of cohort.
TUMOR_BARCODES = ['01', '02', '03','04', '05', '06', '07', '08', '09']
NORMAL_BARCODES = ['11', '12', '13','14', '15', '16', '17', '18', '19']
## TCGA barcode: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
#PC1_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/' #'chr1'
PMD_CPG_DIR = '/data/project/jeewon/research/3D-ITH/pipelines/find-pmd/result' #./{cohort}/pmd_cpg.csv
CHR_LENGTH = pd.read_csv('/data/project/jeewon/research/reference/GRCh37_hg19_chr_length.csv')[['Chromosome', 'Total_length']]
CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
METH_DIR = '/data/project/3dith/data/450k_xena/'#TCGA-[XXXX].HumanMethylation450.tsv'
SAVEDIR = os.path.join(os.getcwd(), 'result')
#print("SAVEDIR: {}".format(SAVEDIR))


# In[ ]:


def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', help='Beta bedgraph file.', required=True)
    #parser.add_argument('-s', '--chrom-size', help='Chromosome size table.', required=True)
    parser.add_argument('-b', '--binsize', type=int, default=int(1e6))
    #parser.add_argument('-c', '--n-min-cpgs', type=int, default=1)
    #parser.add_argument('-o', '--output', help='Output.', required=True)
    parser.add_argument('-ch', '--cohort', help = 'TCGA-cohort', required = True)
    parser.add_argument('-cr', '--chrom', help = 'chromosome', required = True)

    return parser.parse_args()


# In[ ]:


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


# In[ ]:


def preprocess_cpg_annot(CPG_ANNOT):
    CPG_ANNOT2 = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()] #drop any row whose chrom value is null.
    chrom_string = [str(x).split('.')[0] for x in CPG_ANNOT2.CHR.values.flatten()]
    CPG_ANNOT2['CHR'] = chrom_string #convert mixture of int and string types into string type only. 
    return CPG_ANNOT2


# In[ ]:


def get_pmd_cpg(cohort, chrom):
    # return list of CpG probes which belong to PMD in current chromosome. 
    pmd_cpg = pd.read_csv(os.path.join(PMD_CPG_DIR, cohort, 'pmd_cpg.csv'))#chrom, cpg
    pmd_cpg_list = pmd_cpg[pmd_cpg['chrom']==chrom].cpg.values.flatten()
    
    npy_fname = os.path.join(cohort_dir, chrom+'_pmd_cpg_list')
    np.save(npy_fname, np.array(pmd_cpg_list))
    print('pmd_cpg_list npy file: {}'.format(npy_fname+'.npy'))
    
    return pmd_cpg_list


# In[ ]:


def save_dictionary_npz(d, fname):
    #d: dictionary
    #fname: resulting filename
    #i: items to be included in npz file.
    key_list = list(d.keys())
    value_list = []
    for l in range(len(key_list)):
        current_key = key_list[l]
        current_value = d[current_key]
        value_list.append(current_value)
    if len(key_list)!=len(value_list):
        raise ValueError
    np.savez(fname, key=key_list, value=value_list)
    print("dictionary_npz filename: {}".format(fname+'.npz'))
    print("items included in the above npz file: {}".format(['key','value']))


# In[ ]:


def get_pmd_cpg_annot(cohort, chrom, cpg_annot, pmd_cpg_list, chrom_len, binsize, cohort_dir):
    
    # from whole cpg annotation data, extract cpg probes belonging to PMD of current chromsome
    pmd_cpg_annot_df = cpg_annot.loc[pmd_cpg_list]
    
    bin_pmd_cpg = {} #현재 chromosome의 각 bin당 pmd_cpg들의 목록 (key: bin, value: pmd-cpg)
    
    n_bins = (int(chrom_len)//int(binsize))+1
    print("Total {} bins in {}".format(n_bins, chrom))
    
    
    for i in range(n_bins):#for real
    #for i in range(3):#for test
    #if i % 500 == 0:
        #print("----\n")
        #print("processing {}-th bin.".format(i))
        #bin name format -> "chr1:1-1000000"
        bin_start = int(i * binsize)+1
        bin_end = int((i+1) * binsize)
        current_bin_name = chrom+':'+str(bin_start)+'-'+str(bin_end)
        for j in range(len(pmd_cpg_list)):#for real
        #for j in range(10):#for test
            if i % 30 == 0 and j % 100000 == 0:
                print("processing {}-th bin and {}-th CpG probe.".format(i, j))
            current_cpg = pmd_cpg_list[j]
            current_cpg_mapinfo = int(CPG_ANNOT.loc[current_cpg].MAPINFO)
            if current_cpg_mapinfo >= bin_start and current_cpg_mapinfo <= bin_end:
                if current_bin_name not in list(bin_pmd_cpg.keys()):
                    bin_pmd_cpg[current_bin_name] = []
                bin_pmd_cpg[current_bin_name].append(current_cpg)
        
    # save bin_pmd_cpg and pmd_cpg_annot_df
    save_dictionary_npz(bin_pmd_cpg, os.path.join(cohort_dir, chrom+'_bin-pmd-cpg')) #bin_pmd_cpg
    
    pmd_cpg_annot_df.to_csv(os.path.join(cohort_dir, chrom+'_pmd-cpg-annot-df.csv'), index=False) #pmd_cpg_annot_df
    print(os.path.join(cohort_dir, chrom+'_pmd-cpg-annot-df.csv'))
    return bin_pmd_cpg, pmd_cpg_annot_df


# In[ ]:


def compute_binned_diffmat(bin_pmd_cpg, pmd_cpg_annot_df, beta, sample, min_num_cpg):
    
    beta2 = beta[sample].dropna()
    cpg_after_dropna = beta2.index.values
    
    bin_pmd_cpg2 = {}
    for l in range(len(list(bin_pmd_cpg.keys()))):
        current_key = list(bin_pmd_cpg.keys())[l]
        current_value = bin_pmd_cpg[current_key]
        current_value_dropna = np.intersect1d(beta2.index.values, current_value)
        if len(current_value_dropna) >= min_num_cpg:
            bin_pmd_cpg2[current_key] = current_value_dropna
            
    n_pmd_cpg_bins = len(list(bin_pmd_cpg2.keys()))
    binned_diffmat = pd.DataFrame(np.zeros((n_pmd_cpg_bins, n_pmd_cpg_bins), dtype = float), index = list(bin_pmd_cpg2.keys()), columns = list(bin_pmd_cpg2.keys()))
    
    for i in range(n_pmd_cpg_bins):
        current_bin1 = list(bin_pmd_cpg2.keys())[i]
        current_bin_cpg1 = bin_pmd_cpg2[current_bin1]
        beta_current_cpg1 = beta2.loc[current_bin_cpg1]
        for j in range(i, n_pmd_cpg_bins):#binned-diffmat is symmetric matrix.
            current_bin2 = list(bin_pmd_cpg2.keys())[j]
            current_bin_cpg2 = bin_pmd_cpg2[current_bin2]
            beta_current_cpg2 = beta2.loc[current_bin_cpg2]
            
            median1 = np.median(beta_current_cpg1.values.flatten())
            median2 = np.median(beta_current_cpg2.values.flatten())
            
            binned_diffmat[current_bin1][current_bin2] = abs(median1-median2)
            binned_diffmat[current_bin2][current_bin1] = abs(median1-median2)
        
    
    return binned_diffmat


# In[ ]:


if __name__ == '__main__':
    args = parse_arguments()
    
    T, N, S = get_sample_list(args.cohort)
    
    cohort_dir = os.path.join(SAVEDIR, args.cohort)
    
    if not os.path.exists(SAVEDIR):
        os.makedirs(SAVEDIR)
    if not os.path.exists(cohort_dir):
        os.makedirs(cohort_dir)
    
    print("SAVEDIR: {}.format(SAVEDIR)")
    print("cohort: {}".format(args.cohort))
    print("cohort_dir: {}".format(cohort_dir))
    
    print("Import beta value.")
    beta = pd.read_csv(os.path.join(METH_DIR, args.cohort+'.HumanMethylation450.tsv'), sep = '\t', index_col=0) #row: cpg, column: TCGA sample
    
    cpg_annot2 = preprocess_cpg_annot(CPG_ANNOT)
    
    chrom_len = int(CHR_LENGTH[CHR_LENGTH['Chromosome']==args.chrom].Total_length.values[0])
    print("length of {}: {}".format(args.chrom, chrom_len))
    
    if os.path.exists(os.path.join(cohort_dir, args.chrom+'_pmd_cpg_list.npy')):
        pmd_cpg_list = np.load(os.path.join(cohort_dir, args.chrom+'_pmd_cpg_list.npy'), allow_pickle=True)
    else:
        pmd_cpg_list = get_pmd_cpg(args.cohort, args.chrom) #pmd-cpg와 beta data의 cpg들 중 겹치는 것을 찾기
        
    if os.path.exists(os.path.join(cohort_dir, args.chrom+'bin-pmd-cpg.npz')) and os.path.exists(os.path.join(cohort_dir, args.chrom+'_pmd-cpg-annot-df.csv')):
        pmd_cpg_annot_df = pd.read_csv(os.path.join(cohort_dir, args.chrom+'_pmd-cpg-annot-df.csv'))
        
        bin_pmd_cpg_key = np.load(os.path.join(cohort_dir, args.chrom+'bin-pmd-cpg.npz'),allow_pickle=True)['key']
        bim_pmd_cpg_value = np.load(os.path.join(cohort_dir, args.chrom+'bin-pmd-cpg.npz'),allow_pickle=True)['value']
        if len(bin_pmd_cpg_key)!=len(bin_pmd_cpg_value):
            raise ValueError
        else:
            bin_pmd_cpg = {}
            for i in range(len(bin_pmd_cpg_key)):
                bin_pmd_cpg[bin_pmd_cpg_key[i]] = bin_pmd_cpg_value[i]
            
    else: ##pmd-cpg-probe-list들이 binning되어 있는, key가 bin_name이고 value가 그 bin에 속한 pmd-cpg들인 dictionary 파일이 존재하지 않으면
        bin_pmd_cpg, pmd_cpg_annot_df = get_pmd_cpg_annot(args.cohort, args.chrom, cpg_annot2, pmd_cpg_list, chrom_len, int(args.binsize), cohort_dir)
        ##위에서 찾은 cpg들의 mapinfo 데이터를 얻고 binning한 후, resulting df를 csv로 저장.        
    
    for s in S: #for all samples
        if S.index(s) % 50 ==0:
            print('----\nprocessing {}-th sample: {}'.format(S.index(s), s))
        binned_diffmat = compute_binned_diffmat(bin_pmd_cpg, pmd_cpg_annot_df, beta, s, 1)
        ## binned difference matrix of current sample, current chromosome.
        npz_fname = os.path.join(cohort_dir, s+'_'+args.chrom) 
        ## for each sample, 22 binned diffmats corresponding to 22 chromosomes are made.
        np.savez(npz_fname, value = binned_diffmat.values, colnames = binned_diffmat.columns.values, rownames = binned_diffmat.index.values)
        print("npz fname: {}".format(npz_fname+'.npz'))  


# In[ ]:
