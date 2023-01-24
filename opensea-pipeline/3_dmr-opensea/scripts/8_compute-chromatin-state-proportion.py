import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--cpg_type', help = 'island, opensea, shelf_shore', type = str, required = True)
    args.add_argument('--cohort2eid', help = 'cohort2eid fname', type = str, default = '/data/project/3dith/data/etc/cohort2eid.txt', required = False)
    args.add_argument('--chromatin_states', help = 'chromatin state fname', type = str, default = '/data/project/3dith/data/chromatin_states.npy', required = False)
    args.add_argument('--input_fname', help = 'input filename', type = str, default = 'DMR_EPI_features_threshold_mean_std_len.npz', required = True)
    return args.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    working_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/'
    save_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/result/'
    cohort2eid = pd.read_csv(args.cohort2eid, sep = '\t', header = None, names = ['cohort', 'eid'])
    chromatin_states = np.load(args.chromatin_states, allow_pickle = True)
    save_fname_p = f'EPI-category-proportion-stacked-bar-chart-{args.cpg_type}-TN.csv'
    save_fname_l = f'EPI-category-len-stacked-bar-chart-{args.cpg_type}-TN.csv'
    
    all_df = pd.DataFrame(np.zeros((len(cohort2eid.cohort.values), len(chromatin_states)), dtype = float), index = cohort2eid.cohort.values, columns = chromatin_states)
    
    for cohort in cohort2eid.cohort.values:
        print(cohort)
        input_dir = f'/data/project/3dith/pipelines/{args.cpg_type}-pipeline/3_dmr-{args.cpg_type}/result/{cohort}'
        data = np.load(os.path.join(input_dir, args.input_fname))
        for k in list(data.keys()):
            for s in chromatin_states:
                if s in k:
                    all_df.loc[cohort][s] += data[k]
                    
    all_df_proportion = pd.DataFrame(np.zeros((len(cohort2eid.cohort.values), len(chromatin_states)), dtype = float), index = cohort2eid.cohort.values, columns = chromatin_states)
    for cohort in cohort2eid.cohort.values:
        for s in chromatin_states:
            all_df_proportion.loc[cohort][s] = all_df.loc[cohort][s]/all_df.loc[cohort].sum()
            
    all_df.to_csv(os.path.join(save_dir, save_fname_l))
    print(os.path.join(save_dir, save_fname_l))
    
    all_df_proportion.to_csv(os.path.join(save_dir, save_fname_p))
    print(os.path.join(save_dir, save_fname_p))    
