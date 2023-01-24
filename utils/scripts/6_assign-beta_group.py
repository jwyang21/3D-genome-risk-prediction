import pandas as pd
import numpy as np
import os

score_files = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)
num_groups = 2
group_labels = ['Low', 'High']
all_cpg_types = ['opensea', 'island', 'shelf_shore']
score_files.index.values

for cpg_type in all_cpg_types:
    for i in range(score_files.shape[0]):
        cohort = score_files.index.values[i]
        print("===\ncohort: {}".format(cohort))
        avg_beta_fname = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/{cohort}/{cpg_type}_tumors_avg_beta.csv'
        avg_beta_df = pd.read_csv(avg_beta_fname, index_col = 0)
        print("avg_beta_df: ", avg_beta_df.shape)
        tumor_mask = np.array([int(x[13:15]) <= 9 for x in avg_beta_df.index.values]) 
        avg_beta_df = avg_beta_df.iloc[tumor_mask, :]
        print("avg_beta_df after extracting tumors only: ", avg_beta_df.shape)
        avg_beta_df['group'] = pd.qcut(avg_beta_df.avg_beta, q = num_groups, labels = list(range(num_groups)))
        group_name = [group_labels[x] for x in avg_beta_df.group.values]
        avg_beta_df['group_name'] = group_name        
        result_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/{cohort}'
        avg_beta_df_fname = 'avg-beta-group.csv'
        print(os.path.join(result_dir, avg_beta_df_fname))
        avg_beta_df.to_csv(os.path.join(result_dir, avg_beta_df_fname))
        print("avg_beta_df shape: ", avg_beta_df.shape)
