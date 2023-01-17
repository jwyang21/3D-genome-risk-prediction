python3 cox_v1.py --cohort TCGA-BLCA --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BLCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BLCA --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BLCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BLCA --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BLCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BLCA --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BLCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BLCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BRCA --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BRCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-7.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BRCA --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BRCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-7.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BRCA --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BRCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-7.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-BRCA --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-BRCA/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-7.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-BRCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-CHOL --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-CHOL/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_normalized_standardized_to-chrom-19.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-CHOL --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-CHOL/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_normalized_standardized_to-chrom-19.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-CHOL --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-CHOL/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_normalized_standardized_to-chrom-19.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-CHOL --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-CHOL/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_normalized_standardized_to-chrom-19.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-CHOL/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-COAD --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-COAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-20.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-COAD --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-COAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-20.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-COAD --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-COAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-20.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-COAD --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-COAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-20.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-COAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRC --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRC/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_to-chrom-3.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRC --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRC/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_to-chrom-3.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRC --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRC/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_to-chrom-3.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRC --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRC/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_to-chrom-3.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRP --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRP/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRP --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRP/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRP --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRP/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-KIRP --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-KIRP/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-KIRP/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LIHC --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LIHC/stem-closeness_cosine-sim_iebdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-9.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LIHC --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LIHC/stem-closeness_cosine-sim_iebdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-9.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LIHC --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LIHC/stem-closeness_cosine-sim_iebdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-9.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LIHC --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LIHC/stem-closeness_cosine-sim_iebdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-9.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LIHC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUAD --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUAD/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_standardized_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUAD --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUAD/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_standardized_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUAD --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUAD/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_standardized_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUAD --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUAD/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_normalized_standardized_to-chrom-1.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUSC --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUSC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-14.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUSC --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUSC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-14.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUSC --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUSC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-14.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-LUSC --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-LUSC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-14.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-LUSC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PAAD --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PAAD/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PAAD --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PAAD/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PAAD --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PAAD/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PAAD --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PAAD/stem-closeness_cosine-sim_bdm_pc1-avg_simple-avg_half_minmax_to-chrom-8.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PAAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PRAD --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PRAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PRAD --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PRAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PRAD --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PRAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-PRAD --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-PRAD/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-PRAD/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-THCA --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-THCA/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-THCA --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-THCA/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-THCA --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-THCA/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-THCA --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-THCA/stem-closeness_euclidean_bdm_pc1-avg_weighted-avg_half_minmax_standardized_to-chrom-2.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-THCA/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-UCEC --survival_type OS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/cox_v1 --result_file cox-regression-OS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-UCEC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-18.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-UCEC --survival_type DSS --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/cox_v1 --result_file cox-regression-DSS --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-UCEC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-18.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-UCEC --survival_type DFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/cox_v1 --result_file cox-regression-DFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-UCEC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-18.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

python3 cox_v1.py --cohort TCGA-UCEC --survival_type PFI --result_dir /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/cox_v1 --result_file cox-regression-PFI --score_file /data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/TCGA-UCEC/stem-closeness_euclidean_bdm_pc1-avg_simple-avg_half_minmax_standardized_to-chrom-18.csv --avg_beta_file /data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/TCGA-UCEC/opensea_tumors_avg_beta.csv --clinical_file /data/project/3dith/data/TCGA-CDR-SupplementalTableS1.csv

