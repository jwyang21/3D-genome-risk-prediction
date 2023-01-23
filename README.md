# Stem closeness
## 1. overview
### 1-1. qualitative overview
![230117_github_readme_-overview](https://user-images.githubusercontent.com/86412887/212884165-b1908130-92cb-4623-8d48-ebbde1cda9ce.png)
### 1-2. quantitative overview
![30119_figurex_stem-closeness-computation-explanation-BDM-300dpi](https://user-images.githubusercontent.com/86412887/213362375-203bb05a-5253-49bb-b39e-53bd4b8c645f.png)
## 2. Installation       
There are two environments needed: one for processing Hi-C data, and the other for computing stem closeness and running downstream analyses.     
### 2-1. Installing conda environment processing Hi-C data
```shell
conda install mamba -n base -c conda-forge
mamba env create --file hic-processing/environment.yaml
```
### 2-2. Installing conda environment for running stem closeness-related analyses
```shell
conda install mamba -n base -c conda-forge
mamba env create --file stem-closeness.yaml
```
## 3. Process

### 3-1. Process Hi-C data 
```shell
cd hic-processing
snakemake --cores 100 --resources network=1 --restart-time 3
python ab2corrmat.py
cd ../
```

### 3-2. Conduct stem closeness-related analyses
```shell
conda activate stem-closeness
```
#### 3-2-1. Preprocessing data
- Construct metadata and bedgraph file of open sea CpG probes in 450K DNA methylation data
  - Download manifest file for Infinium HumanMethylation450 v1.2 BeadChip, provided by Illumina ([1]) 
```shell
cd data
wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv
cd ../
```
  - Make metadata and bedgraph files
```shell
cd utils
python3 make-450k-probe-metadata.py --cpg_type opensea
python3 make-450k-probe-bedgraph.py --cpg_type opensea
cd ../
```
- Construct list of sample names per cohort
```shell
cd utils/scripts
python3 0_get_sample_name_list.py
cd ../../
```
#### 3-2-2. Make binned difference matrices (BDMs)    
##### 3-2-2-1. Make BDMs for TCGA samples, using open sea CpG probes
```shell
cd binned-difference-matrix-v2-opensea
snakemake -j 10
cd ../
```
##### 3-2-2-2. Make BDMs for PCBC stem cell samples, using open sea CpG probes
```shell
cd binned-difference-matrix-pcbc
snakemake -j 10
cd ../
```
##### 3-2-2-3. Make list of genomic bins used for constructing BDM.
```shell
cd utils/scripts
bash 1_find-bdm-bins.sh
cd ../../
```

#### 3-2-3. Run opensea pipeline
##### 3-2-3-1. Compute stem closeness
```shell
cd opensea-pipeline/1_compute-score-opensea/scripts
bash 1_all-samples-pc1.sh > ../log/1_all-samples-pc1.log
bash 2_compute-reference.sh > ../log/2_compute-reference.log
bash 3_compute-distance.sh > ../log/3_compute-distance.log
bash 4_compute-min-max-distance.sh > ../log/4_compute-min-max-distance.log
bash 5_compute-sc-cosine_and_euclidean-minmax.sh > ../log/5_compute-sc-cosine_and_euclidean-minmax.log
bash 5_compute-sc-euclidean-normalized.sh > ../log/5_compute-sc-euclidean-normalized.log
bash 5_compute-sc-euclidean-raw.sh > ../log/5_compute-sc-euclidean-raw.log
cd ../../../
```

##### 3-2-3-2. Downstream analyses 
- Make representative PC1 vectors of tumor, normal, and stem cells, respectively.
```shell
cd utils/scripts
python3 2_make-TCGA-repr-bdm.py
python3 3_make-PAAD-repr-vector-bdm.py
python3 4_make-SC-repr-bdm.py
python3 5_bdm-iebdm-pc1-pcc.py
cd ../../
```
- Download TCGA clinical data.
```shell
cd data
wget https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
cd ../
```
- Check whether DNA methylation-derived PC1s can reproduce those of HiC-PC1s
```shell
cd opensea-pipeline/2_downstream-opensea/scripts
bash compare-Tumor-hic-450k-pc1-repr.sh
bash compare-Tumor-hic-450k-pc1-individual.sh
bash compare-Stem-hic-450k-pc1-repr.sh
bash compare-Stem-hic-450k-pc1-individual.sh
bash compare-Normal-FIRE-hic-450k-pc1-repr.sh
bash compare-Normal-FIRE-hic-450k-pc1-individual.sh
bash compare-Normal-3DIV-450k-pc1-repr.sh
bash compare-Normal-3DIV-450k-pc1-individual.sh
```
- Check whether binned PC1s are tissue type-specific
```shell
python3 check-tissue-specificity.py
```
- Find out correlation between stem closeness and avg beta, the average DNA methylation level of open sea CpG probes.
```shell
bash pcc-avg_beta-stem_closeness.sh > ../log/pcc-avg_beta-stem_closeness-ALL.log
bash parse-avg_meth-score-corr-log.sh
cd ../../../
```
- Assign score group and beta group
```shell
cd utils/scripts
python3 6_assign-score_group.py
python3 7_assign-beta_group.py
cd ../../
```
- Conduct log-rank test (1): use score group as predictor
```shell
cd opensea-pipeline/2_downstream-opensea/scripts
bash sc-cosine_and_euclidean-minmax.sh > ../log/sc-cosine_and_euclidean-minmax.log
bash sc-euclidean-normalized.sh > ../log/sc-euclidean-normalized.log
bash sc-euclidean-raw.sh > ../log/sc-euclidean-raw.log
```
- Conduct log-rank test (2): use beta group as predictor
```shell
bash survival-analysis-avg_beta.sh > ../log/survival-analysis-avg_beta.log
```
- Write commands to run Cox regression
```shell
python3 write-cox-commands.py --cpg_type opensea --cox_version 5 --command_fname cox-commands --score_fname /data/project/3dith/data/cohort-1-best-score-km.csv
```
- Run Cox regression
```shell
bash cox_commands_v5.sh > ../log/cox_commands_v5.log
```
- Parse Cox regression results
```shell
bash parse_cox_results_v5.sh
```
- Plot Cox regression results (hazard ratios)
```shell
python3 plot-HR.py
cd ../../../
```
##### 3-2-3-3. DMR analysis
```shell
cd opensea-pipeline/3_dmr-opensea/scripts
bash 0_compute_binned_avg_opensea_beta.sh > ../log/0_compute_binned_avg_opensea_beta.log
bash 1_compute-binned-opensea-beta-TN-diff-binsize-1e6.sh > ../log/1_compute-binned-opensea-beta-TN-diff-binsize-1e6.log
bash 2_find_DMR.sh > ../log/2_find_DMR.log
bash 3_find_DMR_features_mean_std_epigenome.sh > ../log/3_find_DMR_features_mean_std_epigenome.log
bash 3_find_DMR_features_mean_std_gene_reg.sh > ../log/3_find_DMR_features_mean_std_gene_reg.log
bash 4_collate_DMR_features.sh > ../log/4_collate_DMR_features.log
bash 5_make-cohort-threshold-feature-list-threshold-mean_std.sh > ../log/5_make-cohort-threshold-feature-list-threshold-mean_std.log
bash 6_write_np2txt.sh > ../log/6_write_np2txt.log
bash 7_gseapy-gene-functional-annot.sh > ../log/7_gseapy-gene-functional-annot.log
bash 8_compute-chromatin-state-proportion.sh > ../log/8_compute-chromatin-state-proportion.log
cd ../../../
```
#### 3-2-4. Run island pipeline
##### 3-2-4-1. DMR analysis
```shell
cd island-pipeline/3_dmr-island/scripts/
bash 0_compute_binned_avg_island_beta.sh > ../log/0_compute_binned_avg_island_beta.log
bash 1_compute-binned-island-beta-TN-diff-binsize-1e6.sh > ../log/1_compute-binned-island-beta-TN-diff-binsize-1e6.log
bash 2_find_DMR.sh > ../log/2_find_DMR.log
bash 3_find_DMR_features_mean_std_epigenome.sh > ../log/3_find_DMR_features_mean_std_epigenome.log
bash 3_find_DMR_features_mean_std_gene_reg.sh > ../log/3_find_DMR_features_mean_std_gene_reg.log
bash 4_collate_DMR_features.sh > ../log/4_collate_DMR_features.log
bash 5_make-cohort-threshold-feature-list-threshold-mean_std.sh > ../log/5_make-cohort-threshold-feature-list-threshold-mean_std.log
bash 6_write_np2txt.sh > ../log/6_write_np2txt.log
bash 7_gseapy-gene-functional-annot.sh > ../log/7_gseapy-gene-functional-annot.log
cd ../../../
```
#### 3-2-5. Run shelf\_shore pipeline
##### 3-2-5-1. DMR analysis
```shell
cd shelf_shore-pipeline/3_dmr-shelf_shore/scripts
bash 0_compute_binned_avg_shelf_shore_beta.sh > ../log/0_compute_binned_avg_shelf_shore_beta.log
bash 1_compute-binned-shelf_shore-beta-TN-diff-binsize-1e6.sh > ../log/1_compute-binned-shelf_shore-beta-TN-diff-binsize-1e6.log
bash 2_find_DMR.sh > ../log/2_find_DMR.log
bash 3_find_DMR_features_mean_std_epigenome.sh > ../log/3_find_DMR_features_mean_std_epigenome.log
bash 3_find_DMR_features_mean_std_gene_reg.sh > ../log/3_find_DMR_features_mean_std_gene_reg.log
bash 4_collate_DMR_features.sh > ../log/4_collate_DMR_features.log
bash 5_make-cohort-threshold-feature-list-threshold-mean_std.sh > ../log/5_make-cohort-threshold-feature-list-threshold-mean_std.log
bash 6_write_np2txt.sh > ../log/6_write_np2txt.log
bash 7_gseapy-gene-functional-annot.sh > ../log/7_gseapy-gene-functional-annot.log
```


## Reference
[1] Bibikova, Marina, et al. "High density DNA methylation array with single CpG site resolution." Genomics 98.4 (2011): 288-295.     
