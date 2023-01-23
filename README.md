# Stem closeness
## 1. Introduction and overview
## 1-1. Introduction    
Stem closeness is a novel metric to measure stem-likeness of a single sample. The overall pipelinen for computing stem closeness consists of 3 main steps: (1) Inferring 3D genome structure of a single sample, using DNA methylation data only (note that inferred 3D genome state is representated as a vector) -> (2) Computing distances between inferred 3D genome states -> (3) Use combinations of distance to quantify stem closeness. The figurew below describe more detailed overview of this pipeline.
### 1-2. Overview
#### 1-2-1. Qualitative overview
![230117_github_readme_-overview](https://user-images.githubusercontent.com/86412887/212884165-b1908130-92cb-4623-8d48-ebbde1cda9ce.png)
#### 1-2-2. Quantitative overview
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
#### 3-1-1. Processing Hi-C data from 3DIV database ([1])
- For this case, we need to process raw Hi-C data to get PC1 vectors, since PC1 values are not provided directly.
  - In case of 'stem-closeness/hic-processing/manifest.csv' in this repository, information of Hi-C data obtained from cancer cell lines are written.
  - When processing Hi-C data of other cell line, <ins>name, library\_layout, and run\_accession</ins> information of that data should be written instead in the aforementioned manifest.csv file.
  - Each name in this manifest.csv file should be unique. No two different samples can have identical name in this file.
  - To use the provided Hi-C processing code, library\_layout should be 'paired'. 
    - Library\_layout information can be checked in [SRA](https://www.ncbi.nlm.nih.gov/sra)
```shell
conda activate hic-processing
cd hic-processing
snakemake --cores 100 --resources network=1 --restart-time 3
python ab2corrmat.py
bash copy_corrmat.sh
cd ../
```
#### 3-2-2. Processing PC1 values computed from Hi-C data of normal tissues ([2])
- For this case, PC1 values derived from Hi-C data of normal tissues are provided by the supplementary information. 
```shell
cd data
python3 download_FIRE_PC1.py 
cd ../
```

### 3-2. Conduct stem closeness-related analyses
```shell
conda activate stem-closeness
```
#### 3-2-1. Preprocessing data
- Download TCGA DNA methylation data from UCSC Xena ([3])
```shell
cd data/450k_xena
snakemake -j 10
cd ../../
```
- Construct metadata and bedgraph file of open sea CpG probes in 450K DNA methylation data
  - Download manifest file for Infinium HumanMethylation450 v1.2 BeadChip, provided by Illumina ([4]) 
```shell
cd data
wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv
cd ../
```
  - Make metadata and bedgraph files
```shell
cd utils
bash make-450k-probe-metadata.sh
bash make-450k-probe-bedgraph.sh
cd ../
```
- Construct list of sample names per cohort
```shell
cd utils/scripts
python3 0_get_sample_name_list.py
cd ../../
```
#### 3-2-2. Make binned difference matrices (BDMs)    
##### 3-2-2-1. Make BDMs for TCGA samples ([5]), using open sea CpG probes
```shell
cd binned-difference-matrix-v2-opensea
snakemake -j 10
cd ../
```
##### 3-2-2-2. Make BDMs for PCBC stem cell samples ([6]), using open sea CpG probes
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
  - Representative PC1s are made by averaging PC1 vectors of 10 samples, computed from the BDM of same chromosome. 
```shell
cd utils/scripts
python3 2_make-TCGA-repr-bdm.py
python3 3_make-PAAD-repr-vector-bdm.py
python3 4_make-SC-repr-bdm.py
cd ../../
```
- Compare BDM-derived PC1 and IEBDM-derived PC1
```shell
cd utils/scripts
python3 5_bdm-iebdm-pc1-pcc.py
cd ../../
```
- Download TCGA clinical data.
```shell
cd data
python3 download_TCGA-CDR.py
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
- Download annotation about gene and regulatory features (version: GRCh37)
```shell
cd data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
gzip -d gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
wget https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz
gzip -d homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz
cd ../
```
- Download chromatin state data ([7])
```shell
cd data/roadmap_epigenomics
bash download.sh
cd ../../
```
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
cd ../../../
```

## Reference
[1] Kim, Kyukwang, et al. "3DIV update for 2021: a comprehensive resource of 3D genome and 3D cancer genome." Nucleic Acids Research 49.D1 (2021): D38-D46.         
[2] Schmitt, Anthony D., et al. "A compendium of chromatin contact maps reveals spatially active regions in the human genome." Cell reports 17.8 (2016): 2042-2059.       
[3] Goldman, Mary J., et al. "Visualizing and interpreting cancer genomics data via the Xena platform." Nature biotechnology 38.6 (2020): 675-678.         
[4] Bibikova, Marina, et al. "High density DNA methylation array with single CpG site resolution." Genomics 98.4 (2011): 288-295.            
[5] Hutter, Carolyn, and Jean Claude Zenklusen. "The cancer genome atlas: creating lasting value beyond its data." Cell 173.2 (2018): 283-285.               
[6] Salomonis, Nathan, et al. "Integrated genomic analysis of diverse induced pluripotent stem cells from the progenitor cell biology consortium." Stem cell reports 7.1 (2016): 110-125.        
[7] Kundaje, Anshul, et al. "Integrative analysis of 111 reference human epigenomes." Nature 518.7539 (2015): 317-330.
