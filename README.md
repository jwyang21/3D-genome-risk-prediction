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
#### 3-2-1. Construct metadata and bedgraph file ofopen sea CpG probes in 450K DNA methylation data
- Download manifest file for Infinium HumanMethylation450 v1.2 BeadChip, provided by Illumina ([1]) 
```shell
cd data
wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv
cd ../

python3 utils/make-450k-probe-metadata.py --cpg_type opensea
python3 utils/make-450k-probe-bedgraph.py --cpg_type opensea
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
```shell
cd opensea-pipeline/1_compute-score-opensea/scripts
```
- Conduct log-rank test.
```shell
cd opensea-pipeline/2_downstream-opensea/scripts/
bash 1_km-sc-cosine_and_euclidean-minmax.sh > ../log/1_km-sc-cosine_and_euclidean-minmax.log
bash 1_km-sc-euclidean-normalized.sh > ../log/1_km-sc-euclidean-normalized.log
bash 1_km-sc-euclidean-raw.sh > ../log/1_km-sc-euclidean-raw.log
```
- Figure out correlation between stem closeness and average open sea DNA methylation level.
```shell
bash 2_pcc-avg_beta-stem_closeness.sh > ../log/2_pcc-avg_beta-stem_closeness-ALL.log
bash 3_parse-avg_meth-score-corr-log.sh
- Write commands to run Cox regression
```
```python
python3 5_write-cox-commands.py --cpg_type opensea --cox_version 5 --command_fname 6_cox-commands --score_fname /data/project/3dith/data/cohort-1-best-score-km.csv
```
- Run Cox regression
```shell
bash 6_cox_commands_v5.sh > ../log/6_cox_commands_v5.log
```
- Parse Cox regression results
```shell
bash 7_parse_cox_results_v5.sh
```


#### 3-2-2-1. Compute score
#### 3-2-2-2. Survival analyses
#### 3-2-2-3. DMR analyses

### 3-4. Run island pipeline
#### 3-4-1. Compute score
#### 3-4-2. Survival analyses
#### 3-4-3. DMR analyses

### 3-5. Run shelf\_shore pipeline
#### 3-5-1. Compute score
#### 3-5-2. Survival analyses
#### 3-5-3. DMR analyses

## Reference
[1] Bibikova, Marina, et al. "High density DNA methylation array with single CpG site resolution." Genomics 98.4 (2011): 288-295.     
