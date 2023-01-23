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
```

### 3-2. Conduct stem closeness-related analyses
```shell
conda activate stem-closeness
```
### 3-2-1. Make binned difference matrices (BDMs)    
#### 3-2-1-1. Make BDMs for TCGA samples, using open sea CpG probes
```shell
cd binned-difference-matrix-v2-opensea
snakemake -j 10
cd ../
```
#### 3-2-1-2. Make BDMs for PCBC stem cell samples, using open sea CpG probes
```shell
cd binned-difference-matrix-pcbc
snakemake -j 10
cd ../
```

### 3-2-2. Run opensea pipeline
```shell
bash opensea-pipeline/1_compute-score-opensea/scripts/1_all-samples-pc1.sh
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
