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
conda activate hic-processing
pip install fanc
```
### 2-2. Installing conda environment for running stem closeness-related analyses
```shell
conda env create --file stem-closeness.yml
```
## 3. Process

### 3-1. Make binned difference matrices (BDMs)
### 3-2. Run opensea pipeline

#### 3-2-1. Compute score
#### 3-2-2. Survival analyses
#### 3-2-3. DMR analyses

### 3-3. Run island pipeline
#### 3-3-1. Compute score
#### 3-3-2. Survival analyses
#### 3-3-3. DMR analyses

### 3-4. Run shelf\_shore pipeline
#### 3-4-1. Compute score
#### 3-4-2. Survival analyses
#### 3-4-3. DMR analyses
