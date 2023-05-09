# A deep learning-based risk prediction using 3D genome-aware epigenetic features
## 1. Introduction and overview
### 1-1. Introduction    
We propose the utilization of 3D genome-aware epigenetic features extracted from the single-sample 450K DNA methylation data.       
The overall pipeline consists of four main steps: (1) inferring 3D genome-aware epigenetic features from the single-sample 450K DNA methylation data (note that inferred 3D genome state is representated as the first principal component; PC1), (2) construction of stem/normal references by averaging the PC1s from multiple samples (normal samples from the identical tissue type or the stem cells), (3) using combinations of distances between PC1s (i.e., vector representation of 3D genome states inferred from the 450K DNA methylation data) to quantify stem closeness of each sample, and (4) risk prediction by the feedforward neural network, using concatenation of 3D genome-aware epigenetic features and survival-related features (age and gender) as input feature.              
The figures in section 1-2 describe a more detailed overview of the whole pipeline. Subplots in the section 1-2-1 illustrate the aforementioned main steps, and section 1-2-2 provides the quantitative explanation of the procedures for extracting the 3D genome-aware epigenetic features from the single-sample 450K DNA methylation data.             
### 1-2. Overview
#### 1-2-1. Qualitative overview
![main fig1  230424_fig1-300dpi](https://user-images.githubusercontent.com/86412887/236715096-3a22e8d1-6d32-4bab-ae55-d401da73d2b9.png)
#### 1-2-2. Quantitative overview
![230423-FigS2-stem-closeness-explanation-300dpi](https://user-images.githubusercontent.com/86412887/236715225-0e5239c7-a5e6-4990-87e2-53ab1bde38c2.png)
## 2. Installation of the conda environments             
There are three environments needed: environments for (1) processing Hi-C data, (2) extracting 3D genome-aware epigenetic features from the 450K DNA methylation data, and (3) survival analysis and DMR analysis.     
### 2-1. Installing conda environment processing Hi-C data
```shell
conda install mamba -n base -c conda-forge
mamba env create --file hic-processing/environment.yaml
```
### 2-2. Installing conda environment for extracting 3D genome-aware epigenetic features
```shell
conda install mamba -n base -c conda-forge
mamba env create --file stem-closeness.yaml
```
### 2-3. Installing conda environment for survival analysis and DMR analysis
```shell
conda install mamba -n base -c conda-forge
mamba env create --file survival-analysis.yaml
python3 -m pip install --upgrade https://github.com/Lasagne/Lasagne/archive/master.zip
git clone https://github.com/jaredleekatzman/DeepSurv.git [destination_folder]
```
In the \[destination\_folder] where the git repository of DeepSurv ([1]) is cloned to, replace the 'deep_surv.py' file in the 'deepsurv' folder with the [deep_surv.py](https://github.com/jwyang21/3D-genome-risk-prediction/blob/main/survival-analysis/deep_surv.py) which is located in the 'survival-analysis' directory of the current repository.
## 3. Processing data and running experiments

### 3-1. Processing Hi-C data 
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
#### 3-1-2. Processing PC1 values computed from Hi-C data of normal tissues ([2])
- For this case, PC1 values derived from Hi-C data of normal tissues are provided by the supplementary information. 
```shell
cd data
python3 download_FIRE_PC1.py 
cd ../
```
### 3-2. Extracting 3D genome-aware epigenetic features from binned difference matrices (BDM)
```shell
conda activate stem-closeness
```

#### 3-2-1. Constructing BDM  
##### 3-2-1-1. Constructing BDMs for TCGA samples ([5]), using open sea CpG probes
```shell
cd binned-difference-matrix-v2-opensea
snakemake -j 10
cd ../
```
##### 3-2-1-2. Constructing BDMs for PCBC stem cell samples ([6]), using open sea CpG probes
```shell
cd binned-difference-matrix-pcbc
snakemake -j 10
cd ../
```
##### 3-2-1-3. Constructing list of genomic bins used for constructing BDM.
```shell
cd utils/scripts
bash 1_find-bdm-bins.sh
cd ../../
```

#### 3-2-2. Extracting 3D genome-aware epigenetic features: BDM PC1s, stem/normal references, stem/normal distances, and stem closeness
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

### 3-3. Investigating the characteristics of BDM 
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
### 3-4. Survival analysis
```shell
conda activate survival-analysis
```        
#### 3-4-1. Risk prediction using a feedforward neural network, followed by the log-rank test       
```shell
cd survival-analysis/
bash risk_prediction.sh
```
#### 3-4-2. Cox regression using the predicted risks
```shell
bash cox_risk.sh
cd ../
```
### 3-5. DMR analysis
```shell
cd dmr-analysis/
bash 0_compute_binned_avg_opensea_beta.sh
bash 1_compute-binned-opensea-beta-HL-diff.sh
bash 2_find_DMR.sh
bash 3_find_DMR_features_epi.sh
bash 3_find_DMR_features_gene_reg.sh
bash 4_collate_DMR_features.sh
bash 5_make-cohort-threshold-feature-list.sh
bash 6_write_np2txt.sh
bash 7_gseapy-gene-functional-annot.sh
bash 8_compute-chromatin-state-proportion.sh
cd ../
```
## Reference
[1] Kim, Kyukwang, et al. "3DIV update for 2021: a comprehensive resource of 3D genome and 3D cancer genome." Nucleic Acids Research 49.D1 (2021): D38-D46.         
[2] Schmitt, Anthony D., et al. "A compendium of chromatin contact maps reveals spatially active regions in the human genome." Cell reports 17.8 (2016): 2042-2059.       
[3] Goldman, Mary J., et al. "Visualizing and interpreting cancer genomics data via the Xena platform." Nature biotechnology 38.6 (2020): 675-678.         
[4] Bibikova, Marina, et al. "High density DNA methylation array with single CpG site resolution." Genomics 98.4 (2011): 288-295.            
[5] Hutter, Carolyn, and Jean Claude Zenklusen. "The cancer genome atlas: creating lasting value beyond its data." Cell 173.2 (2018): 283-285.               
[6] Salomonis, Nathan, et al. "Integrated genomic analysis of diverse induced pluripotent stem cells from the progenitor cell biology consortium." Stem cell reports 7.1 (2016): 110-125.        
[7] Kundaje, Anshul, et al. "Integrative analysis of 111 reference human epigenomes." Nature 518.7539 (2015): 317-330.
