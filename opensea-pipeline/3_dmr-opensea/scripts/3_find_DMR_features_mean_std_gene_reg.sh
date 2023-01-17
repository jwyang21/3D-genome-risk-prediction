binsize=1000000
threshold=mean_std
wd=/data/project/3dith/pipelines/opensea-pipeline/3_dmr-opensea
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 3_find_DMR_features.py --cohort $cohort --binsize $binsize --threshold $threshold --dmr_feature gene_reg -w_dir $wd"
	python3 3_find_DMR_features.py --cohort $cohort --binsize $binsize --threshold $threshold --dmr_feature gene_reg -w_dir $wd
done
