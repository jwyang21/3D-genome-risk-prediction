binsize=1000000
wd=/data/project/3dith/pipelines/island-pipeline/3_dmr-island
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 4_collate_DMR_features.py --cohort $cohort --binsize $binsize -w_dir $wd"
	python3 4_collate_DMR_features.py --cohort $cohort --binsize $binsize -w_dir $wd
done
