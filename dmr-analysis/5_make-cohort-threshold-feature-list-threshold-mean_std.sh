binsize=1000000
threshold=mean_std
cpg_type=opensea
dmr_type=HL
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type
#for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 5_make-cohort-threshold-feature-list.py --cohort $cohort --binsize $binsize --threshold $threshold -w_dir $wd --dmr_type $dmr_type"
	python3 5_make-cohort-threshold-feature-list.py --cohort $cohort --binsize $binsize --threshold $threshold -w_dir $wd --dmr_type $dmr_type
done
