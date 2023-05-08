binsize=1000000
cpg_type=opensea
dmr_type=HL
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type
#for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 4_collate_DMR_features.py --cohort $cohort --binsize $binsize -w_dir $wd --dmr_type $dmr_type"
	python3 4_collate_DMR_features.py --cohort $cohort --binsize $binsize -w_dir $wd --dmr_type $dmr_type
done
