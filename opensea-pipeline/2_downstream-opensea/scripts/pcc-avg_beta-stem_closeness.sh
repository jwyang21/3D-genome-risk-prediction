cpg_type=opensea
w_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/2_downstream-$cpg_type
score_f=/data/project/3dith/data/cohort-1-best-score-km.csv

for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD TCGA-ESCA TCGA-HNSC

do
	echo "python3 pcc-avg_beta-stem_closeness.py -w_dir $w_dir --cohort $cohort --cpg_type $cpg_type --score_fname $score_f"
	python3 pcc-avg_beta-stem_closeness.py -w_dir $w_dir --cohort $cohort --cpg_type $cpg_type --score_fname $score_f
done
