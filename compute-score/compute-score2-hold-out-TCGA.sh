for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	for score_type in avg_pc1 pc1_fluctuation
	do
		echo "compute-score2-score4-hold-out.py --cohort $cohort --reference_type TCGA --score_type $score_type"
		python3 compute-score2-score4-hold-out.py --cohort $cohort --reference_type TCGA --score_type $score_type
	done	
done
