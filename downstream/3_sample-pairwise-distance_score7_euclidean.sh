for cohort in TCGA-BLCA TCGA-CESC TCGA-COAD TCGA-GBM TCGA-KIRC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-READ TCGA-SKCM TCGA-THCA TCGA-UCEC TCGA-BRCA TCGA-CHOL TCGA-ESCA TCGA-HNSC TCGA-KIRP TCGA-LUAD TCGA-PAAD TCGA-PRAD TCGA-SARC TCGA-STAD TCGA-THYM
do
	echo "python3 3_sample-pairwise-distance_score7.py --score score7 --cohort $cohort"
	python3 3_sample-pairwise-distance_score7.py --score score7 --cohort $cohort
done
