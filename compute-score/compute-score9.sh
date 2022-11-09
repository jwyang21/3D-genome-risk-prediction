for cohort in TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD
do
	echo "python3 compute-score9.py --cohort $cohort"
	python3 compute-score9.py --cohort $cohort
done
	
