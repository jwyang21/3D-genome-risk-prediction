# FIRE COHORTS 
# cohrots which have both FIRE calls and 450k data.
# distance from FIRE (normal) reference
# use integration of abs(pc1) of FIRE and score3
for cohort in TCGA-BLCA TCGA-LUAD TCGA-ACC TCGA-OV TCGA-LIHC TCGA-LUSC TCGA-PAAD
do
	echo "python3 compute-score9.py --cohort $cohort"
	python3 compute-score9.py --cohort $cohort
done
	
