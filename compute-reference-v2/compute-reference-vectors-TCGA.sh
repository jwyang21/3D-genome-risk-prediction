for cohort in TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD
do
	for score_type in avg_pc1 pc1_fluctuation
	do
		echo "python3 compute-reference-vectors.py --cohort $cohort --reference_type TCGA --score_type $score_type"
		python3 compute-reference-vectors.py --cohort $cohort --reference_type TCGA --score_type $score_type
	done
done
