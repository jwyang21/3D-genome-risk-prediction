# To-do: compute binned difference matrix using CpG probes which (1) belong to PMD and (2) has notnull beta value
for cohort in TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD
do
	for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
	do
		echo "python3 beta-pmd-binned-diffmat.py --cohort $cohort --chrom $chrom --binsize 1000000"
		python3 beta-pmd-binned-diffmat.py --cohort $cohort --chrom $chrom --binsize 1000000
	done
done
