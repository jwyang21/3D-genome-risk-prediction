binsize=1000000
cpg_type=island
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 1_compute-binned-island-beta-TN-diff.py -b $binsize -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/3_dmr-$cpg_type -c $cohort --cpg_metadata_fname /data/project/3dith/data/450k_metadata.$cpg_type.sorted.bed --cpg_type $cpg_type --dmr_type TN"
	python3 1_compute-binned-island-beta-TN-diff.py -b $binsize -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/3_dmr-$cpg_type -c $cohort --cpg_metadata_fname /data/project/3dith/data/450k_metadata.$cpg_type.sorted.bed --cpg_type $cpg_type --dmr_type TN
done
