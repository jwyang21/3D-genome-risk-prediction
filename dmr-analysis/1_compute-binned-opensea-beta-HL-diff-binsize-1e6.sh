binsize=1000000
cpg_type=opensea
w_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/5_risk-HL-DMR-$cpg_type
cpg_metadata_fname=/data/project/3dith/data/450k_metadata.open_sea.sorted.bed

version=v7.1 # optimal feature version for TCGA cohorts
lr=0.001 # optimal lr for TCGA cohorts

dmr_type=risk_HL
for event in OS DSS DFI PFI
do
	for fold in fold1 fold2 fold3 fold4 fold5
	do
		for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
		do
			echo "python3 1_compute-binned-opensea-beta-HL-diff.py -b $binsize -w_dir $w_dir --cpg_metadata_fname $cpg_metadata_fname --cpg_type $cpg_type --version $version --lr $lr --dmr_type $dmr_type -c $cohort --event $event --fold $fold"
			python3 1_compute-binned-opensea-beta-HL-diff.py -b $binsize -w_dir $w_dir --cpg_metadata_fname $cpg_metadata_fname --cpg_type $cpg_type --version $version --lr $lr --dmr_type $dmr_type -c $cohort --event $event --fold $fold
		done
	done
done
