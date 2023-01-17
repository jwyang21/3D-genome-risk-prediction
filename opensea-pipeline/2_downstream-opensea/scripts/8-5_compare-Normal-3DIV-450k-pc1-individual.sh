cpg_type=opensea
matching_df_fname=/data/project/3dith/data/etc/3div-normal_hic-tcga-matching.csv
hic_corrmat_dir=/data/project/3dith/data/hic_corrmat
hg19_len_fname=/data/project/jeewon/research/reference/hg19.fa.sizes
binsize=1000000
for cohort in TCGA-PAAD TCGA-LUAD TCGA-LUSC
do

	echo "python3 8-5_compare-Normal-3DIV-450k-pc1-individual.py --cpg_type $cpg_type --matching_df_fname $matching_df_fname --hic_corrmat_dir $hic_corrmat_dir --hg19_len_fname $hg19_len_fname --binsize $binsize --cohort $cohort"
	python3 8-5_compare-Normal-3DIV-450k-pc1-individual.py --cpg_type $cpg_type --matching_df_fname $matching_df_fname --hic_corrmat_dir $hic_corrmat_dir --hg19_len_fname $hg19_len_fname --binsize $binsize --cohort $cohort
done

