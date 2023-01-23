cpg_type=opensea
matching_df_fname=/data/project/3dith/data/etc/3div-cancer_hic-tcga-matching.csv
hic_corrmat_dir=/data/project/3dith/data/hic_corrmat
hg19_len_fname=/data/project/jeewon/research/reference/hg19.fa.sizes
binsize=1000000
m=bdm

for cohort in TCGA-BRCA TCGA-COAD TCGA-ESCA TCGA-KIRC TCGA-LIHC TCGA-LUSC TCGA-PAAD TCGA-PRAD
do

    echo "python3 compare-Tumor-hic-450k-pc1-individual-plot.py --cpg_type $cpg_type --matching_df_fname $matching_df_fname --hic_corrmat_dir $hic_corrmat_dir --hg19_len_fname $hg19_len_fname --binsize $binsize --cohort $cohort --matrix_type $m"
    python3 compare-Tumor-hic-450k-pc1-individual-plot.py --cpg_type $cpg_type --matching_df_fname $matching_df_fname --hic_corrmat_dir $hic_corrmat_dir --hg19_len_fname $hg19_len_fname --binsize $binsize --cohort $cohort --matrix_type $m
    
done
