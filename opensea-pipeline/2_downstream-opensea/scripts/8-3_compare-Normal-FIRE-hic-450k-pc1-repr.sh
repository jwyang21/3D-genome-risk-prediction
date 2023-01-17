cpg_type=opensea
hg19_len_fname=/data/project/jeewon/research/reference/hg19.fa.sizes
binsize=1000000
tcga_fire_cohort_fname=/data/project/3dith/data/etc/tcga-fire-cohorts.csv
hic_pc1_fname=/data/project/3dith/data/fire_pc1.csv

echo "python3 8-3_compare-Normal-FIRE-hic-450k-pc1-repr.py --cpg_type $cpg_type --hg19_len_fname $hg19_len_fname --binsize $binsize --tcga_fire_cohort_fname $tcga_fire_cohort_fname --hic_pc1_fname $hic_pc1_fname"
python3 8-3_compare-Normal-FIRE-hic-450k-pc1-repr.py --cpg_type $cpg_type --hg19_len_fname $hg19_len_fname --binsize $binsize --tcga_fire_cohort_fname $tcga_fire_cohort_fname --hic_pc1_fname $hic_pc1_fname
