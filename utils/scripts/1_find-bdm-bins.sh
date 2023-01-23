cpg_type=opensea

echo "python3 utils/scripts/1_find-bdm-bins.py -w_dir utils --tcga_bdm_dir binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir binned-difference-matrix-pcbc/result --hg19_chr_len ../data/hg19.fa.sizes --binsize 1000000 --save_dir data/bdm_bins --cpg_type $cpg_type"

python3 utils/scripts/1_find-bdm-bins.py -w_dir utils --tcga_bdm_dir binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir binned-difference-matrix-pcbc/result --hg19_chr_len ../data/hg19.fa.sizes --binsize 1000000 --save_dir data/bdm_bins --cpg_type $cpg_type
