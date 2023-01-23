cpg_type=opensea

do
    echo "python3 1_find-bdm-bins.py -w_dir /data/project/3dith/pipelines/utils --tcga_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-pcbc/result --hg19_chr_len /data/project/3dith/data/hg19.fa.sizes --binsize 1000000 --save_dir /data/project/3dith/data/bdm_bins --cpg_type $cpg_type"

    python3 1_find-bdm-bins.py -w_dir /data/project/3dith/pipelines/utils --tcga_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-pcbc/result --hg19_chr_len /data/project/3dith/data/hg19.fa.sizes --binsize 1000000 --save_dir /data/project/3dith/data/bdm_bins --cpg_type $cpg_type
done
