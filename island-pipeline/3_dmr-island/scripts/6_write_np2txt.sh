cpg_type=island
threshold=mean_std
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/3_dmr-$cpg_type

echo "python3 6_write_np2txt.py --threshold $threshold -w_dir $wd"
python3 6_write_np2txt.py --threshold $threshold -w_dir $wd
