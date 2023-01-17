threshold=mean_std
wd=/data/project/3dith/pipelines/opensea-pipeline/3_dmr-opensea

echo "python3 6_write_np2txt.py --threshold $threshold -w_dir $wd"
python3 6_write_np2txt.py --threshold $threshold -w_dir $wd
