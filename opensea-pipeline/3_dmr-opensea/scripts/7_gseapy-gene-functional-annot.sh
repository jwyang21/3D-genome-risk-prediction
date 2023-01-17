wd=/data/project/3dith/pipelines/opensea-pipeline/3_dmr-opensea
threshold=mean_std
echo "python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold"
python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold
