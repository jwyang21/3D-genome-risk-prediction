cpg_type=shelf_shore
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/3_dmr-$cpg_type
threshold=mean_std
echo "python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold"
python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold
