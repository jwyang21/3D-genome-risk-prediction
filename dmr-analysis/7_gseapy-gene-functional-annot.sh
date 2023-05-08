cpg_type=opensea
dmr_type=HL
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type
threshold=mean_std
echo "python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold --dmr_type $dmr_type"
python3 7_gseapy-gene-functional-annot.py -w_dir $wd --threshold $threshold --dmr_type $dmr_type
