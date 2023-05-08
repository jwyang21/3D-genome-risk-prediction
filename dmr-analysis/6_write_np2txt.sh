threshold=mean_std
cpg_type=opensea
dmr_type=HL
wd=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type

echo "python3 6_write_np2txt.py --threshold $threshold -w_dir $wd --dmr_type $dmr_type"
python3 6_write_np2txt.py --threshold $threshold -w_dir $wd --dmr_type $dmr_type
