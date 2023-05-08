binsize=1000000
cpg_type=opensea
dmr_type=HL

echo "python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type --cpg_type $cpg_type --dmr_type $dmr_type"
python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type --cpg_type $cpg_type --dmr_type $dmr_type
