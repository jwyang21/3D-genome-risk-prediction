binsize=1000000
cpg_type=opensea
echo "python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/opensea-pipeline/3_dmr-opensea --cpg_type $cpg_type"
python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/opensea-pipeline/3_dmr-opensea --cpg_type $cpg_type
