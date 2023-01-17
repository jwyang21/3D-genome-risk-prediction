binsize=1000000
cpg_type=shelf_shore
echo "python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/shelf_shore-pipeline/3_dmr-shelf_cpg --cpg_type $cpg_type"
python3 2_find_DMR.py --binsize $binsize -w_dir /data/project/3dith/pipelines/shelf_shore-pipeline/3_dmr-shelf_shore --cpg_type $cpg_type
