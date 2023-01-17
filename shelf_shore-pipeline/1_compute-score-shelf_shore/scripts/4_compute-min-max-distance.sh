cpg_type=shelf_shore
manifest=/data/project/3dith/data/manifest_normal7.csv
echo "python3 4_compute-min-max-distance.py --manifest_fname $manifest -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --score_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result/ --result_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result/ --usage_option half"
python3 4_compute-min-max-distance.py --manifest_fname $manifest -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --score_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result/ --result_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result/ --usage_option half
