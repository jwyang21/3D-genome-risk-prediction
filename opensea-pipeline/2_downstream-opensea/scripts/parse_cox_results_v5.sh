cpg_type=opensea
log_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/2_downstream-$cpg_type/log
version=5
log_fname=parse_cox_results-v$version.txt

echo "python3 parse_cox_results.py --version $version --log_dir $log_dir --log_fname $log_fname"
python3 parse_cox_results.py --version $version --log_dir $log_dir --log_fname $log_fname

