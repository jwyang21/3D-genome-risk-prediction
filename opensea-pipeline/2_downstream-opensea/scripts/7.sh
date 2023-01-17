cpg_type=opensea
for i in 1 2 3 4 5
do
	python3 7_parse_cox_results.py --version $i --log_dir /data/project/3dith/pipelines/$cpg_type-pipeline/2_downstream-$cpg_type/log --log_fname 7_parse_cox_results-v$i.txt
done
