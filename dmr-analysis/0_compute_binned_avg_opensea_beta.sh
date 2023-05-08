cpg_type=opensea
w_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/5_risk-HL-DMR-$cpg_type
log_fname=0_compute_binned_avg_opensea_beta.log

echo "python3 0_compute_binned_avg_opensea_beta.py -w_dir $w_dir --cpg_type $cpg_type > ../log/$log_fname"
python3 0_compute_binned_avg_opensea_beta.py -w_dir $w_dir --cpg_type $cpg_type > ../log/$log_fname
