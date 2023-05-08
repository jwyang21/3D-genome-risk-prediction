cpg_type=opensea
dmr_type=HL
working_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type
save_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/4_$dmr_type-DMR-$cpg_type/result
cohort2eid=/data/project/3dith/data/etc/cohort2eid_score.txt #/data/project/3dith/data/etc/cohort2eid.txt

echo "python3 8_compute-chromatin-state-proportion.py --cpg_type $cpg_type --cohort2eid $cohort2eid --chromatin_states /data/project/3dith/data/chromatin_states.npy --input_fname DMR_EPI_features_threshold_mean_std_len.npz --dmr_type $dmr_type --working_dir $working_dir --save_dir $save_dir > ../log/8_compute-chromatin-state-proportion.log"
python3 8_compute-chromatin-state-proportion.py --cpg_type $cpg_type --cohort2eid $cohort2eid --chromatin_states /data/project/3dith/data/chromatin_states.npy --input_fname DMR_EPI_features_threshold_mean_std_len.npz --dmr_type $dmr_type --working_dir $working_dir --save_dir $save_dir > ../log/8_compute-chromatin-state-proportion.log
