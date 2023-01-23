CPG_TYPE=opensea
score_name=avg_beta
score_result_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/2_downstream-$CPG_TYPE/result
score_fname=opensea_tumors_avg_beta.csv
w_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/2_downstream-$CPG_TYPE
result_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/2_downstream-$CPG_TYPE/result
lowest_save_dir=kaplan-meier-$score_name-v2


for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	echo "python3 survival-analysis-avg_beta.py --score_name $score_name --cohort $cohort --score_result_dir $score_result_dir --score_fname $score_fname -w_dir $w_dir --result_dir $result_dir --lowest_save_dir $lowest_save_dir"
	python3 survival-analysis-avg_beta.py --score_name $score_name --cohort $cohort --score_result_dir $score_result_dir --score_fname $score_fname -w_dir $w_dir --result_dir $result_dir --lowest_save_dir $lowest_save_dir
done
