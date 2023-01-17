CPG_TYPE=island
usage=half
minmax=Y
normalize=N

w_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/2_downstream-$CPG_TYPE
result_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/2_downstream-$CPG_TYPE/result

score=stem-closeness
score_result_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/1_compute-score-$CPG_TYPE/result
score_type=pc1-avg
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	for standardize in Y N
	do
		for distance in cosine-sim euclidean
		do
			for matrix_type in bdm iebdm
			do
				for use_weighted_avg in Y N
				do
					for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
                    do
                        echo "python3 1_survival-analysis.py --num_chrom $num --standardize $standardize -w_dir $w_dir --result_dir $result_dir --usage_option $usage --score_type $score_type --normalize $normalize --minmax $minmax --matrix_type $matrix_type --distance $distance --use_weighted_avg $use_weighted_avg --cohort $cohort --score $score --score_result_dir $score_result_dir"
                        python3 1_survival-analysis.py --num_chrom $num --standardize $standardize -w_dir $w_dir --result_dir $result_dir --usage_option $usage --score_type $score_type --normalize $normalize --minmax $minmax --matrix_type $matrix_type --distance $distance --use_weighted_avg $use_weighted_avg --cohort $cohort --score $score --score_result_dir $score_result_dir
                    done
                done
			done
		done
	done
done
