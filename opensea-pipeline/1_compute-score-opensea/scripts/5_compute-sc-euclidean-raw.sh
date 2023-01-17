CPG_TYPE=opensea
usage=half
minmax=N
normalize=N
minmax_file_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/1_compute-score-$CPG_TYPE/result
w_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/1_compute-score-$CPG_TYPE
result_dir=/data/project/3dith/pipelines/$CPG_TYPE-pipeline/1_compute-score-$CPG_TYPE/result
distance=euclidean
score_type=pc1-avg

for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD
do
	for standardize in Y N
	do
		for matrix_type in bdm iebdm
		do
			for use_weighted_avg in Y N
			do
				for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
                do
                    echo "python3 5_compute-stem-closeness.py --num_chrom $num --standardize $standardize -w_dir $w_dir --result_dir $result_dir --usage_option $usage --score_type $score_type --normalize $normalize --minmax $minmax --matrix_type $matrix_type --distance $distance --use_weighted_avg $use_weighted_avg --minmax_file_dir $minmax_file_dir --cohort $cohort"
                    python3 5_compute-stem-closeness.py --num_chrom $num --standardize $standardize -w_dir $w_dir --result_dir $result_dir --usage_option $usage --score_type $score_type --normalize $normalize --minmax $minmax --matrix_type $matrix_type --distance $distance --use_weighted_avg $use_weighted_avg --minmax_file_dir $minmax_file_dir --cohort $cohort
                done
            done
		done
	done
done
