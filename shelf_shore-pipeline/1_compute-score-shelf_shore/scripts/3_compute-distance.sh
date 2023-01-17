cpg_type=shelf_shore
usage_option=half
w_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type
pc1_upper_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result
result_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result
reference_dir=/data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result
hg19_chr_len=/data/project/3dith/data/hg19.fa.sizes
bdm_bin_dir=/data/project/3dith/data/bdm_bins/$cpg_type
score_type=pc1-avg

for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD PCBC
do
	for standardize in Y N
	do
		for use_weighted_avg in Y N
		do
			for distance in euclidean cosine-sim
			do
				for m_type in iebdm bdm
				do
					for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
                    do
                        echo "python3 3_compute-distance.py --num_chrom $num --standardize $standardize --cohort $cohort --usage_option $usage_option -w_dir $w_dir --pc1_upper_dir $pc1_upper_dir --result_dir $result_dir --reference_dir $reference_dir --hg19_chr_len $hg19_chr_len --bdm_bin_dir $bdm_bin_dir --score_type $score_type --use_weighted_avg $use_weighted_avg --distance $distance --matrix_type $m_type --cpg_type $cpg_type"
                        python3 3_compute-distance.py --num_chrom $num --standardize $standardize --cohort $cohort --usage_option $usage_option -w_dir $w_dir --pc1_upper_dir $pc1_upper_dir --result_dir $result_dir --reference_dir $reference_dir --hg19_chr_len $hg19_chr_len --bdm_bin_dir $bdm_bin_dir --score_type $score_type --use_weighted_avg $use_weighted_avg --distance $distance --matrix_type $m_type --cpg_type $cpg_type
                    done
				done
			done
		done
	done
done
