for score_type in avg_pc1 pc1_fluctuation
do
	echo "compute-score2-hold-out.py --cohort PCBC --reference_type PCBC --score_type $score_type"
	python3 compute-score2-hold-out.py --cohort PCBC --reference_type PCBC --score_type $score_type	
done
