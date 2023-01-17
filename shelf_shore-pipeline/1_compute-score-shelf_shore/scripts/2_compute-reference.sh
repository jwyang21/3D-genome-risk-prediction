cpg_type=shelf_shore
score_type=pc1-avg
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD PCBC
do
    for m_type in iebdm bdm
    do
        for standardize in Y N
        do
            echo "python3 2_compute-reference.py --standardize $standardize --cohort $cohort -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --pc1_upper_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result -r_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result --score_type $score_type --usage_option half -m_type $m_type"
            python3 2_compute-reference.py --standardize $standardize --cohort $cohort -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --pc1_upper_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result -r_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result --score_type $score_type --usage_option half -m_type $m_type
        done
    done
done
