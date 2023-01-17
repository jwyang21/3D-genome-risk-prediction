cpg_type=island
for cohort in TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD PCBC
do
    for m_type in iebdm bdm
    do
        echo "python3 1_all-samples-pc1.py --cohort $cohort -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --tcga_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-pcbc/result -r_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result -m_type $m_type"

        python3 1_all-samples-pc1.py --cohort $cohort -w_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type --tcga_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-v2-$cpg_type/result --pcbc_bdm_dir /data/project/3dith/pipelines/binned-difference-matrix-pcbc/result -r_dir /data/project/3dith/pipelines/$cpg_type-pipeline/1_compute-score-$cpg_type/result -m_type $m_type
    done
done
