import pandas as pd
import numpy as np
import os

ALL_TCGA = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD'.split(' ')

ALL_TCGA_PCBC = 'TCGA-LGG TCGA-UCS TCGA-BLCA TCGA-LUAD TCGA-THYM TCGA-PRAD TCGA-DLBC TCGA-ACC TCGA-KICH TCGA-GBM TCGA-READ TCGA-KIRC TCGA-LAML TCGA-ESCA TCGA-STAD TCGA-UCEC TCGA-KIRP TCGA-OV TCGA-SARC TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-PCPG TCGA-SKCM TCGA-TGCT TCGA-CESC TCGA-CHOL TCGA-PAAD TCGA-UVM TCGA-MESO TCGA-BRCA TCGA-COAD PCBC'.split(' ')

if __name__ == '__main__':
    all_samples = {}
    for cohort in ALL_TCGA_PCBC:
        if 'TCGA' in cohort:
            cohort_dir = os.path.join('/data/project/3dith/data/450k_bdgs_v2/', cohort+'-opensea')
        else: 
            cohort_dir = '/data/project/3dith/data/pcbc/bdgs_opensea'
        
        all_samples[cohort] = []
        
        for f in os.listdir(cohort_dir):
            if f.endswith('.sorted.bedGraph'):
                current_sample = f.split('.sorted.bedGraph')[0]
                all_samples[cohort].append(current_sample)
            else:
                pass
    
    result_fname = '/data/project/3dith/data/samplenames'
    np.savez(result_fname, **all_samples)
    print("result file: {}.npz".format(result_fname))
