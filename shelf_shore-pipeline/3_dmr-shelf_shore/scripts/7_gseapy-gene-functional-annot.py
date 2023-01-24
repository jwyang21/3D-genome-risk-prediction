import numpy as np
import pandas as pd
import os
import gseapy as gp
from pybiomart import Dataset
import argparse
dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
res = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'])

ensg2symbol = {r['Gene stable ID']:r['Gene name'] for r in res[res['Gene type'] == 'protein_coding'].to_records()}

NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

SMALL_CATEGORY_FNAME = '/data/project/3dith/data/etc/dmr-feature-small-category.npz'
if os.path.exists(SMALL_CATEGORY_FNAME):
    SMALL_CATEGORY = np.load(SMALL_CATEGORY_FNAME, allow_pickle = True)
else:
    SMALL_CATEGORY = {}
    SMALL_CATEGORY['GENE'] = ['gene','transcript']
    SMALL_CATEGORY['REG'] = ['open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', 'enhancer', 'promoter', 'promoter_flanking_region']
    SMALL_CATEGORY['EPI'] = ['18_Quies', '13_Het', '17_ReprPCWk', '16_ReprPC', '14_TssBiv', '2_TssFlnk','12_ZNF/Rpts', '11_EnhWk', '1_TssA', '6_TxWk', '5_Tx', '9_EnhA1', '7_EnhG1',     '4_TssFlnkD' ,'15_EnhBiv', '10_EnhA2', '3_TssFlnkU', '8_EnhG2']
    np.savez(SMALL_CATEGORY_FNAME,**SMALL_CATEGORY)
    
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
cohort2eid = pd.read_csv('/data/project/3dith/data/etc/cohort2eid.txt', sep = '\t', header = None)
cohort2eid.columns = ['cohort', 'eid']
eid_cohorts = cohort2eid.cohort.values
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]

def parse_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    args.add_argument('--threshold', help = 'threshold for DMR', type = str, required = True)
    return args.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    os.chdir(args.working_dir)
    ensg_fname = 'ENSG_'+args.threshold+'.txt'

    all_df = pd.DataFrame(np.zeros((len(NORMAL7_COHORT), 1), dtype = int), index = NORMAL7_COHORT, columns = [args.threshold])
    item_list = []
    for cohort in NORMAL7_COHORT:
        print("===\n{}".format(cohort))
        cohort_dir = os.path.join(os.getcwd(), 'result', cohort)

        ensg = pd.read_csv(os.path.join(cohort_dir, ensg_fname), sep = '\t', header = None)

        all_df.loc[cohort][args.threshold] = ensg.shape[0]
        current_item = cohort+'_'+args.threshold+'_'+'gene_symbol'
        if current_item not in item_list:
            item_list.append(current_item)
        globals()[cohort+'_'+args.threshold+'_'+'gene_symbol'] = []

        for g in ensg.values.flatten():
            if g in list(ensg2symbol.keys()):
                if str(ensg2symbol[g]).lower() != 'nan':
                    globals()[cohort+'_'+args.threshold+'_'+'gene_symbol'].append(ensg2symbol[g])

    all_df_fname = os.path.join(os.getcwd(), 'result', 'Tumor-Normal-DMR-genes-num.csv')
    all_df.to_csv(all_df_fname, index = True)
    print("num_DMR_genes_file (tumor-normal): {}".format(all_df_fname))

    all_item_dictionary={}
    for item_ in item_list:
        all_item_dictionary[item_] = globals()[item_]
    result_fname = os.path.join(os.getcwd(), 'result', 'Tumor-Normal-DMR-genes-GeneSymbols')
    np.savez(result_fname, **all_item_dictionary)
    print(result_fname+'.npz')

    all_df = pd.DataFrame(np.zeros((len(NORMAL7_COHORT), 1), dtype = int), index = NORMAL7_COHORT, columns = [args.threshold])
    print(item_list)
    print("total {} items".format(len(item_list)))
    for item_ in item_list:
        cohort = item_.split('_')[0].strip()
        cohort_dir = os.path.join(os.getcwd(), 'result', cohort)
        print("---\n{}".format(item_))
        print("num_DMR_genes: {}".format(len(globals()[item_])))
        result = gp.enrichr(globals()[item_], 'GO_Biological_Process_2015').res2d.sort_values('Adjusted P-value')
        result_sig = result[result['Adjusted P-value'] < 5e-2].copy()
        print("num_significant_genes: {}".format(result_sig.shape[0]))
        result_fname = os.path.join(cohort_dir, 'DMR-gene-GSEAPY-'+args.threshold+'.csv')
        all_df.loc[cohort][args.threshold] = result_sig.shape[0]
        print("significant_DMR_genes: {}".format(result_fname))
        result_sig.to_csv(result_fname)
    all_df_fname = os.path.join(os.getcwd(), 'result', 'Tumor-Normal-DMR-genes-Significant-Genes.csv')
    all_df.to_csv(all_df_fname)
    print("===\n"+all_df_fname)
