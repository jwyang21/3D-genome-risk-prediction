import numpy as np
import pandas as pd
import os
import gseapy as gp
from pybiomart import Dataset
import argparse
import csv

dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
res = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'])

ensg2symbol = {r['Gene stable ID']:r['Gene name'] for r in res[res['Gene type'] == 'protein_coding'].to_records()}
SCORE_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')

SMALL_CATEGORY_FNAME = '/data/project/3dith/data/etc/dmr-feature-small-category.npz'
if os.path.exists(SMALL_CATEGORY_FNAME):
    SMALL_CATEGORY = np.load(SMALL_CATEGORY_FNAME, allow_pickle = True)
else:
    SMALL_CATEGORY = {}
    SMALL_CATEGORY['GENE'] = ['gene','transcript']
    SMALL_CATEGORY['REG'] = ['open_chromatin_region', 'TF_binding_site', 'CTCF_binding_site', 'enhancer', 'promoter', 'promoter_flanking_region']
    SMALL_CATEGORY['EPI'] = ['18_Quies', '13_Het', '17_ReprPCWk', '16_ReprPC', '14_TssBiv', '2_TssFlnk',     '12_ZNF/Rpts', '11_EnhWk', '1_TssA', '6_TxWk', '5_Tx', '9_EnhA1', '7_EnhG1',     '4_TssFlnkD' ,'15_EnhBiv', '10_EnhA2', '3_TssFlnkU', '8_EnhG2']
    np.savez(SMALL_CATEGORY_FNAME,**SMALL_CATEGORY)
    
BIG_CATEGORY = ['GENE', 'REG', 'EPI']
cohort2eid = pd.read_csv('/data/project/3dith/data/etc/cohort2eid.txt', sep = '\t', header = None)
cohort2eid.columns = ['cohort', 'eid']
eid_cohorts = cohort2eid.cohort.values
CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w_dir', '--working_dir', help = 'working directory', type = str, required = True)
    parser.add_argument('--threshold', help = 'threshold for DMR', type = str, required = True)
    parser.add_argument('--dmr_type', type = str, help = 'TN or HL', required = True)
    
    parser.add_argument('--event', help = 'survival event', type = str, required = True)
    parser.add_argument('--fold', help = 'fold number', type = str, required = True, default = 'fold1')
    parser.add_argument('--version', help = 'feature version', type = str, required = True, default = 'v7.1')
    parser.add_argument('--lr', help = 'learning rate. 0.001, 0.0005, 0.0001', type = float, required = True, default = 0.001)
    
    parser.add_argument('-c', '--cohort', help = 'TCGA cohort', type = str, required = True)
    return parser.parse_args()

if __name__=='__main__':
    args = parse_arguments()
    
    os.chdir(args.working_dir)
    
    if args.dmr_type != 'risk_HL':
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort)
    else:
        SAVEDIR = os.path.join(os.getcwd(), 'result', args.cohort, f'{args.version}_lr_{args.lr}_{args.event}_{args.fold}')
    SAVEDIR_basename = os.path.basename(SAVEDIR)
        
    ensg_fname = 'ENSG_'+args.threshold+'.txt'
    
    result_dir = os.path.join(os.getcwd(), 'result')
    
    if args.dmr_type == 'TN':
        all_df_fname = os.path.join(result_dir, 'Tumor-Normal-DMR-genes-num.csv')
    elif args.dmr_type == 'HL':
        all_df_fname = os.path.join(result_dir, 'High-Low-DMR-genes-num.csv')
    elif args.dmr_type == 'risk_HL':
        all_df_fname = os.path.join(result_dir, 'Risk_High-Low-DMR-genes-num.csv')
    else:
        pass
    
    if os.path.exists(all_df_fname):
        f_ = open(all_df_fname, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        f_ = open(all_df_fname, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        f_writer.writerow(['cohort', 'category', 'mean_std'])
      
    item_list = []
    
    print("===\n{}".format(args.cohort))
    
    ensg = pd.read_csv(os.path.join(SAVEDIR, ensg_fname), sep = '\t', header = None)    
    contents = [args.cohort, SAVEDIR_basename, ensg.shape[0]]
    f_writer.writerow(contents)
    
    current_item = args.cohort+'_'+args.threshold+'_'+'gene_symbol'
    if current_item not in item_list:
        item_list.append(current_item)
    globals()[args.cohort+'_'+args.threshold+'_'+'gene_symbol'] = []

    for g in ensg.values.flatten():
        if g in list(ensg2symbol.keys()):
            if str(ensg2symbol[g]).lower() != 'nan':
                globals()[args.cohort+'_'+args.threshold+'_'+'gene_symbol'].append(ensg2symbol[g])
    f_.close()

    print("num_DMR_genes_file: {}".format(all_df_fname))

    all_item_dictionary={}
    for item_ in item_list:
        all_item_dictionary[item_] = globals()[item_]
    if args.dmr_type == 'TN':
        result_fname = os.path.join(os.getcwd(), 'result', 'Tumor-Normal-DMR-genes-GeneSymbols')
    elif args.dmr_type == 'HL':
        result_fname = os.path.join(os.getcwd(), 'result', 'High-Low-DMR-genes-GeneSymbols')
    elif args.dmr_type == 'risk_HL':
        result_fname = os.path.join(os.getcwd(), 'result', 'Risk_High-Low-DMR-genes-GeneSymbols')
    else:
        pass
    
    np.savez(result_fname, **all_item_dictionary)
    print(result_fname+'.npz')

    result_dir = os.path.join(os.getcwd(), 'result')
    if args.dmr_type == 'TN':
        all_df_fname = os.path.join(result_dir, 'Tumor-Normal-DMR-genes-Significant-Genes.csv')
    elif args.dmr_type == 'HL':
        all_df_fname = os.path.join(result_dir, 'High-Low-DMR-genes-Significant-Genes.csv')
    elif args.dmr_type == 'risk_HL':
        all_df_fname = os.path.join(result_dir, 'Risk_High-Low-DMR-genes-Significant-Genes.csv')
    else:
        pass    
            
    if os.path.exists(all_df_fname):
        f_ = open(all_df_fname, 'a', encoding = 'utf-8')
        f_writer = csv.writer(f_)
    else:
        f_ = open(all_df_fname, 'w', encoding = 'utf-8')
        f_writer = csv.writer(f_)
        f_writer.writerow(['cohort', 'category', 'mean_std'])   
    
    print(item_list)
    print("total {} items".format(len(item_list)))
    for item_ in item_list:
        cohort = item_.split('_')[0].strip()
        print("---\n{}".format(item_))
        print("num_DMR_genes: {}".format(len(globals()[item_])))
        result = gp.enrichr(globals()[item_], 'GO_Biological_Process_2015').res2d.sort_values('Adjusted P-value')
        result_sig = result[result['Adjusted P-value'] < 5e-2].copy()
        print("num_significant_genes: {}".format(result_sig.shape[0]))
        result_fname = os.path.join(SAVEDIR, 'DMR-gene-GSEAPY-'+args.threshold+'.csv')
        print("significant_DMR_genes: {}".format(result_fname))        
        content = [cohort, SAVEDIR_basename, result_sig.shape[0]]
        f_writer.writerow(content)        
        result_sig.to_csv(result_fname)

    f_.close()
    print("===\n"+all_df_fname)
