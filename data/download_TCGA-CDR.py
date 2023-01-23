import pandas as pd
import numpy as np
import os

# download TCGA clinical data from https://gdc.cancer.gov/about-data/publications/pancanatlas
cmd_ = 'wget https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81 -O TCGA-CDR.xlsx'

os.system(cmd_)

df = pd.read_excel('TCGA-CDR.xlsx', sheet_name = 'TCGA-CDR', index_col = 0)

df.to_csv('TCGA-CDR-SupplementalTableS1.csv', index = False)
