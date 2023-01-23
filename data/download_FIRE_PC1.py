import pandas as pd
import numpy as np
import os
# Download PC1 calls derived from Hi-C data of normal tissues, publicly provided by Schmitt, Anthony D., et al. "A compendium of chromatin contact maps reveals spatially active regions in the human genome." Cell reports 17.8 (2016): 2042-2059.

cmd_ = 'wget https://ars.els-cdn.com/content/image/1-s2.0-S2211124716314814-mmc3.xlsx'
os.system(cmd_)

df = pd.read_excel('1-s2.0-S2211124716314814-mmc3.xlsx', sheet_name = 'PCA')
df.to_csv('fire_pc1.csv', index = False)
