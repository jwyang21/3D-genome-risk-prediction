#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os

def process_chrom(c):
    if isinstance(c, float):
        return f'chr{int(c)}'
    elif isinstance(c, int):
        return f'chr{c}'
    else:
        return 'chr'+c

meta = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = 7)

raw_meta = meta.copy()

cpg_type = 'opensea'
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & meta.Relation_to_UCSC_CpG_Island.notnull()].copy()
print(globals()['meta_notnull_'+cpg_type].shape)
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))

meta_notnull = globals()['meta_notnull_resort'].copy()

meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom)  
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)

meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates

meta_notnull[['chrom', 'start', 'end', 'name']].to_csv('/data/project/3dith/data/450k_metadata.opensea.bed', sep = '\t', index = False, header = False)
cmd_ = 'bedtools sort -i /data/project/3dith/data/450k_metadata.opensea.bed > /data/project/3dith/data/450k_metadata.opensea.sorted.bed && rm /data/project/3dith/data/450k_metadata.opensea.bed'
os.system(cmd_)
