#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os


# In[2]:


meta = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = 7)


# In[3]:


raw_meta = meta.copy()


# In[70]:


meta.Relation_to_UCSC_CpG_Island.unique()


# ---

# ## Preliminary
# - 지금 있는 resort bedfile 그대로 만들어 본 다음에, 똑같이 재현되는지 확인 후 island, shore, shelf cpg들 뽑기

# - 일단 이 metadata 쓰는게 맞는지 확인하기 위해, https://www.notion.so/dohlee/450K-probe-metadata-CpG-resort-only-processed-137dd157c1034e65b5713ff17f2c82ce 의 코드 그대로 따라해 보고 노션과 같은 결과 숫자들이 나오는지 확인.

# In[10]:


cpg_type = 'resort'
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & meta.Relation_to_UCSC_CpG_Island.notnull()].copy()
print(globals()['meta_notnull_'+cpg_type].shape)
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))


# - https://www.notion.so/dohlee/450K-probe-metadata-processed-fc38c2d1c01d4de3a74fe269ea5726e1 의 코드를 그대로 따라했을 때, 결과 bed file '/data/project/3dith/data/450k_metadata.resort.sorted.bed'과 동일한 파일이 나오는지 확인

# In[30]:


meta_notnull = globals()['meta_notnull_resort'].copy()


# In[4]:


def process_chrom(c):
    if isinstance(c, float):
        return f'chr{int(c)}'
    elif isinstance(c, int):
        return f'chr{c}'
    else:
        return 'chr'+c


# In[32]:


meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom) #CHR column이 int, str 섞여있어서 정리. 
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)


# In[33]:


meta_notnull.CHR.unique()


# In[34]:


meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BEd should be 0-based, subtract 1 from 1-based probe coordinates


# In[37]:


meta_notnull[['chrom', 'start', 'end', 'name']] #원래 metadata가 깔끔하게 정리됨


# In[38]:


meta_notnull[['chrom', 'start', 'end', 'name']].to_csv('/data/project/jeewon/test_450k_resort.bed', sep = '\t', index = False, header = False)
get_ipython().system("bedtools sort -i '/data/project/jeewon/test_450k_resort.bed' > '/data/project/jeewon/test_450k_resort.sorted.bed' && rm '/data/project/jeewon/test_450k_resort.bed'")


# In[39]:


test_df = pd.read_csv('/data/project/jeewon/test_450k_resort.sorted.bed', sep = '\t', header = None)


# In[40]:


desired_df = pd.read_csv('/data/project/3dith/data/450k_metadata.resort.sorted.bed', sep = '\t', header = None)


# In[41]:


display(test_df.head(3))
print("---")
display(desired_df.head(3))


# In[42]:


print(test_df.shape)
print("---")
print(desired_df.shape)


# - 이렇게 하면 될 듯!

# ---

# ## (1) make island CpG bedfile

# ### (1-1) 먼저, metadata 내의 cpg_type들 확인

# In[6]:


print(meta.Relation_to_UCSC_CpG_Island.unique())


# ### (1-2) 현재 타겟하는 cpg_type을 specify한 후 metadata에서 그 type의 cpg probe들만 추출.

# In[29]:


cpg_type = 'Island'
#globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & meta.Relation_to_UCSC_CpG_Island.notnull()].copy() # resort ver.
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & (meta.Relation_to_UCSC_CpG_Island==cpg_type)].copy() #island ver.
 
print("meta_notnull_{} shape: {}: ".format(cpg_type, globals()['meta_notnull_'+cpg_type].shape))
print("proportion of meta_notnull_{} in whole metadata: ".format(cpg_type), end = '')
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))
meta_notnull = globals()['meta_notnull_'+cpg_type].copy()

meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom) #CHR column이 int, str 섞여있어서 정리. 
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)

meta_notnull.CHR.unique()

meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BEd should be 0-based, subtract 1 from 1-based probe coordinates

meta_notnull[['chrom', 'start', 'end', 'name']] #원래 metadata가 깔끔하게 정리됨


# ### (1-3) 결과파일을 저장하고자 하는 디렉토리 및 파일명 확인

# In[39]:


result_dir = '/data/project/3dith/data/'

result_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.bed')
print(result_fname)


# ### (1-4) 결과파일 저장

# In[40]:


meta_notnull[['chrom', 'start', 'end', 'name']].to_csv(result_fname, sep = '\t', index = False, header = False)


# ### (1-5) 결과파일을 sorting한 후 어떤 파일명으로 저장할지 확인

# In[40]:


result_sorted_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.sorted.bed')
print(result_sorted_fname)


# ### (1-6) bedtools command 확인

# In[42]:


command_ = "bedtools sort -i {} > {} && rm {}".format(result_fname, result_sorted_fname, result_fname)
print(command_)


# ### (1-7) bedtools 돌려서 결과파일을 sorting한 bedfile 만들기

# In[43]:


os.system(command_) #!bedtools sort -i '/data/project/jeewon/test_450k_resort.bed' > '/data/project/jeewon/test_450k_resort.sorted.bed' && rm '/data/project/jeewon/test_450k_resort.bed'


# ### (1-8) sorting 완료된 bedfile 불러와서, cpg probe 개수 맞는지 확인

# In[46]:


assert pd.read_csv(result_sorted_fname, header = None, sep = '\t').shape[0] == meta_notnull.shape[0]


# ## (2) make shore CpG bedfile

# In[151]:


meta = raw_meta.copy()


# ### (2-1) 먼저 각 cpg_type의 cpg probe 개수 확인

# In[152]:


print(meta.Relation_to_UCSC_CpG_Island.unique())


# In[153]:


print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# ### (2-2) shore의 라벨은 S_Shelf와 N_Shelf 두 개 존재하기 때문에, process_chrom 함수와 비슷한 process_cpg_type 함수를 만들어서 'Relation_to_UCSC_CpG_Island' column을 process하자

# In[5]:


def process_cpg_type(c):
    if str(c) in ['N_Shore', 'S_Shore']:
        return 'shore'
    elif str(c) in ['N_Shelf', 'S_Shelf']:
        return 'shelf'
    else:
        return c


# ### (2-3) 원래 meta는 raw_meta에 저장해뒀기 때문에, meta를 그대로 변형해서 사용.
# - 대신, 다음에 shelf 시작하기 전에 raw_meta로 다시 meta를 다시 초기화

# In[155]:


meta['Relation_to_UCSC_CpG_Island'] = meta.Relation_to_UCSC_CpG_Island.apply(process_cpg_type)


# ### (2-4) cpg_type이 제대로 바뀌었는지 확인

# In[156]:


print(meta.Relation_to_UCSC_CpG_Island.unique())


# In[157]:


print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# - 제대로 바뀜!

# In[159]:


cpg_type = 'shore'
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & (meta.Relation_to_UCSC_CpG_Island==cpg_type)].copy() 


# In[160]:


print("meta_notnull_{} shape: {}: ".format(cpg_type, globals()['meta_notnull_'+cpg_type].shape))
print("proportion of meta_notnull_{} in whole metadata: ".format(cpg_type), end = '')
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))
meta_notnull = globals()['meta_notnull_'+cpg_type].copy()


# In[161]:


meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom) #CHR column이 int, str 섞여있어서 정리. 
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)


# In[162]:


meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BEd should be 0-based, subtract 1 from 1-based probe coordinates


# In[163]:


result_dir = '/data/project/3dith/data/'

result_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.bed')
print(result_fname)


# In[164]:


meta_notnull[['chrom', 'start', 'end', 'name']].to_csv(result_fname, sep = '\t', index = False, header = False)


# In[165]:


result_sorted_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.sorted.bed')
print(result_sorted_fname)


# In[166]:


command_ = "bedtools sort -i {} > {} && rm {}".format(result_fname, result_sorted_fname, result_fname)
print(command_)


# In[167]:


os.system(command_) 


# In[168]:


assert pd.read_csv(result_sorted_fname, header = None, sep = '\t').shape[0] ==  meta_notnull.shape[0]


# ## (3) make shelf CpG bedfile

# In[6]:


meta = raw_meta.copy()

print(meta.Relation_to_UCSC_CpG_Island.unique())

print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# In[7]:


meta['Relation_to_UCSC_CpG_Island'] = meta.Relation_to_UCSC_CpG_Island.apply(process_cpg_type)

print(meta.Relation_to_UCSC_CpG_Island.unique())

print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# In[8]:


cpg_type = 'shelf'
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & (meta.Relation_to_UCSC_CpG_Island==cpg_type)].copy() 

print("meta_notnull_{} shape: {}: ".format(cpg_type, globals()['meta_notnull_'+cpg_type].shape))
print("proportion of meta_notnull_{} in whole metadata: ".format(cpg_type), end = '')
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))
meta_notnull = globals()['meta_notnull_'+cpg_type].copy()


# In[9]:


meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom) #CHR column이 int, str 섞여있어서 정리. 
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)

meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BEd should be 0-based, subtract 1 from 1-based probe coordinates


# In[10]:


result_dir = '/data/project/3dith/data/'

result_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.bed')
print(result_fname)


# In[11]:


meta_notnull[['chrom', 'start', 'end', 'name']].to_csv(result_fname, sep = '\t', index = False, header = False)


# In[12]:


result_sorted_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.sorted.bed')
print(result_sorted_fname)


# In[13]:


command_ = "bedtools sort -i {} > {} && rm {}".format(result_fname, result_sorted_fname, result_fname)
print(command_)


# In[14]:


os.system(command_) 


# In[15]:


assert pd.read_csv(result_sorted_fname, header = None, sep = '\t').shape[0] ==  meta_notnull.shape[0]


# ## (3) make shelf_shore CpG bedfile

# In[16]:


meta = raw_meta.copy()

print(meta.Relation_to_UCSC_CpG_Island.unique())

print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# In[17]:


def process_cpg_type_shelf_shore(c):
    if str(c) in ['N_Shore', 'S_Shore', 'N_Shelf', 'S_Shelf']:
        return 'shelf_shore'
    else:
        return c


# In[18]:


meta['Relation_to_UCSC_CpG_Island'] = meta.Relation_to_UCSC_CpG_Island.apply(process_cpg_type_shelf_shore)

print(meta.Relation_to_UCSC_CpG_Island.unique())

print(meta.Relation_to_UCSC_CpG_Island.value_counts())


# In[19]:


cpg_type = 'shelf_shore'
globals()['meta_notnull_'+cpg_type] = meta[(meta.Genome_Build == 37.0) & meta.CHR.notnull() & meta.MAPINFO.notnull() & (meta.Relation_to_UCSC_CpG_Island==cpg_type)].copy() 

print("meta_notnull_{} shape: {}: ".format(cpg_type, globals()['meta_notnull_'+cpg_type].shape))
print("proportion of meta_notnull_{} in whole metadata: ".format(cpg_type), end = '')
print(len(globals()['meta_notnull_'+cpg_type])/len(meta))
meta_notnull = globals()['meta_notnull_'+cpg_type].copy()


# In[20]:


meta_notnull['CHR'] = meta_notnull.CHR.apply(process_chrom) #CHR column이 int, str 섞여있어서 정리. 
meta_notnull['MAPINFO'] = meta_notnull.MAPINFO.astype(int)

meta_notnull = meta_notnull.rename({'IlmnID':'name', 'MAPINFO':'start', 'CHR':'chrom'}, axis = 1)
meta_notnull['end'] = meta_notnull['start'] +2 -1 # since BED should be 0-based, subtract 1 from 1-based probe coordinates
meta_notnull['start'] = meta_notnull['start'] -1 # since BEd should be 0-based, subtract 1 from 1-based probe coordinates


# In[21]:


result_dir = '/data/project/3dith/data/'

result_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.bed')
print(result_fname)


# In[22]:


meta_notnull[['chrom', 'start', 'end', 'name']].to_csv(result_fname, sep = '\t', index = False, header = False)


# In[23]:


result_sorted_fname = os.path.join(result_dir, '450k_metadata.'+cpg_type.lower()+'.sorted.bed')
print(result_sorted_fname)


# In[24]:


command_ = "bedtools sort -i {} > {} && rm {}".format(result_fname, result_sorted_fname, result_fname)
print(command_)


# In[25]:


os.system(command_) 


# In[26]:


assert pd.read_csv(result_sorted_fname, header = None, sep = '\t').shape[0] ==  meta_notnull.shape[0]


# In[ ]:




