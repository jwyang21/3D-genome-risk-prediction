#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))


# In[2]:


df = pd.read_csv('/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/pcc_sc_avg-beta.csv', index_col = 0)


# In[5]:


fig = plt.figure(figsize = (2, 1.5))

ax = fig.add_subplot(111)
ax.hist(df.pcc_all.values)
#ax.set_title('All samples', fontsize = 8)

ax.set_title("histogram of pcc(stem closeness, avg beta)", fontsize = 7)
fig.tight_layout()
plt.savefig('/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/pcc_sc_avg-beta-ALL_samples.png')
print('/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/pcc_sc_avg-beta-ALL_samples.png')


# In[22]:


fig = plt.figure(figsize = (3, 1.5))
for i in range(df.shape[1]):
    ax = fig.add_subplot(1, 2, i+1)
    ax.hist(df.iloc[:,i].values)
    if 'all' in df.columns[i]:
        ax.set_title('All samples', fontsize = 8)
    else:
        ax.set_title('Tumor samples', fontsize = 8)
plt.suptitle("histogram of pcc(stem closeness, avg_beta) values")
fig.tight_layout()
plt.savefig('/data/project/3dith/pipelines/opensea-pipeline/2_downstream-opensea/result/pcc_sc_avg-beta.png')


# In[17]:





# In[ ]:




