#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
import os


# In[5]:


log_dir = '/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/log'


# In[6]:


#log_files = ['6_score_meth_corr-score2.log', '6_score_meth_corr-score3.log','6_score_meth_corr-score4.log']
log_files = ['6_score_meth_corr_score7.log']


# In[7]:


#cohort:
#PCC:
#pvalue:


# In[8]:


'''
# readline_all.py
f = open("C:/doit/새파일.txt", 'r')
while True:
    line = f.readline()
    if not line: break
    print(line)
f.close()
'''


# In[9]:


os.chdir('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses')


# In[10]:


SAVEDIR = os.path.join(os.getcwd(), 'result')


# In[11]:


'''
#test
l = 'cohort: TCGA-LIHC'
l2 = 'PCC: -0.8789970389101859, pvalue: 2.0061069605680253e-139'
l.split(':')[-1].strip() #cohort
float(l2.split(':')[1].split(',')[0].strip()) #pcc
float(l2.split('pvalue:')[-1].strip()) #pvalue

test = '6_score_meth_corr-score2.log'
test.split('.')[0].split('-')[-1]
'''


# In[12]:




# In[13]:


for f in log_files:
    #score = f.split('.')[0].split('-')[-1]
    score = f.split('_')[1].split('.')[0]
    cohorts = []
    pcc = []
    p = []
    #read_parse
    f_ = open(os.path.join(log_dir, f), 'r')   
    while True:
        line = f_.readline()
        if not line: break
        if line.startswith('cohort:'):
            cohorts.append(line.split(':')[-1].strip())
        elif line.startswith('PCC:'):
            pcc.append(float(line.split(':')[1].split(',')[0].strip()))
            p.append(float(line.split('pvalue:')[-1].strip()))
    f_.close()
    globals()[score] = pd.DataFrame(zip(pcc, p), index = cohorts, columns = ['pcc', 'p-value'])


# In[14]:


print(score)
print(cohorts)
print(pcc)
print(p)


# In[ ]:


'''
# for score2, score3 and score4
#df_name = 'opensea_avg_beta-'+score+'.csv'
for s in ['score2', 'score3', 'score4']:
    df_name = 'opensea_avg_beta-'+s+'.csv'
    print("df name: {}".format(os.path.join(SAVEDIR, df_name)))
    globals()[s].to_csv(os.path.join(SAVEDIR, df_name), index = True)
'''
#'opensea_avg_beta-'+'score2'+'.csv'


# In[15]:


df = pd.DataFrame(zip(pcc, p), index = cohorts, columns = ['pcc', 'p-value'])


# In[13]:


# for score5


# In[16]:


# for score7
df.head(3)
score = 'score7'
fname = 'opensea_avg_beta-'+score+'.csv'
print(fname)


# In[17]:


df.to_csv(os.path.join(SAVEDIR, fname))
print(os.path.join(SAVEDIR, fname))
#pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/opensea_avg_beta-score5.csv', index_col=0)


# In[ ]:


# 각 TCGA cohort의 score과 beta value 간의 pearson correlation.


# In[19]:
'''
# check the file by re-loading.
pd.read_csv('/data/project/jeewon/research/3D-ITH/pipelines/downstream-analyses/result/opensea_avg_beta-score7.csv', index_col = 0)
'''




