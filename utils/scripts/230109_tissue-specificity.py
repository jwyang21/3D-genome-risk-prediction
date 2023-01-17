#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os
import glob
from collections import Counter, defaultdict
from tqdm import tqdm
from scipy import stats
from scipy.stats import pearsonr
from tqdm.auto import tqdm
from time import sleep

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))
pd.set_option("display.max_columns", None)
def save_figures(f, exts=['png', 'pdf']):
    for ext in exts:
        plt.savefig(f + f'.{ext}', dpi=300, bbox_inches='tight', transparent=True)



NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[3]:


cpg_type = 'opensea'
assert cpg_type in ['opensea', 'island', 'shelf_shore']


# In[4]:


CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


# In[5]:


# fix random seed
np.random.seed(42)


# ## functions

# In[6]:


def smoothen(v, window):
    return pd.Series(v).rolling(window=window, center=True).agg('mean').dropna().values

def standardize(v):
    return (v - v.mean()) / v.std()


# In[7]:


tcga_fire_cohort = pd.read_csv('/data/project/3dith/data/etc/tcga-fire-cohorts.csv', index_col = 0)
# tcga_fire_cohort.loc['TCGA-LIHC','FIRE'] # -> 'LI'
# print(tcga_fire_cohort.columns) #Index(['FIRE', 'TCGA_disease', 'FIRE_tissue'], dtype='object')


# In[8]:


tcga_fire_cohort


# In[9]:


tcga_cohort_info = pd.read_csv('/data/project/3dith/data/etc/manifest.csv', index_col = 0)


# In[10]:


tcga_cohort_info.head(3)


# ## 1. tissue specificity



# In[12]:


save_dir = f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors'

cohort_few_normal = ['TCGA-CHOL', 'TCGA-ESCA', 'TCGA-PAAD']


# 각 cohort마다, tumor/normal group 각각을 random shuffle해서 index list를 얻은 후, 맨 앞에서부터 10개씩 잘라서 tumor/normal representative vector를 만든다. 

# In[13]:

'''
for chrom in CHR_LIST:
    print(f"==\nchrom: {chrom}")
    for cohort in NORMAL7_COHORT:
        files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*_inv_exp.npz') 
        if chrom == CHR_LIST[0]:
            #print(f"{len(files)} samples in {cohort}")
            
            if not os.path.exists(os.path.join(save_dir, cohort)):
                os.makedirs(os.path.join(save_dir, cohort))
            
        # make representative vector
        data_T = []
        data_N = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            #print(int(sample[13:15]) <= 9)
            if int(sample[13:15])<=9:
                data_T.append(standardize(np.load(f)[f'{chrom}_pc1']))
            else:
                data_N.append(standardize(np.load(f)[f'{chrom}_pc1']))

            #print(len(data_T), len(data_N))
        data_T = np.vstack(data_T)
        data_N = np.vstack(data_N)

        idx_T = np.arange(len(data_T))
        np.random.shuffle(idx_T)

        idx_N = np.arange(len(data_N))
        np.random.shuffle(idx_N)

        data_T = data_T[idx_T]
        data_N = data_N[idx_N]

        num_repr_T_vectors = (len(data_T) // 10)
        num_repr_N_vectors = (len(data_N) // 10)

        for i in range(num_repr_T_vectors):
            current_T_repr = data_T[(i) * 10 : (i+1)*10]
            assert len(current_T_repr) == 10
            np.save(f'{save_dir}/{cohort}/repr_T_vector_{chrom}_{str(i)}.npy', current_T_repr)
            #print(f'{save_dir}/{cohort}/repr_T_vector_{chrom}_{str(i)}.npy')

        for i in range(num_repr_N_vectors):
            if cohort not in cohort_few_normal:
                current_N_repr = data_N[(i) * 10 : (i+1)*10]
                assert len(current_N_repr) == 10
                np.save(f'{save_dir}/{cohort}/repr_N_vector_{chrom}_{str(i)}.npy', current_T_repr)
                #print(f'{save_dir}/{cohort}/repr_N_vector_{chrom}_{str(i)}.npy')

'''
# In[41]:





# In[38]:



# In[14]:



#------
### bhi4 tmux1 (230110)
# 모든 cohort pair 간 T - T 비교

# initialize
for target_cohort in NORMAL7_COHORT:
    target_type = target_cohort.split('-')[1]
    for probe_cohort in NORMAL7_COHORT:
        probe_type = probe_cohort.split('-')[1]
        globals()[f'{target_type}_{probe_type}'] = {}
        for chrom in CHR_LIST:
            globals()[f'{target_type}_{probe_type}'][f'{chrom}'] = []


# In[ ]:

print("Tumor-Tumor comparison between all cohort pair") 

type_ = 'Tumor'
type_flag = type_[0]

pbar_chrom = tqdm(CHR_LIST, total = len(CHR_LIST), desc = '===\nchrom', ncols = 100, ascii = ' =', leave = True)#real
#pbar_chrom = tqdm(CHR_LIST[:1], total = len(CHR_LIST), desc = '===\nchrom', ncols = 100, ascii = ' =', leave = True)#debug
#for chrom in tqdm(CHR_LIST):

for chrom in pbar_chrom:
#for chrom in tqdm(chrom_barcodes, desc = 'chrom')
    pbar_chrom.set_description(f'{chrom}')
    #sleep(10)
    #print(f"===\nchrom: {chrom}")
    pbar_target = tqdm(NORMAL7_COHORT, total = len(NORMAL7_COHORT), desc = 'target cohort', ncols = 100, ascii = ' =', leave = True)#real
    #pbar_target = tqdm(NORMAL7_COHORT[:2], total = len(NORMAL7_COHORT), desc = 'target cohort', ncols = 100, ascii = ' =', leave = True)#debug
    for target_cohort in pbar_target:
    #for target_cohort in tqdm(NORMAL7_COHORT):
        sleep(1000)
        pbar_target.set_description(f'---\ntarget: {target_cohort}')
        target_type = target_cohort.split('-')[1]
        target_files = glob.glob(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/{target_cohort}/repr_{type_flag}_vector_{chrom}_*.npy')
        target_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{target_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
        #print(f"---\ntarget cohort: {target_cohort}")
        for target_f in tqdm(target_files, desc = 'target files'):#real
        #for target_f in tqdm(target_files[:2], desc = 'target files'):#debug
            current_target = np.load(target_f).mean(axis = 0)
            assert len(current_target) == len(target_bins)
            
            pbar_probe = tqdm(NORMAL7_COHORT, total = len(NORMAL7_COHORT), desc = 'probe cohort', ncols = 100, ascii = ' =', leave = True)#real
            #pbar_probe = tqdm(NORMAL7_COHORT[NORMAL7_COHORT.index(target_cohort):2], total = len(NORMAL7_COHORT), desc = 'probe cohort', ncols = 100, ascii = ' =', leave = True)#debug
            for probe_cohort in pbar_probe:
            #for probe_cohort in tqdm(NORMAL7_COHORT):
                sleep(1000)
                pbar_probe.set_description(f'probe: {probe_cohort}')
                #print(f"probe cohort: {probe_cohort}")
                probe_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{probe_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
                probe_type = probe_cohort.split('-')[1]
                probe_files = glob.glob(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/{probe_cohort}/repr_{type_flag}_vector_{chrom}_*.npy')
                for probe_f in tqdm(probe_files, desc = 'probe files'):#real
                #for probe_f in tqdm(probe_files[:2], desc = 'probe files'):#debug
                    current_probe = np.load(probe_f).mean(axis = 0)
                    assert len(current_probe) == len(probe_bins)
                    
                    if target_type != probe_type:                        
                        intersecting_bins = np.intersect1d(target_bins, probe_bins)
                        
                        target_mask = [x in intersecting_bins for x in target_bins]
                        probe_mask = [x in intersecting_bins for x in probe_bins]
                        
                        masked_target = current_target[target_mask] 
                        masked_probe = current_probe[probe_mask]
                    else:
                        masked_target = current_target.copy()
                        masked_probe = current_probe.copy()
                        
                        
                    assert len(masked_target) == len(masked_probe)  
                    #print(pearsonr(masked_target, masked_probe)[0])
                    #print(f'{target_type}_{probe_type}')
                    globals()[f'{target_type}_{probe_type}'][f'{chrom}'].append(pearsonr(masked_target, masked_probe)[0])


# In[ ]:


# start from here.
# Tumor representative vector간의 pcc 계산. 
# 각 chromosome 별로, 모든 가능한 (target cohort, probe cohort) pair에 대해 '이 pair의 이 chrom에서의 target cohort 의 대표벡터와 probe cohort의 대표벡터 간 pcc들의 median 값'을 계산 
## -> target cohort와 probe cohort 간 대표 pcc
for chrom in CHR_LIST:
    globals()[f'all_df_{chrom}'] = pd.DataFrame(np.zeros((len(NORMAL7_COHORT), len(NORMAL7_COHORT)), dtype = float), index = NORMAL7_COHORT, columns = NORMAL7_COHORT)
    for target_cohort in NORMAL7_COHORT:
        target_type = target_cohort.split('-')[1]
        for probe_cohort in NORMAL7_COHORT[NORMAL7_COHORT.index(target_cohort):]:
            probe_type = probe_cohort.split('-')[1]
            globals()[f'all_df_{chrom}'].loc[target_cohort][probe_cohort] = np.median(globals()[f'{target_type}_{probe_type}'][f'{chrom}'])
            globals()[f'all_df_{chrom}'].loc[probe_cohort][target_cohort] = np.median(globals()[f'{target_type}_{probe_type}'][f'{chrom}'])


# In[ ]:


for chrom in CHR_LIST:
    globals()[f'all_df_{chrom}'].to_csv(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/pcc/{type_.lower()}-{type_.lower()}_{chrom}.csv')
    del(globals()[f'all_df_{chrom}'])


# ---
### bhi4 tmux2 (230110)
# 모든 cohort pair 간 N - N 비교
# In[ ]:


NORMAL7_COHORT_w_N = list(set(NORMAL7_COHORT) - set(cohort_few_normal))
del(NORMAL7_COHORT)


# In[ ]:


# initialize
for target_cohort in NORMAL7_COHORT_w_N:
    target_type = target_cohort.split('-')[1]
    for probe_cohort in NORMAL7_COHORT_w_N:
        probe_type = probe_cohort.split('-')[1]
        globals()[f'{target_type}_{probe_type}'] = {}
        for chrom in CHR_LIST:
            globals()[f'{target_type}_{probe_type}'][f'{chrom}'] = []


# In[ ]:



type_ = 'Normal'
type_flag = type_[0]
print(f"type: {type_}, type_flag: {type_flag}")

pbar_chrom = tqdm(CHR_LIST, total = len(CHR_LIST), desc = '===\nchrom', ncols = 100, ascii = ' =', leave = True)#real
#pbar_chrom = tqdm(CHR_LIST[:1], total = len(CHR_LIST), desc = '===\nchrom', ncols = 100, ascii = ' =', leave = True)#debug
#for chrom in tqdm(CHR_LIST):


for chrom in pbar_chrom:
#for chrom in tqdm(chrom_barcodes, desc = 'chrom')
    pbar_chrom.set_description(f'{chrom}')
    #sleep(10)
    #print(f"===\nchrom: {chrom}")
    pbar_target = tqdm(NORMAL7_COHORT_w_N, total = len(NORMAL7_COHORT_w_N), desc = 'target cohort', ncols = 100, ascii = ' =', leave = True)#real
    #pbar_target = tqdm(NORMAL7_COHORT_w_N[:2], total = len(NORMAL7_COHORT_w_N), desc = 'target cohort', ncols = 100, ascii = ' =', leave = True)#debug
    for target_cohort in pbar_target:
    #for target_cohort in tqdm(NORMAL7_COHORT_w_N):
        sleep(1000)
        pbar_target.set_description(f'---\ntarget: {target_cohort}')
        target_type = target_cohort.split('-')[1]
        target_files = glob.glob(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/{target_cohort}/repr_{type_flag}_vector_{chrom}_*.npy')
        target_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{target_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
        #print(f"---\ntarget cohort: {target_cohort}")
        for target_f in tqdm(target_files, desc = 'target files'):#real
        #for target_f in tqdm(target_files[:2], desc = 'target files'):#debug
            current_target = np.load(target_f).mean(axis = 0)
            assert len(current_target) == len(target_bins)
            
            pbar_probe = tqdm(NORMAL7_COHORT_w_N, total = len(NORMAL7_COHORT_w_N), desc = 'probe cohort', ncols = 100, ascii = ' =', leave = True)#real
            #pbar_probe = tqdm(NORMAL7_COHORT_w_N[NORMAL7_COHORT_w_N.index(target_cohort):2], total = len(NORMAL7_COHORT_w_N), desc = 'probe cohort', ncols = 100, ascii = ' =', leave = True)#debug
            for probe_cohort in pbar_probe:
            #for probe_cohort in tqdm(NORMAL7_COHORT_w_N):
                sleep(1000)
                pbar_probe.set_description(f'probe: {probe_cohort}')
                #print(f"probe cohort: {probe_cohort}")
                probe_bins = np.load(f'/data/project/3dith/data/bdm_bins/{cpg_type}/{probe_cohort}_diffmat_bins.npz')[f'{chrom}_bins']
                probe_type = probe_cohort.split('-')[1]
                probe_files = glob.glob(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/{probe_cohort}/repr_{type_flag}_vector_{chrom}_*.npy')
                for probe_f in tqdm(probe_files, desc = 'probe files'):#real
                #for probe_f in tqdm(probe_files[:2], desc = 'probe files'):#debug
                    current_probe = np.load(probe_f).mean(axis = 0)
                    assert len(current_probe) == len(probe_bins)
                    
                    if target_type != probe_type:                        
                        intersecting_bins = np.intersect1d(target_bins, probe_bins)
                        
                        target_mask = [x in intersecting_bins for x in target_bins]
                        probe_mask = [x in intersecting_bins for x in probe_bins]
                        
                        masked_target = current_target[target_mask] 
                        masked_probe = current_probe[probe_mask]
                    else:
                        masked_target = current_target.copy()
                        masked_probe = current_probe.copy()
                        
                        
                    assert len(masked_target) == len(masked_probe)  
                    #print(pearsonr(masked_target, masked_probe)[0])
                    #print(f'{target_type}_{probe_type}')
                    globals()[f'{target_type}_{probe_type}'][f'{chrom}'].append(pearsonr(masked_target, masked_probe)[0])


# In[ ]:


# Tumor representative vector간의 pcc 계산. 
# 각 chromosome 별로, 모든 가능한 (target cohort, probe cohort) pair에 대해 '이 pair의 이 chrom에서의 target cohort 의 대표벡터와 probe cohort의 대표벡터 간 pcc들의 median 값'을 계산 
## -> target cohort와 probe cohort 간 대표 pcc
for chrom in CHR_LIST:
    globals()[f'all_df_{chrom}'] = pd.DataFrame(np.zeros((len(NORMAL7_COHORT_w_N), len(NORMAL7_COHORT_w_N)), dtype = float), index = NORMAL7_COHORT_w_N, columns = NORMAL7_COHORT_w_N)
    for target_cohort in NORMAL7_COHORT_w_N:
        target_type = target_cohort.split('-')[1]
        for probe_cohort in NORMAL7_COHORT_w_N[NORMAL7_COHORT_w_N.index(target_cohort):]:
            probe_type = probe_cohort.split('-')[1]
            globals()[f'all_df_{chrom}'].loc[target_cohort][probe_cohort] = np.median(globals()[f'{target_type}_{probe_type}'][f'{chrom}'])
            globals()[f'all_df_{chrom}'].loc[probe_cohort][target_cohort] = np.median(globals()[f'{target_type}_{probe_type}'][f'{chrom}'])


# In[ ]:


for chrom in CHR_LIST:
    globals()[f'all_df_{chrom}'].to_csv(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/2_downstream-{cpg_type}/result/repr_vectors/pcc/{type_.lower()}-{type_.lower()}_{chrom}.csv')
    del(globals()[f'all_df_{chrom}'])


# In[24]:


BRCA_LUAD['chr1']


# In[35]:


target_cohort


# In[36]:


probe_cohort


# In[37]:


chrom


# In[34]:


len(target_bins)


# In[ ]:


current_tar


# In[30]:


current_target.shape


# In[31]:


current_probe.shape


# In[33]:


len(target_mask)


# In[21]:


target_cohort


# In[23]:


probe_cohort


# In[25]:


len(target_mask)


# In[27]:


len(probe_mask)


# In[ ]:





# In[ ]:


# 동일 cohort 내에서 T - T 비교


# In[ ]:


# 서로 다른 cohort끼리 Normal representative vector 비교


# In[ ]:


# 동일 cohort 내에서 Tumor representative vector - Normal representative vector 비교


# In[ ]:


arr = np.arange(10)
np.random.shuffle(arr)
arr


# In[ ]:


data_T


# In[ ]:


data_T[idx_T]


# In[ ]:


idxs = np.random.shuffle(np.arange(len(data_T)))


# In[ ]:


idxs


# In[ ]:


for chrom in CHR_LIST:
    print("<{}>".format(chrom))
    for i in range(len(NORMAL7_COHORT)):


        cohort = NORMAL7_COHORT[i]
        print("===\ncohort: {}".format(cohort))

        #cohort = 'TCGA-LUSC'
        #cohort = 'TCGA-LUAD'
        #files = glob.glob(f'/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/{cohort}/*_inv_exp.npz') #이전 버전
        #-------------
        # 1. Prepare data before plotting.
        files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*_inv_exp.npz') #새로 계산한 버전
        tmp = np.load(files[11])
        len(files)
        print("1. collect 10 tumor samples' pc1 vectors")
        data = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            if sample.endswith('11'):
                continue

            #data.append(np.load(f)[chrom])#이전 버전
            data.append(np.load(f)[chrom+'_pc1'])#현재 버전

        data = np.vstack(data)

        np.random.seed(42)
        idxs = np.random.choice(range(len(data)), size=10)
        data = data[idxs]

        np.save('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy', data)
        print('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')

        print("---\n2. collect 10 normal samples' data.")
        data = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            if not sample.endswith('11'):#tcga barcode가 11로 끝나지 않으면. 즉 tumor sample이면 -> 건너뜀.
                continue

            #data.append(np.load(f)[chrom])#이전 버전
            data.append(np.load(f)[chrom+'_pc1'])

        data = np.vstack(data)

        np.random.seed(42)
        idxs = np.random.choice(range(len(data)), size=10)
        data = data[idxs]

        np.save('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy', data)
        print('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')
        #------------------------------
        # 2. plot 10 tumor samples' PC1 vectors
        # Plot PC1 vectors from 10 randomly picked tumor samples of current cohort. 
        # draw 10 tumor samples
        fig = plt.figure(figsize = (4, 0.3*22))

        data = np.load('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')
        #ax = fig.add_subplot(1, 2, 1)
        #fig = plt.figure(figsize=(4, 3))

        for i in range(10):
            ax = fig.add_subplot(22, 1, i+1)
            ax.spines['left'].set_linewidth(0.6)#선 두께 #https://zephyrus1111.tistory.com/241
            if i == 0:
                ax.set_title("{}".format(cohort))
            if i == 4:
                ax.set_ylabel('Tumor PC1', fontsize = 6.5)

            y = smoothen(standardize(data[i]), window=3)
            x = np.arange(len(y))

            ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
            ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)

            for direction in ['top', 'right', 'bottom']:
                ax.spines[direction].set_visible(False)
            ax.spines['left'].set_position(('outward', 3))

            ax.set_xlim([0, len(y)])
            ax.set_xticks([])
            ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)

            ax.tick_params(axis='y', labelsize=4.5, length=2, pad=1, width = 0.6)
            '''
            if i == 9:
                ax.set_xlabel('Position ({})'.format(chrom))

                mb_length = 1 / len(y)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.6, -0.6], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.6, '10Mb', transform=ax.transAxes, ha='left', va='center')
            '''
        #fig.supylabel('PC1', fontsize=7)

        #save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')
        #print('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')

        #------------------------------
        # 3. plot 10 normal samples' PC1 vectors. 
        # Plot PC1 vectors from 10 randomly picked normal samples of current cohort. 
        data = np.load('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')

        #fig = plt.figure(figsize=(4, 3))

        for i in range(10):
            ax = fig.add_subplot(22, 1, i+12)
            ax.spines['left'].set_linewidth(0.6)
            #if i==0:
                #ax.set_title('Normal')
            if i == 4:
                ax.set_ylabel('Normal PC1', fontsize = 6.5)

            y = smoothen(standardize(data[i]), window=3)
            x = np.arange(len(y))

            ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
            ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)

            for direction in ['top', 'right', 'bottom']:
                ax.spines[direction].set_visible(False)
            ax.spines['left'].set_position(('outward', 3))

            ax.set_xlim([0, len(y)])
            ax.set_xticks([])
            ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)

            ax.tick_params(axis='y', labelsize=4.5, length=2, pad=1, width = 0.6)

            if i == 9:
                ax.set_xlabel('Position ({})'.format(chrom))

                mb_length = 1 / len(y)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.8, -0.8], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.8, '1Mb', transform=ax.transAxes, ha='left', va='center')

        #fig.supylabel('PC1', fontsize=7)
        #plt.suptitle('{}'.format(cohort), va='center')
        fig.tight_layout()
        #



        #save_figures('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')
        #print('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')
        save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_normal_10samples_pc1_'+chrom)
        print('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_normal_10samples_pc1_'+chrom)
        #plt.clf()


# In[ ]:


# select target cohort
cohort = 'TCGA-LIHC'
#cohort = 'TCGA-BLCA'
assert cohort in NORMAL7_COHORT


# In[ ]:


chrom = 'chr2'
assert chrom in CHR_LIST


# In[ ]:


#tissue_type = 'kidney'#previous version
tissue_type = tcga_cohort_info.loc[cohort, 'tissue']
print("Make sure that tissue type {} matches the selected cohort {}!!".format(tissue_type, cohort))


# In[ ]:


files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-{diffmat_cpg_type}/result/{cohort}/pc1/*_inv_exp.npz') #새로 계산한 버전
# files = glob.glob(f'/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/{cohort}/*_inv_exp.npz') #기존 버전. 
tmp = np.load(files[11]) #11번째 샘플의 pc1 데이터. 
print('number of samples in {}: {}'.format(cohort, len(files))) #이 코호트에 있는 전체 샘플 수.

data = []# tumor samples' data
data_n = [] #normal samples' data
for f in files:
    sample = os.path.basename(f).split('_')[0]
    if sample.endswith('11'): #barcode가 11로 끝나는 normal sample이면 건너뜀.
        #continue
        data_n.append(np.load(f)[chrom+'_pc1'])
        
    else:
        #data.append(np.load(f)['chr2']) #tumor sample의 chr2 데이터만 수집. #기존 버전
        data.append(np.load(f)[chrom+'_pc1']) #tumor sample의 chr2 데이터만 수집 #새로 계산한 버전. 
    
data = np.vstack(data)
data_n = np.vstack(data_n)
print("pc1 from {} tumor samples of {} are vertically stacked into data".format(data.shape[0], cohort))
print("pc1 from {} normal samples of {} are vertically stacked into data_n".format(data_n.shape[0], cohort))


# In[ ]:


fig = plt.figure(figsize=(10, 1))
ax = fig.add_subplot(111)

print("Plot pc1 of first 20 tumor samples in {} ({})".format(cohort, chrom))
for i in range(20):
    ax.plot(data[i], lw=1, alpha=0.1)


# In[ ]:


#선택된 cohort의 tumor sample의 chr2 pc1 벡터들에서, 각 genomic bin 당 median 값만 고른 것. 
print("Plot median pc1 from first 20 tumor samples in {} ({})".format(cohort, chrom))
fig = plt.figure(figsize=(10, 1))
ax = fig.add_subplot(111)

ax.plot(np.median(data, axis=0))


# In[ ]:


'''
# import fire pc1 (computed from hi-c data) and compare this with binnd diffmat-derived pc1. 
## previous version
diffmat_cpg_type = 'opensea-binned-diffmat-bins'

bins = np.load(f'/data/project/jeewon/research/3D-ITH/pipelines/find-diffmat-bins/{diffmat_cpg_type}/{cohort}_diffmat_bins.npz')
fire_pc1 = pd.read_csv('/data/project/3dith/data/fire_pc1.csv')

fire_pc1['chr'] = 'chr' + fire_pc1.chr.astype(str)
fire_pc1 = fire_pc1.rename({'chr': 'chrom'}, axis=1)
'''


# ## 1. Compare hi-c pc1 (fire pc1) and iebdm pc1.

# In[ ]:


# import fire pc1 (computed from hi-c data) and compare this with binnd diffmat-derived pc1.
fire_pc1 = pd.read_csv('/data/project/3dith/data/fire_pc1.csv')
fire_pc1['chr'] = 'chr' + fire_pc1.chr.astype(str)
fire_pc1 = fire_pc1.rename({'chr': 'chrom'}, axis=1)


# In[ ]:


# import bins used in iebdm construction of current cohort
bins = np.load(f'/data/project/3dith/data/bdm_bins/{diffmat_cpg_type}/{cohort}_diffmat_bins.npz')


# In[ ]:


# 전체 fire data 중, 현재 선택된 cohort의 tissue type과 매칭되는 fire pc1 데이터만 고르기
if cohort in tcga_fire_cohort.index.values:
    fire_type = tcga_fire_cohort.loc[cohort].values[0]
    
    # tumor data
    d = defaultdict(list)
    for bin_id, pc1 in zip(bins[f'{chrom}_bins'], np.mean(data, axis=0)):
        d['bin_id'].append(bin_id)
        d['bmdm_pc1'].append(pc1) #{bin_id}번째에 해당하는 bin에 대해, 현재 cohort의 모든 tumor sample들의 PC1 vector들을 모은 뒤 이 bin에 corresponding하는 entry의 값들만 모아서 평균낸 값 (np.mean)
    d = pd.DataFrame(d)

    d['chrom'] = d.bin_id.str.split(':', expand=True)[0]
    d['start'] = d.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int) + 1
    d['end'] = d.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    d = d[['chrom', 'start', 'end', 'bmdm_pc1']]
    d = fire_pc1.merge(d, on=['chrom', 'start', 'end'])

    d = d[['bmdm_pc1', fire_type]]


    idxs = np.random.choice(len(data), size=10)
    for i, idx in enumerate(idxs):
        d[f'bmdm_sample_{i}'] = data[idx]

    d = d.dropna()
    d.to_csv('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_tumor.csv', index=False)  
    print('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_tumor.csv')
    
    # normal data
    d_n = defaultdict(list)
    for bin_id, pc1 in zip(bins[f'{chrom}_bins'], np.mean(data_n, axis=0)):
        d_n['bin_id'].append(bin_id)
        d_n['bmdm_pc1'].append(pc1) #{bin_id}번째에 해당하는 bin에 대해, 현재 cohort의 모든 normal sample들의 PC1 vector들을 모은 뒤 이 bin에 corresponding하는 entry의 값들만 모아서 평균낸 값 (np.mean)
    d_n = pd.DataFrame(d_n)

    d_n['chrom'] = d_n.bin_id.str.split(':', expand=True)[0]
    d_n['start'] = d_n.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int) + 1
    d_n['end'] = d_n.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    d_n = d_n[['chrom', 'start', 'end', 'bmdm_pc1']]
    d_n = fire_pc1.merge(d_n, on=['chrom', 'start', 'end'])

    d_n = d_n[['bmdm_pc1', fire_type]]


    idxs = np.random.choice(len(data_n), size=10)
    for i, idx in enumerate(idxs):
        d_n[f'bmdm_sample_{i}'] = data_n[idx]

    d_n = d_n.dropna()
    d_n.to_csv('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_normal.csv', index=False)
    print('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_normal.csv')   
    
    print("tumor data shape: ", d.shape)
    print("normal data shpae: ", d_n.shape)
    
    
else:
    print("{}(tissue type: {}) has no matching type in pc1 data.".format(cohort, tissue_type))


# In[ ]:


print("Compare FIRE pc1 and iebdm-tumor-pc1")
tn_type = 'tumor'
if cohort in tcga_fire_cohort.index.values:
    fire_type = tcga_fire_cohort.loc[cohort].values[0]
    d = pd.read_csv('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_'+tn_type+'.csv')

    fig = plt.figure(figsize=(6, 0.75))
    ax = fig.add_subplot(111)

    ax.axhline(0, lw=0.75, ls='--', c='0.8')

    y1 = (d[fire_type] - d[fire_type].mean()) / d[fire_type].std() #fire pc1
    ax.plot(smoothen(y1, window=3), lw=1, c='0.5', label='Hi-C ({})'.format(tissue_type))

    y2 = d['bmdm_pc1'] 
    y2 = -(y2 - y2.mean()) / y2.std() #standardized iebdm-derived pc1
    ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimations ({}, averaged)'.format(cohort))

    for i in range(10):
        y3 = d[f'bmdm_sample_{i}']
        y3 = -(y3 - y3.mean()) / y3.std()
        ax.plot(smoothen(y3, window=3), c='C3', lw=0.75, alpha=0.2)

    for direction in ['top', 'right', 'bottom']:
        ax.spines[direction].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    
    assert len(y1)==len(y2) and len(y2)==len(y3)
    ax.set_xlim([0, len(y1)])
    ax.set_ylim([-2.5, 2.5])

    ax.set_xlabel('Position ({})'.format(chrom))
    ax.set_ylabel('PC1\n(Standardized)')

    ax.set_xticks([])
    mb_length = 1 / len(y1)
    ax.plot([0.75, 0.75 + 10*mb_length], [-0.1, -0.1], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
    ax.text(0.75 + 10*mb_length + 0.01, -0.1, '10Mb', transform=ax.transAxes, ha='left', va='center')

    ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.4))

    #save_figures('source_data_figures/'+cohort.split('-')[-1]+'_lung_pc1_comparison')
    #save_figures('source_data_figures/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1')
    print('source_data_figures/'+cohort.split('-')[-1]+'_'+tissue_type+'_'+tn_type+'_pc1')
    
    print("pcc between hi-c ({}) pc1 and {} {} pc1: {}".format(tissue_type, cohort, tn_type, pearsonr(y1, y2)[0]))
else:
    print("No matching FIRE cohort whose tissue type matches with the selected cohort {}".format(cohort))


# In[ ]:


print("Compare FIRE PC1 and iebdm-normal-pc1")
tn_type = 'normal'
if cohort in tcga_fire_cohort.index.values:
    fire_type = tcga_fire_cohort.loc[cohort].values[0]
    d = pd.read_csv('source_data/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1_comparison_'+tn_type+'.csv')

    fig = plt.figure(figsize=(6, 0.75))
    ax = fig.add_subplot(111)

    ax.axhline(0, lw=0.75, ls='--', c='0.8')

    y1 = (d[fire_type] - d[fire_type].mean()) / d[fire_type].std() #fire pc1
    ax.plot(smoothen(y1, window=3), lw=1, c='0.5', label='Hi-C ({})'.format(tissue_type))

    y2 = d['bmdm_pc1'] 
    y2 = -(y2 - y2.mean()) / y2.std() #standardized iebdm-derived pc1
    ax.plot(smoothen(y2, window=3), lw=1, c='C3', label='Single-sample estimations ({}, averaged)'.format(cohort))

    for i in range(10):
        y3 = d[f'bmdm_sample_{i}']
        y3 = -(y3 - y3.mean()) / y3.std()
        ax.plot(smoothen(y3, window=3), c='C3', lw=0.75, alpha=0.2)

    for direction in ['top', 'right', 'bottom']:
        ax.spines[direction].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    
    assert len(y1)==len(y2) and len(y2)==len(y3)
    ax.set_xlim([0, len(y1)])
    ax.set_ylim([-2.5, 2.5])

    ax.set_xlabel('Position ({})'.format(chrom))
    ax.set_ylabel('PC1\n(Standardized)')

    ax.set_xticks([])
    mb_length = 1 / len(y1)
    ax.plot([0.75, 0.75 + 10*mb_length], [-0.1, -0.1], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
    ax.text(0.75 + 10*mb_length + 0.01, -0.1, '10Mb', transform=ax.transAxes, ha='left', va='center')

    ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.4))

    #save_figures('source_data_figures/'+cohort.split('-')[-1]+'_lung_pc1_comparison')
    #save_figures('source_data_figures/'+cohort.split('-')[-1]+'_'+tissue_type+'_pc1')
    print('source_data_figures/'+cohort.split('-')[-1]+'_'+tissue_type+'_'+tn_type+'_pc1')
    
    print("pcc between hi-c ({}) pc1 and {} {} pc1: {}".format(tissue_type, cohort, tn_type, pearsonr(y1, y2)[0]))
else:
    print("No matching FIRE cohort whose tissue type matches with the selected cohort {}".format(cohort))


# ## 2. 동일 cohort 내에서 tumor와 normal의 pc1 그래프 개형이 차이 난다는 것을 보이기

# ### re-load data

# In[ ]:


NORMAL7_COHORT = 'TCGA-BLCA TCGA-LUAD TCGA-PRAD TCGA-KIRC TCGA-ESCA TCGA-UCEC TCGA-KIRP TCGA-THCA TCGA-HNSC TCGA-LIHC TCGA-LUSC TCGA-CHOL TCGA-PAAD TCGA-BRCA TCGA-COAD'.split(' ')


# In[ ]:


CHR_LIST = ['chr'+str(i) for i in np.arange(1, 23)]


# In[ ]:


for chrom in CHR_LIST:
    print("<{}>".format(chrom))
    for i in range(len(NORMAL7_COHORT)):
        #tissue_type = 'kidney'#previous version

        cohort = NORMAL7_COHORT[i]
        print("===\ncohort: {}".format(cohort))
        tissue_type = tcga_cohort_info.loc[cohort, 'tissue']
        print("Make sure that tissue type {} matches the selected cohort {}\n---".format(tissue_type, cohort))


        #cohort = 'TCGA-LUSC'
        #cohort = 'TCGA-LUAD'
        #files = glob.glob(f'/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/{cohort}/*_inv_exp.npz') #이전 버전
        #-------------
        # 1. Prepare data before plotting.
        files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*_inv_exp.npz') #새로 계산한 버전
        tmp = np.load(files[11])
        len(files)
        print("1. collect 10 tumor samples' pc1 vectors")
        data = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            if sample.endswith('11'):
                continue

            #data.append(np.load(f)[chrom])#이전 버전
            data.append(np.load(f)[chrom+'_pc1'])#현재 버전

        data = np.vstack(data)

        np.random.seed(42)
        idxs = np.random.choice(range(len(data)), size=10)
        data = data[idxs]

        np.save('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy', data)
        print('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')

        print("---\n2. collect 10 normal samples' data.")
        data = []
        for f in files:
            sample = os.path.basename(f).split('_')[0]
            if not sample.endswith('11'):#tcga barcode가 11로 끝나지 않으면. 즉 tumor sample이면 -> 건너뜀.
                continue

            #data.append(np.load(f)[chrom])#이전 버전
            data.append(np.load(f)[chrom+'_pc1'])

        data = np.vstack(data)

        np.random.seed(42)
        idxs = np.random.choice(range(len(data)), size=10)
        data = data[idxs]

        np.save('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy', data)
        print('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')
        #------------------------------
        # 2. plot 10 tumor samples' PC1 vectors
        # Plot PC1 vectors from 10 randomly picked tumor samples of current cohort. 
        # draw 10 tumor samples
        fig = plt.figure(figsize = (4, 0.3*22))

        data = np.load('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')
        #ax = fig.add_subplot(1, 2, 1)
        #fig = plt.figure(figsize=(4, 3))

        for i in range(10):
            ax = fig.add_subplot(22, 1, i+1)
            ax.spines['left'].set_linewidth(0.6)#선 두께 #https://zephyrus1111.tistory.com/241
            if i == 0:
                ax.set_title("{}".format(cohort))
            if i == 4:
                ax.set_ylabel('Tumor PC1', fontsize = 6.5)

            y = smoothen(standardize(data[i]), window=3)
            x = np.arange(len(y))

            ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
            ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)

            for direction in ['top', 'right', 'bottom']:
                ax.spines[direction].set_visible(False)
            ax.spines['left'].set_position(('outward', 3))

            ax.set_xlim([0, len(y)])
            ax.set_xticks([])
            ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)

            ax.tick_params(axis='y', labelsize=4.5, length=2, pad=1, width = 0.6)
            '''
            if i == 9:
                ax.set_xlabel('Position ({})'.format(chrom))

                mb_length = 1 / len(y)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.6, -0.6], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.6, '10Mb', transform=ax.transAxes, ha='left', va='center')
            '''
        #fig.supylabel('PC1', fontsize=7)

        #save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')
        #print('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')

        #------------------------------
        # 3. plot 10 normal samples' PC1 vectors. 
        # Plot PC1 vectors from 10 randomly picked normal samples of current cohort. 
        data = np.load('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')

        #fig = plt.figure(figsize=(4, 3))

        for i in range(10):
            ax = fig.add_subplot(22, 1, i+12)
            ax.spines['left'].set_linewidth(0.6)
            #if i==0:
                #ax.set_title('Normal')
            if i == 4:
                ax.set_ylabel('Normal PC1', fontsize = 6.5)

            y = smoothen(standardize(data[i]), window=3)
            x = np.arange(len(y))

            ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
            ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)

            for direction in ['top', 'right', 'bottom']:
                ax.spines[direction].set_visible(False)
            ax.spines['left'].set_position(('outward', 3))

            ax.set_xlim([0, len(y)])
            ax.set_xticks([])
            ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)

            ax.tick_params(axis='y', labelsize=4.5, length=2, pad=1, width = 0.6)

            if i == 9:
                ax.set_xlabel('Position ({})'.format(chrom))

                mb_length = 1 / len(y)
                ax.plot([0.75, 0.75 + 10*mb_length], [-0.8, -0.8], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
                ax.text(0.75 + 10*mb_length + 0.01, -0.8, '1Mb', transform=ax.transAxes, ha='left', va='center')

        #fig.supylabel('PC1', fontsize=7)
        #plt.suptitle('{}'.format(cohort), va='center')
        fig.tight_layout()
        #



        #save_figures('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')
        #print('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')
        save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_normal_10samples_pc1_'+chrom)
        print('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_normal_10samples_pc1_'+chrom)
        #plt.clf()


# In[ ]:


# single cohort (1)
#tissue_type = 'kidney'#previous version
cohort = NORMAL7_COHORT[0]
tissue_type = tcga_cohort_info.loc[cohort, 'tissue']
print("Make sure that tissue type {} matches the selected cohort {}!!".format(tissue_type, cohort))


# In[ ]:


# single cohort (2)
#cohort = 'TCGA-LUSC'
#cohort = 'TCGA-LUAD'
#files = glob.glob(f'/data/project/jeewon/research/3D-ITH/pipelines/all-samples-pc1/result/{cohort}/*_inv_exp.npz') #이전 버전

files = glob.glob(f'/data/project/3dith/pipelines/opensea-pipeline/1_compute-score-opensea/result/{cohort}/pc1/*_inv_exp.npz') #새로 계산한 버전
tmp = np.load(files[11])
len(files)
print("1. collect 10 tumor samples' data")
data = []
for f in files:
    sample = os.path.basename(f).split('_')[0]
    if sample.endswith('11'):
        continue
        
    #data.append(np.load(f)[chrom])#이전 버전
    data.append(np.load(f)[chrom+'_pc1'])#현재 버전
    
data = np.vstack(data)

np.random.seed(42)
idxs = np.random.choice(range(len(data)), size=10)
data = data[idxs]

np.save('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy', data)
print('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')

print("---\n2. collect 10 normal samples' data.")
data = []
for f in files:
    sample = os.path.basename(f).split('_')[0]
    if not sample.endswith('11'):#tcga barcode가 11로 끝나지 않으면. 즉 tumor sample이면 -> 건너뜀.
        continue
        
    #data.append(np.load(f)[chrom])#이전 버전
    data.append(np.load(f)[chrom+'_pc1'])
    
data = np.vstack(data)

np.random.seed(42)
idxs = np.random.choice(range(len(data)), size=10)
data = data[idxs]

np.save('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy', data)
print('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')


# In[ ]:


# single cohort (3)
# Plot PC1 vectors from 10 randomly picked tumor samples of current cohort. 
# draw 10 tumor samples
data = np.load('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')

fig = plt.figure(figsize=(4, 3))

for i in range(10):
    ax = fig.add_subplot(10, 1, i+1)
    
    y = smoothen(standardize(data[i]), window=3)
    x = np.arange(len(y))
    
    ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
    ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)
    
    for direction in ['top', 'right', 'bottom']:
        ax.spines[direction].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    
    ax.set_xlim([0, len(y)])
    ax.set_xticks([])
    ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)
    
    ax.tick_params(axis='y', labelsize=6, length=2, pad=1)
    
    if i == 9:
        ax.set_xlabel('Position ({})'.format(chrom))

        mb_length = 1 / len(y)
        ax.plot([0.75, 0.75 + 10*mb_length], [-0.6, -0.6], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
        ax.text(0.75 + 10*mb_length + 0.01, -0.6, '10Mb', transform=ax.transAxes, ha='left', va='center')
    
fig.supylabel('PC1', fontsize=7)

save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')
print('source_data_figures/'+cohort.split('-')[1]+'_tumor_10samples_pc1')


# In[ ]:


# single cohort (4)
# Plot PC1 vectors from 10 randomly picked normal samples of current cohort. 
data = np.load('source_data/'+cohort.split('-')[1]+'_normal_10samples.npy')

fig = plt.figure(figsize=(4, 3))

for i in range(10):
    ax = fig.add_subplot(10, 1, i+1)
    
    y = smoothen(standardize(data[i]), window=3)
    x = np.arange(len(y))
    
    ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
    ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)
    
    for direction in ['top', 'right', 'bottom']:
        ax.spines[direction].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    
    ax.set_xlim([0, len(y)])
    ax.set_xticks([])
    ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)
    
    ax.tick_params(axis='y', labelsize=6, length=2, pad=1)
    
    if i == 9:
        ax.set_xlabel('Position ({})'.format(chrom))

        mb_length = 1 / len(y)
        ax.plot([0.75, 0.75 + 10*mb_length], [-0.6, -0.6], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
        ax.text(0.75 + 10*mb_length + 0.01, -0.6, '10Mb', transform=ax.transAxes, ha='left', va='center')
    
fig.supylabel('PC1', fontsize=7)

save_figures('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')
print('source_data_figures/'+cohort.split('-')[1]+'_normal_10samples_pc1')


# In[ ]:


'''
# 라벨, 축 없애고 그래프만 뽑고 싶을 때. 
# draw 1 tumor sample (additional. not necessary to run)
data = np.load('source_data/'+cohort.split('-')[1]+'_tumor_10samples.npy')

fig = plt.figure(figsize=(4, 3/10))


ax = fig.add_subplot(111)

y = smoothen(standardize(data[0]), window=3)#fig1b
#y = smoothen(standardize(data[5]), window=3)#fig1a
x = np.arange(len(y))

ax.fill_between(x, y, y2=0, where=(y >= 0), fc='C3', lw=0)
ax.fill_between(x, y, y2=0, where=(y < 0), fc='C0', lw=0)

for direction in ['top', 'right', 'bottom']:
    ax.spines[direction].set_visible(False)
ax.spines['left'].set_position(('outward', 3))

ax.set_xlim([0, len(y)])
ax.set_xticks([])
ax.axhline(0, c='0.8', lw=0.75, ls='--', zorder=-20)

ax.tick_params(axis='y', labelsize=6, length=2, pad=1)


ax.set_xlabel('Position ({})'.format(chrom))

mb_length = 1 / len(y)
ax.plot([0.75, 0.75 + 10*mb_length], [-0.6, -0.6], lw=0.75, c='k', clip_on=False, transform=ax.transAxes)
ax.text(0.75 + 10*mb_length + 0.01, -0.6, '10Mb', transform=ax.transAxes, ha='left', va='center')
    
fig.supylabel('PC1', fontsize=7)
ax.yaxis.set_visible(False)#y축 눈금 없애기
ax.spines['left'].set_visible(False)
save_figures('source_data_figures/'+cohort.split('-')[1]+'_tumor_1sample_pc1-index5')
print('source_data_figures/'+cohort.split('-')[1]+'_tumor_1sample_pc1-index5')
'''


# In[ ]:


print(os.getcwd())


# In[ ]:





# In[ ]:


fig = plt.figure(figsize=(10, 1))
ax = fig.add_subplot(111)

for i in range(9):
    ax.plot(smoothen(data[i], window=3), lw=1, alpha=0.1)


# In[ ]:


fig = plt.figure(figsize=(10, 1))
ax = fig.add_subplot(211)

ax.plot((d['LG'] - d['LG'].mean()) / d['LG'].std())

ax = fig.add_subplot(212)
y = np.mean(data, axis=0)
y = -(y - y.mean()) / y.std()

ax.fill_between(x=np.arange(len(y)), y1=y, y2=0)

# for y in data[:20]:
#     ax.plot(-y, lw=1, alpha=0.1)


# In[ ]:


for tissue in tissues:
    mask = (np.abs(d[tissue]) > 0.5) & (np.abs(d['bmdm_pc1']) > 0.75)
    print(tissue, stats.spearmanr(d[mask][tissue], d[mask]['bmdm_pc1']))


# In[ ]:


data


# In[ ]:


corr_data = defaultdict(list)

for f in tqdm(files[:50]):
    data = defaultdict(list)
    
    sample = os.path.basename(f).split('_')[0]
    if sample.endswith('11'):
        continue
    
    bmdm_pc1 = np.load(f)
    
    for chrom in [f'chr{i}' for i in range(1, 23)]:
        for bin_id, pc1 in zip(bins[f'{chrom}_bins'], bmdm_pc1[chrom]):
            data['bin_id'].append(bin_id)
            data['bmdm_pc1'].append(pc1)

    data = pd.DataFrame(data)

    data['chrom'] = data.bin_id.str.split(':', expand=True)[0]
    data['start'] = data.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int) + 1
    data['end'] = data.bin_id.str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    data = data[['chrom', 'start', 'end', 'bmdm_pc1']]
    data = fire_pc1.merge(data, on=['chrom', 'start', 'end'])
    
    tissues = ['GM12878', 'H1', 'IMR90', 'MES', 'MSC', 'NPC', 'TRO', 'AD', 'AO', 'BL', 'CO', 'HC', 'LG', 'LI', 'LV', 'OV', 'PA', 'PO', 'RV', 'SB', 'SX']

    for tissue in tissues:
        for chrom in [f'chr{i}' for i in range(1, 23)]:
            mask = (data['bmdm_pc1'] > 0) & (np.abs(data[tissue]) > 0.1) & (data['chrom'] == chrom)
            r, p = stats.spearmanr(data[mask]['bmdm_pc1'], data[mask][tissue])

            corr_data['sample'].append(sample)
            corr_data['tissue'].append(tissue)
            corr_data['chrom'].append(chrom)
            corr_data['r'].append(np.abs(r))
            corr_data['p'].append(p)
        
corr_data = pd.DataFrame(corr_data)


# In[ ]:


fig = plt.figure(figsize=(10, 1))
ax = fig.add_subplot(111)

sns.boxplot(data=corr_data[(corr_data.chrom == 'chr2')], x='tissue', y='r', ax=ax)


# In[ ]:


corr_data[(corr_data.chrom == 'chr12')].groupby('tissue').agg({'r': 'mean'}).sort_values('r')


# In[ ]:


diffmat_cpg_type = 'opensea-binned-diffmat-bins'

bins = np.load(f'/data/project/jeewon/research/3D-ITH/pipelines/find-diffmat-bins/{diffmat_cpg_type}/{cohort}_diffmat_bins.npz')


# In[ ]:






# In[ ]:


data = fire_pc1.merge(data, on=['chrom', 'start', 'end'])


# In[ ]:





# In[ ]:


meta = pd.read_csv('/data/project/3dith/data/450k_metadata.resort.sorted.bed', sep='\t', names=['chrom', 'start', 'end', 'name'])


# In[ ]:





# In[ ]:


fire_pc1


# In[ ]:


fire_pc1[fire_pc1['chr'] == 1]


# In[ ]:


meta.shape


# In[ ]:


meta.head(3)

