#!/usr/bin/env python
# coding: utf-8

# # Exploration of hormonal receptor genes in Polioudakis et al human fetal brain dataset

# [Reference paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6831089/#SD1) : *We focus on the cortical anlage at mid-gestation (gestation week (GW) 17 to 18) because this period contains the major germinal zones and the developing cortical laminae containing migrating and newly born neurons, and neurodevelopmental processes occurring during this epoch are implicated in neuropsychiatric disease. To optimize detection of distinct cell types, prior to single-cell isolation we separated the cortex into:
# the germinal zones (ventricular zone (VZ) and subventricular zone (SVZ)) and developing cortex (subplate (SP) and cortical plate (CP)).*

# ## 1. Environment Set Up

# ### 1.1 Library upload

# In[2]:


import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import igraph as ig
import matplotlib.pyplot as plt 
from scipy.sparse import csr_matrix, isspmatrix
from datetime import datetime
sys.path.append('../')
import functions as fn

print(np.__version__)
print(pd.__version__)
print(sc.__version__)


# In[3]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)


# ### 1.2 Starting computations: timestamp

# In[4]:


print(datetime.now())


# ## 2. Read input files

# ### 2.1 adata loading

# GENERAL INFO before filtering:
# - Samples: for all donors profiled DNA was acquired from fetal brain tissue (obtained from the UCLA Gene and Cell Therapy Core)
# - Sequencing method: Drop-seq (Macosko et al., 2015)
# - Obtained number of cells: ~40,000

# In[5]:


adata = sc.read('../../../../Polioudakis/3_FiltNormAdata.h5ad')


# In[6]:


adata


# In[7]:


isspmatrix(adata.X)


# In[8]:


print('Loaded Normalizes AnnData object: number of cells', adata.n_obs)
print('Loaded Normalizes AnnData object: number of genes', adata.n_vars)

# To see the columns of the metadata (information available for each cell)  
print('Available metadata for each cell: ', adata.obs.columns)


# ### 2.2 Receptors signature loading

# Loading of hormonal receptor gene signature.

# In[9]:


signatures = '../../../../DataDir/ExternalData/Receptors/ReceptorsComplete.txt'


# In[10]:


sig = pd.read_csv(signatures, sep="\t", keep_default_na=False)  #keep_default_na=False: remove Na values
print(sig.shape)
sig


# In[11]:


genes = sig["GeneName"].values.tolist()


# ## 3. Visualizations

# ### 3.1 Counts from adata

# In[12]:


adata.obsm


# In[13]:


sc.pl.embedding(adata, basis="X_umap_harmony", color=['n_genes_by_counts',"total_counts", 'pct_counts_mito', 'pct_counts_ribo'])


# ### 3.2 Clusters annotation

# In[14]:


sc.pl.embedding(adata,  basis="X_umap_harmony", color=['sample_id', 'cell_label'], ncols=1)


# In[15]:


sc.pl.embedding(adata,  basis="X_fa_harmony", color=['sample_id', 'cell_label'], ncols=1)


# In[16]:


col1 = {'Off-Target':'#c5b0d5',
        'ExcitatoryNeu':'#aa40fc',
        'Progenitors&RG':'#ff7f0e',
        'InhibitoryNeu':'#17becf'}


# In[17]:


sup_dictG  = {'End':'Off-Target',
               'ExDp1': 'ExcitatoryNeu',
               'ExDp2': 'ExcitatoryNeu', 
               'ExM': 'ExcitatoryNeu', 
               'ExM-U': 'ExcitatoryNeu', 
               'ExN': 'ExcitatoryNeu', 
               'IP': 'Progenitors&RG', 
               'InCGE': 'InhibitoryNeu', 
               'InMGE': 'InhibitoryNeu', 
               'Mic': 'Off-Target', 
               'OPC': 'Off-Target', 
               'Per':'Off-Target', 
               'PgG2M':'Progenitors&RG', 
               'PgS':'Progenitors&RG', 
               'oRG':'Progenitors&RG', 
               'vRG':'Progenitors&RG'}

#Crate aggregated annotation
adata.obs['super_cell_label'] = adata.obs['cell_label'].replace(sup_dictG)


# In[18]:


sc.pl.embedding(adata,  basis="X_umap_harmony", color=['super_cell_label'], ncols=1, palette=col1)


# ### 3.3 Visualization of receptors on UMAP

# In[19]:


fn.CustomUmap(adata, genes, embedding="X_umap_harmony") 


# In[20]:


fn.CustomUmap(adata, genes, embedding="X_fa_harmony") 


# In[25]:


available_genes = [gene for gene in genes if gene in adata.var_names]

if available_genes:
    sc.pl.dotplot(adata, available_genes, groupby='cell_label')
else:
    print("None of the specified genes are found in adata.var_names.")


# ## 4. Save Notebooks

# Adata is not saved, since no new computations have been performed. I just save the notebooks.

# ### 4.1 Timestamp finished computations

# In[26]:


print(datetime.now())


# ### 4.2 Save python and html version of notebook

# In[27]:


get_ipython().run_cell_magic('bash', '', '\n# save also html and python versions for git\njupyter nbconvert ExplorationPolioudakis.ipynb --to="python" --output="ExplorationPolioudakis"\njupyter nbconvert ExplorationPolioudakis.ipynb --to="html" --output="ExplorationPolioudakis"\n')


# In[ ]:




