#!/usr/bin/env python
# coding: utf-8

# # Exploration of hormonal receptor genes in Castaldi et al organoid dataset

# [Reference paper](https://doi.org/10.1038/s41592-024-02555-5)

# ## 1. Environment Set Up

# ### 1.1 Library upload

# In[3]:


import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import igraph as ig
import matplotlib.pyplot as plt 
from scipy.sparse import csr_matrix, isspmatrix
from datetime import datetime
import sys
sys.path.append('../')
import functions as fn

print(np.__version__)
print(pd.__version__)
print(sc.__version__)


# In[4]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)


# ### 1.2 Starting computations: timestamp

# In[5]:


print(datetime.now())


# ## 2. Read input files

# ### 2.1 adata loading

# In[8]:


adata = sc.read('../../../../Castaldi_multiplexingCBO/adataPaga.h5ad')


# In[9]:


adata


# In[11]:


print('Loaded Normalizes AnnData object: number of cells', adata.n_obs)
print('Loaded Normalizes AnnData object: number of genes', adata.n_vars)

# To see the columns of the metadata (information available for each cell)  
print('Available metadata for each cell: ', adata.obs.columns)


# In[12]:


np.unique(adata.obs.values[:,14])


# ### 2.2 Receptors signature loading

# Loading of hormonal receptor gene signature.

# In[13]:


signatures = '../../../../DataDir/ExternalData/Receptors/ReceptorsComplete.txt'


# In[14]:


sig = pd.read_csv(signatures, sep="\t", keep_default_na=False)  
print(sig.shape)
sig


# In[15]:


genes = sig["GeneName"].values.tolist()


# ## 3. Visualizations

# ### 3.1 Counts from adata

# In[16]:


adata.obsm


# In[17]:


sc.pl.embedding(adata, basis="X_umap", color=['n_genes_by_counts',"total_counts", 'pct_counts_mt', 'pct_counts_ribo'])


# ### 3.2 Clusters annotation

# In[18]:


sc.pl.embedding(adata,  basis="X_umap", color=['leidenAnnotated'], ncols=1)


# In[19]:


sc.pl.embedding(adata,  basis="X_draw_graph_fa", color=['leidenAnnotated'], ncols=1)


# ### 3.3 Visualization of receptors on UMAP

# In[20]:


fn.CustomUmap(adata, genes, embedding="X_umap") 


# In[21]:


fn.CustomUmap(adata, genes, embedding="X_draw_graph_fa") 


# In[26]:


available_genes = [gene for gene in genes if gene in adata.var_names]

if available_genes:
    sc.pl.dotplot(adata, available_genes, groupby='cell_label')
else:
    print("None of the specified genes are found in adata.var_names.")
sc.pl.dotplot(adata, genes, groupby='leidenAnnotated')


# ## 4. Save Notebooks

# Adata is not saved, since no new computations have been performed. I just save the notebooks.

# ### 4.1 Timestamp finished computations

# In[24]:


print(datetime.now())


# ### 4.2 Save python and html version of notebook

# In[25]:


get_ipython().run_cell_magic('bash', '', '\n# save also html and python versions for git\njupyter nbconvert ExplorationCastaldi.ipynb --to="python" --output="ExplorationCastaldi"\njupyter nbconvert ExplorationCastaldi.ipynb --to="html" --output="ExplorationCastaldi"\n')


# In[ ]:




