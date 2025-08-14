#!/usr/bin/env python
# coding: utf-8

# # Exploration of hormonal receptor genes in Polioudakis et al human fetal brain dataset - Metacell calculation by SeaCells

# [Reference paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6831089/#SD1) : *We focus on the cortical anlage at mid-gestation (gestation week (GW) 17 to 18) because this period contains the major germinal zones and the developing cortical laminae containing migrating and newly born neurons, and neurodevelopmental processes occurring during this epoch are implicated in neuropsychiatric disease. To optimize detection of distinct cell types, prior to single-cell isolation we separated the cortex into:
# the germinal zones (ventricular zone (VZ) and subventricular zone (SVZ)) and developing cortex (subplate (SP) and cortical plate (CP)).*

# __Upstream Steps__
# 
# * Assemble adata
# * QC filter on cells; Expression filter on genes
# * Normalization and log10 transformation by Scanpy functions
# * Feature selection (HVG) by Scanpy functions
# * Dimensionality reduction
# * Cluster identification
# 
# 
# __This Notebook__
# 
# Metacell calculation by SEACells (following [this tutorial](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb))
# 

# -----

# # 1. Environment

# ## 1.1 Modules

# In[3]:


import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

import pickle

#Plotting
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

#utils
import ipynbname
from datetime import datetime

# SeaCell
import SEACells

#import custom functions
sys.path.append('../')
import functions as fn


# In[2]:


# Some plotting aesthetics
get_ipython().run_line_magic('matplotlib', 'inline')

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [3.5, 3.5]
matplotlib.rcParams['figure.dpi'] = 100


# In[3]:


print("Scanpy version: ", sc.__version__)
print("Pandas version: ", pd.__version__)
print("SEACell version: ", SEACells.__version__)


# ## 1.2 Settings
# 

# In[4]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)


# ## 1.3 Files and parameters

# In[5]:


path = '../../../../Polioudakis/'
input_file = path + '3_FiltNormAdata.h5ad'
output_file = '../../../../DataDir/ExternalData/SingleCellData/PolioudakisAdata_Metacells.h5ad'


# ## 1.4 Start computations

# In[6]:


print(datetime.now())


# ----

# # 2. Data Load
# 

# ## 2.1 Read adata file
# 
# * SEACells uses as input filtered unnormalized counts from scRNASeq data, that are subjected to normalization and log-transformation before metacell assignment
# * Here we start from data that has already undergone the steps of QC, filtering, normalization, log transformation
# * The sum of expression levels from cells to metacells will occur at the count level  

# In[7]:


adata = sc.read(input_file)


# In[8]:


adata


# In[9]:


adata.obsm


# In[10]:


adata.X


# In[11]:


print(adata.X[40:45, 40:45])


# In[12]:


adata.layers['counts']


# In[13]:


print(adata.layers['counts'][40:45, 40:45])


# In[14]:


sc.pl.scatter(adata, basis='umap_harmony', color='cell_label', frameon=False)


# In[15]:


col1 = {'Off-Target':'#c5b0d5',
        'ExcitatoryNeu':'#aa40fc',
        'Progenitors&RG':'#ff7f0e',
        'InhibitoryNeu':'#17becf'}

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


# In[16]:


sc.pl.scatter(adata, basis='pca', color='sample_id', frameon=False)


# In[17]:


sc.pl.scatter(adata, basis='pca_harmony', color='sample_id', frameon=False)


# In[18]:


sc.pl.scatter(adata, basis='umap_nocorr', color='sample_id', frameon=False)


# In[19]:


sc.pl.scatter(adata, basis='umap_harmony', color='sample_id', frameon=False)


# In[20]:


adata.obs.head(3)


# ## 2.2 Pre-processing
# 
# Normalization, log-transformation, HVG and dimensionality reduction have been already performed.

# # 3. Run SEACells
# 
# We follow the workflow described in the [SEACells tutorial](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb)
# 
# 

# ## 3.1 Define parameters
# 
# Parameters to be defined: 
# 
# * __Number of metacells__ (graining level): choosing one metacell for every 60 single-cells.
# * __X_pca_harmony__ as the key in .obsm for computing in metacells. PCs already corrected by harmony.
# * __Number of eigenvalues__ for metacell initialization
# 
# We then employ the __SEACells.core.SEACells__ function to initialize the model. 

# In[21]:


## Core parameters 
n_SEACells = round(adata.n_obs/75)
print(n_SEACells)

build_kernel_on = 'X_pca_harmony' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells


# In[22]:


#help(SEACells.core.SEACells)

model = SEACells.core.SEACells(adata, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)


# ## 3.2 Construct kernel matrix 
# 
# [__construct_kernel_matrix__](https://github.com/dpeerlab/SEACells/blob/main/SEACells/cpu.py) function constructs the kernel matrix from data matrix using PCA/SVD and nearest neighbors. Key parameters are: 
# 
# * n_neighbors: (int) number of nearest neighbors to use for graph construction (default 15). Increasing this number results in a more homogeneous distribution of metacell size, but may impact on the preservation of rare cell types. Usually stable performance in the range 5-30. 
# * graph_construction: method for graph construction.

# In[23]:


model.construct_kernel_matrix(n_neighbors=20)
M = model.kernel_matrix


# ## 3.3 Initialize archetypes
# 
# * [__initialize_archetypes__](https://github.com/dpeerlab/SEACells/blob/main/SEACells/cpu.py) function initializes the B matrix which defines cells as SEACells.
# * The function uses waypoint analysis for initialization into to fully cover the phenotype space, and then greedily selects the remaining cells.
# 

# In[24]:


# Initialize archetypes
model.initialize_archetypes()


# In[25]:


adata


# In[26]:


adata.obsm['X_umap'] = adata.obsm['X_umap_harmony'].copy()


# In[27]:


# Plot the initilization to ensure they are spread across phenotypic space
SEACells.plot.plot_initialization(adata, model)


# ## 3.4 Fit Model
# 
# We apply the __model.fit__ function. 

# In[28]:


model.fit(min_iter=10, max_iter=50)
model


# ## 4. Access SEACells results 

# ### 4.1 Check model convergence

# In[29]:


# Plot model convergence
model.plot_convergence()


# ### 4.2 Inspect SEACells soft assignment
# 
# * SEACells analysis returns soft assignments of cells to SEACells (full assignment matrix can be accessed at model.A_)
# * The majority of single-cells are assigned to no more than 4 archetypes with non-trivial weight
# * The algorithm returns the top 5 metacell assignments as well as the corresponding assignment weights in the function __model.get_soft_assignments()__
# 

# In[30]:


#model.A_

plt.figure(figsize=(4,3))
sns.distplot((model.A_.T > 0.1).sum(axis=1), kde=False)
plt.title(f'Non-trivial (> 0.1) assignments per cell')
plt.xlabel('# Non-trivial SEACell Assignments')
plt.ylabel('# Cells')
plt.show()

plt.figure(figsize=(4,3))
b = np.partition(model.A_.T, -5)    
sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
plt.title('Strength of top 5 strongest assignments')
plt.xlabel('$n^{th}$ strongest assignment')
plt.show()


# In[31]:


labels,weights = model.get_soft_assignments()
#labels.head()


# ### 4.3 Inspect and retrieve SEACells hard assignment
# 
# We will use hard assignment of each cell to a single metacell to define our metacells for downstream steps. Hard assingment can be retreived in two ways:
# 
# * in the modified anndata object in .obs['SEACell']
# * from the model using .get_hard_assignments()
# 

# In[32]:


adata.obs[['SEACell']].head()

# Alternatively: 
# model.get_hard_assignments().head()


# ## 5. Summarizing data by hard assignment
# 
# * The original dataset can be summarized by SEACell by aggregating cells within each SEACell, summing over all raw data for all cells belonging to a SEACell.  
# * The output of this function is an anndata object of shape n_metacells x original_data_dimension.
# * Data is unnormalized and raw aggregated counts are stored in X.
# * Attributes associated with variables (.var) are copied over, but relevant per SEACell attributes must be manually copied, since certain attributes may need to be summed, or averaged etc, depending on the attribute.

# In[33]:


# By default, ad.raw is used for summarization. Other layers present in the anndata can be specified using the parameter summarize_layer parameter
SEA_adata = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='counts', celltype_label='cell_label')
SEA_adata


# In[34]:


SEA_adata.obs.head()


# In[35]:


print(SEA_adata.X[40:43, 40:43])


# ## 6. Summarize metadata (for hard assignment)
# 
# * Metadata that is available at the cell level can be summarized as an annotation at the metacell level
# * Countinous variables are averaged across the cells assigned to each metacell
# * Categorical variables are determined by the most prevalent category amongst the cells assigned to the metacell

# In[36]:


top_sampleid = adata.obs['sample_id'].groupby(adata.obs['SEACell']).value_counts().groupby(level=0, group_keys=False).head(1) 
SEA_adata.obs['agg_sample_id'] = top_sampleid[SEA_adata.obs_names].index.get_level_values(1)


# ## 7. Visualize results
# 
# * __plot_2D()__:  visualize metacell assignments on UMAP or any 2-dimensional embedding froma adata.obsm. Plots can also be coloured by metacell assignment. The visualization allows to assess the __representativeness__ of the calculated metacells: ability to represent the global structure of the original dataset. Better representation corresponds to a more uniform coverage of the dataset.
# 
# * __plot_SEACell_sizes()__: distribution of number of cells assigned to each metacell. Given to the non-uniform distribution of cell type density, metacells are expected to vary in size: larger metacells are retrieved from denser regions, and smaller metacells from sparser ones. However outlier values may require further inspection. 
# 
# 

# In[37]:


SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=False)


# In[38]:


SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=True)


# In[39]:


SEACells.plot.plot_SEACell_sizes(adata, bins=10, figsize=(10, 5))


# ## 8. Quantifying results
# 
# * __Purity__: compute_celltype_purity(ad, col_name) computes the purity of different celltype labels within a SEACell metacell. 
#   
# * __Compactness__: per-SEAcell variance in diffusion components (typically 'X_pca' for RNA). Lower values of compactness suggest more compact/lower variance metacells.
# 
# * __Separation__: distance between a SEACell and its nth_nbr. If cluster is provided as a string, e.g. 'celltype', nearest neighbors are restricted to have the same celltype value. Higher values of separation suggest better distinction between metacells. 

# ### 8.1 Purity
# 
# #### Cell Population

# In[40]:


SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata, 'cell_label')

plt.figure(figsize=(6,4))
sns.boxplot(data=SEACell_purity, y='cell_label_purity')
plt.title('cell_label Purity')
sns.despine()
plt.show()
plt.close()

SEACell_purity.head()


# In[41]:


plt.figure(figsize=(6,4))
sns.boxplot(data=SEACells.evaluate.compute_celltype_purity(adata, 'super_cell_label'), y='super_cell_label_purity')
plt.title('super_cell_label Purity')
sns.despine()
plt.show()
plt.close()


# #### Diagnosis
plt.figure(figsize=(6,4))
sns.boxplot(data=SEACells.evaluate.compute_celltype_purity(adata, 'diagnosis'), y='diagnosis_purity')
plt.title('Diagnosis Purity')
sns.despine()
plt.show()
plt.close()
# #### Brain Region
plt.figure(figsize=(6,4))
sns.boxplot(data=SEACells.evaluate.compute_celltype_purity(adata, 'brainRegion'), y='brainRegion_purity')
plt.title('Brain Region Purity')
sns.despine()
plt.show()
plt.close()
# ### 8.2 Compactness

# In[42]:


compactness = SEACells.evaluate.compactness(adata, 'X_pca')

plt.figure(figsize=(6,4))
sns.boxplot(data=compactness, y='compactness')
plt.title('Compactness')
sns.despine()
plt.show()
plt.close()

compactness.head()


# ### 8.3 Separation

# In[43]:


separation = SEACells.evaluate.separation(adata, 'X_pca')

plt.figure(figsize=(6,4))
sns.boxplot(data=separation, y='separation')
plt.title('Separation')
sns.despine()
plt.show()
plt.close()

separation.head()


# ## 9. Hard assignment: downstream processing
# 
# * normalization and log-transformation
# * dimensionality reduction: PCA and UMAP on SEACells

# In[44]:


SEA_adata.layers['counts'] = SEA_adata.X.copy()

sc.pp.normalize_total(SEA_adata)
sc.pp.log1p(SEA_adata)
SEA_adata.layers['lognorm'] = SEA_adata.X.copy()

sc.pp.highly_variable_genes(SEA_adata, inplace=True, batch_key="agg_sample_id")
sc.tl.pca(SEA_adata, use_highly_variable=True)
sc.pp.neighbors(SEA_adata, use_rep='X_pca')
sc.tl.umap(SEA_adata)


# In[45]:


sc.pl.umap(SEA_adata, color=['cell_label'], s=75)


# In[46]:


sc.pl.umap(SEA_adata, color=['cell_label_purity'], s=75)


# In[47]:


SEA_adata.obs.columns


# ----------------
# 
# # 11. Saving

# ## 11.1 Save SEACell Adata

# In[53]:


SEA_adata


# In[54]:


SEA_adata.X


# In[55]:


SEA_adata.write(output_file)


# ## 11.2 Save in other formats

# In[4]:


nb_fname = ipynbname.name()
nb_fname


# In[5]:


get_ipython().run_cell_magic('bash', '-s "$nb_fname"', 'jupyter nbconvert "$1".ipynb --to="python"\njupyter nbconvert "$1".ipynb --to="html"\n')


# ## 11.3 Finished computations: timestamp

# In[56]:


print(datetime.now())

