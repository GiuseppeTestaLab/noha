#!/usr/bin/env python
# coding: utf-8

# # Heatmaps showing EndPoints DEGs 

# ## 1. Environment Set Up

# In[1]:


import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from datetime import datetime


# In[2]:


print(pd.__version__)
print(sns.__version__)


# In[3]:


import os 
cwd = os.getcwd()


# ## 2. Heatmap plots

# ### 2.1 CTL04

# All comparisons

# In[4]:


heatmap = pd.read_excel("~/DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Input/DEGs_CTL04.xlsx", index_col = 0)
heatmap


# In[5]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'top'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 60)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL04_withAgvsInh_h.svg", dpi = 300, transparent = False)
plt.show()


# In[6]:


plt.figure(figsize=(35,15))
ax = sns.heatmap(heatmap.T, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'right'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 0)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL04_withAgvsInh_v.svg", dpi = 300, transparent = False)
plt.show()


# Without AgvsInh

# In[7]:


columns_to_keep = [col for col in heatmap.columns if "AgvsInh" not in col]
heatmap2 = heatmap[columns_to_keep]
heatmap2


# In[8]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap2, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'top'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 60)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL04_h.svg", dpi = 300, transparent = False)
plt.show()


# In[9]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap2.T, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'right'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 0)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL04_v.svg", dpi = 300, transparent = False)
plt.show()


# ### 2.2 CTL08

# In[10]:


heatmap3 = pd.read_excel("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Input/DEGs_CTL08.xlsx", index_col = 0)
heatmap3


# In[11]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap3, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'top'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 60)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL08_withAgvsInh_h.svg", dpi = 300, transparent = False)
plt.show()


# In[12]:


plt.figure(figsize=(35,15))
ax = sns.heatmap(heatmap3.T, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'right'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 0)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL08_withAgvsInh_v.svg", dpi = 300, transparent = False)
plt.show()


# Without AgvsInh

# In[13]:


columns_to_keep = [col for col in heatmap3.columns if "AgvsInh" not in col]
heatmap4 = heatmap3[columns_to_keep]
heatmap4


# In[14]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap4, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'top'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 60)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL08_h.svg", dpi = 300, transparent = False)
plt.show()


# In[15]:


plt.figure(figsize=(16,10))
ax = sns.heatmap(heatmap4.T, linewidth = 0.1, annot = True, annot_kws = {'fontsize': 12}, cmap = 'YlOrBr', square = True, fmt = 'g', cbar_kws = dict(use_gridspec = False, location = 'right'), vmax = 1000)
ax.xaxis.set_tick_params(labelsize=14, rotation = 0)
ax.yaxis.set_tick_params(labelsize=14, rotation = 0)
plt.savefig("../../../DataDir/bulkRNASeq/5.DifferentialExpression/Heatmaps/Output/Heatmap_degs_CTL08_v.svg", dpi = 300, transparent = False)
plt.show()


# ## 3. Save Notebooks

# ### 3.1 Timestamp finished computations

# In[16]:


print(datetime.now())


# ### 3.2 Save python and html version of notebook

# In[17]:


get_ipython().run_cell_magic('bash', '', '\n# save also html and python versions for git\njupyter nbconvert Heatmap_degs.ipynb --to="python" --output="Heatmap_degs"\njupyter nbconvert Heatmap_degs.ipynb --to="html" --output="Heatmap_degs"\n')

