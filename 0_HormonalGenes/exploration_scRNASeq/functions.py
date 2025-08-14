import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

### MARKER GENES               
def selectMarkers(adataObj, mList):
    """  From a list of gene names select only the genes that are present in adata.var
    """
    #Select markers present in adata
    p = adataObj.var_names[adataObj.var_names.isin(mList) == True]
    #Keep the same order as input list
    p = [x for x in mList if x in p]   
    
    #Select missing genes
    ab = set(mList).difference(set(adataObj.var_names))
    
    #Print message SISTEMA
    if len(ab) == len(mList):
        print('\nAll markers are missing')
    else:
        print('\nThe following marker genes are missing: ', ab)
        
    return(p)

def CustomUmap(adata, genes, embedding, s=10):
    genes = selectMarkers(adata, genes)
    #sc.pl.umap(adata, color=genes, size=10, frameon=False,vmin='p1',  vmax='p99')
    sc.pl.embedding(adata,  basis=embedding, color=genes, ncols=1, size=s, frameon=False, vmin='p1',  vmax='p99')