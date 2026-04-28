import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

### MARKER GENES               
def selectMarkers(adata, mList, var_col=None):

    p = []
    found = set()

    for g in mList:
        if g in adata.var_names:
            p.append(g)
            found.add(g)

        elif var_col is not None and var_col in adata.var.columns:
            if (adata.var[var_col] == g).any():
                p.append(g) 
                found.add(g)

    missing = set(mList) - found
    if len(missing) == len(mList):
        print("\nAll markers are missing")
    elif missing:
        print("\nMissing:", missing)

    return p


def CustomUmap(adata, genes, embedding, s=10, var_col=None, **kwargs):
    genes = selectMarkers(adata, genes, var_col=var_col)
    sc.pl.embedding(adata, basis=embedding, color=genes, ncols=1,
                    size=s, frameon=False, vmin='p1', vmax='p99', **kwargs)