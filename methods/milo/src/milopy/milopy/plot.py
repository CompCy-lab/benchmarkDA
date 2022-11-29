import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData

import matplotlib.pyplot as plt
import seaborn as sns

def plot_nhood_graph(
    adata, 
    alpha=0.1,
    min_logFC=0,
    min_size=10, 
    plot_edges=False,
    title = "DA log-Fold Change",
    **kwargs
    ):
    '''
    Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)
    
    - adata: AnnData object
    - alpha: significance threshold
    - min_logFC: minimum absolute log-Fold Change to show results (default: 0, show all significant neighbourhoods)
    - min_size: minimum size of nodes in visualization (default: 10)
    - plot_edges: boolean indicating if edges for neighbourhood overlaps whould be plotted (default: False)
    - title: plot title (default: 'DA log-Fold Change')
    - **kwargs: other arguments to pass to scanpy.pl.embedding
    '''
    nhood_adata = adata.uns["nhood_adata"].copy()

    nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] > alpha, "graph_color"] = np.nan
    nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
    nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"] < min_logFC, "graph_color"] = np.nan
    
    ## Plotting order - extreme logFC on top
    nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(), "abs_logFC"] = np.nan
    ordered = nhood_adata.obs.sort_values('abs_logFC', na_position='first').index
    nhood_adata = nhood_adata[ordered]

    vmax = np.max([nhood_adata.obs["graph_color"].max(), abs(nhood_adata.obs["graph_color"].min())])
    vmin = - vmax

    sc.pl.embedding(nhood_adata, "X_milo_graph", 
                    color="graph_color", cmap="RdBu_r",
                    size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size, 
                    edges=plot_edges, neighbors_key="nhood",
                    # edge_width = 
                    sort_order=False,
                    frameon=False,
                    vmax=vmax, vmin=vmin,
                    title=title,
                    **kwargs
                   )
    
def plot_nhood(adata, ix, basis="X_umap"):    
    '''
    Visualize cells in a neighbourhood
    '''
    adata.obs["Nhood"] = adata.obsm["nhoods"][:,ix].toarray().ravel()
    sc.pl.embedding(adata, basis, color="Nhood", size=30, title="Nhood" + str(ix))
    
### Plot nhood beeswarm
def plot_DA_beeswarm(
    adata, 
    anno_col = 'nhood_annotation', 
    alpha=0.1,
    subset_nhoods=None
):
    '''
    Plot beeswarm plot of logFC against nhood labels
    '''
    try:
        nhood_adata = adata.uns["nhood_adata"].copy()
    except KeyError:
        raise KeyError(
                'Cannot find "nhood_adata" slot in adata.uns -- please run milopy.make_nhoods_adata(adata)'
            )

    if subset_nhoods is not None:
        nhood_adata = nhood_adata[subset_nhoods]
    
    try:
        nhood_adata.obs[anno_col]
    except KeyError:
        raise KeyError(
                'Cannot find {a} in adata.uns["nhood_adata"].obs -- \
                please run milopy.utils.annotate_nhoods(adata, anno_col) first'.format(a=anno_col)
            )
    
    try:
        nhood_adata.obs["logFC"]
    except KeyError:
        raise KeyError(
                'Cannot find `logFC` in adata.uns["nhood_adata"].obs -- \
                please run milopy.core.DA_nhoods(adata) first'
            )
    
    sorted_annos = nhood_adata.obs[[anno_col,"logFC"]].\
        groupby(anno_col).\
        median().\
        sort_values("logFC", ascending=True).index

    anno_df = nhood_adata.obs[[anno_col, "logFC", "SpatialFDR"]]
    anno_df['is_signif'] = anno_df['SpatialFDR'] < alpha
    # anno_df['across_organs'] = ["Significant across organs" if x else "" for x in (keep_nhoods & signif_nhoods)]
    anno_df = anno_df[anno_df[anno_col]!="nan"]

    sns.violinplot(data=anno_df, y=anno_col, x="logFC", order=sorted_annos, 
                   size=190, cmap="Set1", inner=None, orient="h",
                   linewidth=0,
                  scale="width")
    sns.stripplot(data=anno_df, y=anno_col, x="logFC", order=sorted_annos, size=2, 
                  color="black",
                  orient="h", alpha=0.5);
    plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")