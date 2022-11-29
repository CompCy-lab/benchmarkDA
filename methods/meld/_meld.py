import anndata
import numpy as np
import pandas as pd
import scanpy as sc

import graphtools as gt
import meld
from sklearn.preprocessing import normalize


def replicate_normalize_densities(sample_densities, replicate):
    sample_likelihoods = sample_densities.copy()
    for rep in replicate:
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        sample_likelihoods[curr_cols] = normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods


def normalize_densities_real(sample_densities, samplem, replicate, replicate_col):
    sample_likelihoods = sample_densities.copy()
    for rep in replicate:
        curr_cols = sample_densities.columns[
            [
                # col.endswith(rep)
                rep == samplem.loc[col, replicate_col]
                for col in sample_densities.columns
            ]
        ]
        sample_likelihoods[curr_cols] = normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods


def runMELD_real(
    adata: anndata.AnnData,
    k: int,
    tgt_label: str,
    sample_col: str,
    label_col: str,
    replicate_col: str = None,
    beta=20
):
    # add sample and label dataframe to adata
    samplem = pd.DataFrame(index=pd.Series(adata.obs[sample_col]).unique())
    samplem.loc[:,label_col] = \
        adata.obs[[sample_col, label_col]].groupby(by=sample_col).aggregate(lambda x: x[0])
    if replicate_col != None:
        samplem.loc[:,replicate_col] = \
            adata.obs[[sample_col, replicate_col]].groupby(by=sample_col).aggregate(lambda x: x[0])
    adata.uns['samplem'] = samplem

    # run meld method
    if adata.n_vars <= 50:
        G = gt.Graph(adata.X, knn=k, use_pygsp=True)
    else:
        if 'X_pca' not in adata.obsm:
            # perform pca and use it to generate graph
            sc.tl.pca(adata, n_comps=50)
        G = gt.Graph(adata.obsm['X_pca'], knn=k, use_pygsp=True)
    
    meld_op = meld.MELD(beta=beta)
    # generate the densities of each sample
    if replicate_col:
        # normalize the densities for each replicate
        replicates = samplem[replicate_col].unique()
        # w/ replicates, run MELD on sample level
        sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[sample_col])
        sample_likelihoods = normalize_densities_real(sample_densities, samplem, replicates, replicate_col)
        tgt_columns = samplem.loc[samplem[label_col] == tgt_label,:].index.to_list()
        sample_likelihoods = sample_likelihoods[tgt_columns].mean(axis=1)
    else:
        # w/o replicates, run MELD on label level
        sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[label_col])
        sample_likelihoods = sample_densities.div(sample_densities.sum(axis=1), axis=0)
        sample_likelihoods = sample_likelihoods[tgt_label]
    return sample_likelihoods.values


def runMELD(adata: anndata.AnnData, k: int, sample_col: str, label_col: str, beta=20):
    # add sample and label dataframe to adata
    samplem = pd.DataFrame(index=pd.Series(adata.obs[sample_col]).unique())
    samplem.loc[:,label_col] = \
        adata.obs[[sample_col, label_col]].groupby(by=sample_col).aggregate(lambda x: x[0])
    adata.uns['samplem'] = samplem

    # run meld method
    if adata.n_vars <= 50:
        G = gt.Graph(adata.X, knn=k, use_pygsp=True)
    else:
        if 'X_pca' not in adata.obsm:
            # perform pca and use it to generate graph
            sc.tl.pca(adata, n_comps=50)
        G = gt.Graph(adata.obsm['X_pca'], knn=k, use_pygsp=True)
    
    meld_op = meld.MELD(beta=beta)
    # generate the densities of each sample
    sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[sample_col])
    # normalize the densities for each replicate
    replicates = samplem.index.map(lambda x: x.split('_')[-1]).unique()
    sample_likelihoods = replicate_normalize_densities(sample_densities, replicates)
    # average the likelihoods w.r.t conditions
    obj_cond = sorted(samplem[label_col].unique())[-1]
    obj_cond_columns = samplem.loc[samplem[label_col] == obj_cond,:].index.to_list()
    sample_likelihoods = sample_likelihoods[obj_cond_columns].mean(axis=1)
    return sample_likelihoods.values


def meld2output(meld_res, out_type="continuous", thresholds=None):
    if out_type == "continuous":
        da_cell = meld_res
    else:
        def get_da_cell(thres):
            isPos = meld_res > thres
            isNeg = meld_res < 1 - thres
            
            da = np.array(["NotDA"] * len(meld_res), dtype=object)
            da[isPos] = "PosLFC"
            da[isNeg] = "NegLFC"
            return da
        
        if isinstance(thresholds, float):
            da_cell = get_da_cell(thresholds)
        elif isinstance(thresholds, list) or isinstance(thresholds, np.ndarray):
            da_cell = []
            for _, thres in enumerate(thresholds):
                da_cell.append(get_da_cell(thres))
        else:
            raise RuntimeError("param: alphas can only support list or float")
    
    return da_cell
