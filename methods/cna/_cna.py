import sys
import warnings
import anndata
import numpy as np
import scanpy as sc
import os.path as osp

sys.path.append(
    osp.join(osp.expanduser("~"),
    "Documents/proj/benchmarkDA/methods/cna/src/cna/src")
)
import cna
from multianndata import MultiAnnData


def runCNA(
    adata: anndata.AnnData,
    k: int,
    sample_col: str,
    label_col: str,
    encode_label_dict: dict = {'Condition1': 0, "Condition2": 1},
    batch_col: str = None,
):
    # build the kNN graph
    sc.pp.neighbors(adata, n_neighbors=k)
    # create multi-anndata and convert condition and batch labels to numeric vars
    adata.obs['sample_id'] = adata.obs[sample_col].astype('category').cat.codes + 1
    adata.obs['label_id'] = adata.obs[label_col].map(encode_label_dict).astype(int)
    if batch_col:
        adata.obs['batch_id'] = adata.obs[batch_col].astype('category').cat.codes
    md = MultiAnnData(adata, sampleid='sample_id', dtype=np.float64)
    md.obs_to_sample(["label_id", "batch_id"] if batch_col else ["label_id"])
    # association test
    cna_res = cna.tl.association(md,
                                 getattr(md.samplem, "label_id"),
                                 covs=None,
                                 batches=getattr(md.samplem, "batch_id", None))
    return cna_res, md


def cna2output(md, cna_res, out_type="continuous", alphas=None):
    """
    Get the DA cells from the results of cna object
    """
    if out_type == "continuous":
        da_cell = np.repeat(0., len(md))
        da_cell[cna_res.kept] = cna_res.ncorrs
    else:
        def get_da_cell(alpha):
            passed = cna_res.fdrs[cna_res.fdrs.fdr <= alpha]
            if len(passed) == 0:
                warnings.warn("no neighborhoods were significant at FDR < {}".format(alpha))
                thresh = np.infty
            else:
                thresh = passed.threshold.iloc[0]
            
            isPos = np.repeat(False, len(md))
            isNeg = np.repeat(False, len(md))
            isPos[cna_res.kept] = cna_res.ncorrs > thresh
            isNeg[cna_res.kept] = cna_res.ncorrs < - thresh
            
            da = np.array(["NotDA"] * len(md), dtype=object)
            da[isPos] = "PosLFC"
            da[isNeg] = "NegLFC"
            return da
        
        if isinstance(alphas, float):
            da_cell = get_da_cell(alphas)
        elif isinstance(alphas, list) or isinstance(alphas, np.ndarray):
            da_cell = []
            for _, alpha in enumerate(alphas):
                da_cell.append(get_da_cell(alpha))
        else:
            raise RuntimeError("param: alphas can only support list or float")

    return da_cell