import argparse
import anndata
import numpy as np
import pandas as pd
import os.path as osp

from _meld import runMELD_real, meld2output


def parse_args():
    parser = argparse.ArgumentParser("run meld method on the real data.")

    parser.add_argument("--adata", type=str, help="path to .h5ad data file")
    parser.add_argument("--data_id", type=str, help="dataset id (name)")
    parser.add_argument("--k", type=int, default=30, help="K parameter to build KNN graph")
    parser.add_argument("--sample_col", type=str, help="column to specify different samples")
    parser.add_argument("--label_col", type=str, help="column to specify different labels")
    parser.add_argument("--tgt_label", type=str, help="target label used to compute the likelihood")
    parser.add_argument("--rep_col", type=str, default=None, help="column to specify different replicates")
    parser.add_argument("--beta", type=float, default=20, help="the smoothing parameter for density estimation")
    parser.add_argument("--thres", type=float, default=0.05, help="prob threshold to determine DAs")
    parser.add_argument("--seed", type=int, default=43, help="random seed")
    parser.add_argument("--outdir", type=str, help="dir to save benchmark results")

    args = parser.parse_args()
    if args.rep_col == 'NULL':
        args.rep_col = None
    return args


def runDA(args):
    # set random seed for reproducibility
    np.random.seed(args.seed)
    # read the data
    try:
        adata = anndata.read(args.adata)
    except FileNotFoundError as e:
        print(e)

    # run analysis using CNA method
    meld_res = runMELD_real(adata,
                            k=args.k,
                            sample_col=args.sample_col,
                            label_col=args.label_col,
                            tgt_label=args.tgt_label,
                            beta=args.beta)
    # get da cells with the threshold, alpha.
    da_score = meld2output(meld_res, 'continuous')
    da_cell  = meld2output(meld_res, 'label', args.thres)
    
    # prepare the output dataframe
    DA_df = pd.DataFrame(
        {
            'DA': da_cell,
            'Cond_Prob': da_score
        },
        index=adata.obs.index
    )
    prefix = "benchmark_{}_k={}_beta={}_thres={}_seed={}".format(args.data_id,
                                                                 args.k,
                                                                 args.beta,
                                                                 args.thres,
                                                                 args.seed)
    DA_file = osp.join(args.outdir, prefix + '.DAresults.meld.csv')
    DA_df.to_csv(DA_file, index=True)


if __name__ == "__main__":
    args = parse_args()
    runDA(args)