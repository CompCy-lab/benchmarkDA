import warnings
import argparse
import anndata
import numpy as np
import pandas as pd
import os.path as osp

from _cna import runCNA, cna2output


def parse_args():
    parser = argparse.ArgumentParser("run cna method on the real data.")

    parser.add_argument("--adata", type=str, help="path to .h5ad data file")
    parser.add_argument("--data_id", type=str, help="dataset id (name)")
    parser.add_argument("--k", type=int, default=30, help="K parameter to build KNN graph")
    parser.add_argument("--sample_col", type=str, help="column to specify different samples")
    parser.add_argument("--label_col", type=str, help="column to specify different labels")
    parser.add_argument("--batch_col", type=str, default=None, help="whether add batch in the model")
    parser.add_argument(
        '--encode_label_dict',
        type=lambda x: {k: int(v) for k,v in (i.split('=') for i in x.split(','))},
        help='comma-separated key:value pairs to specify the condition encoders, e.g. Condition1=0,Condition2=1'
    )
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR threshold to determine DAs")
    parser.add_argument("--seed", type=int, default=43, help="random seed")
    parser.add_argument("--outdir", type=str, help="dir to save benchmark results")

    args = parser.parse_args()
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
    cna_res, md = runCNA(adata,
                         k=args.k,
                         sample_col=args.sample_col,
                         label_col=args.label_col,
                         encode_label_dict=args.encode_label_dict,
                         batch_col=args.batch_col if args.batch_col else None)
    if cna_res.p > .05:
        warnings.warn("Global association p-value: {} > .05".format(cna_res.p))
    # get da cells with the threshold, alpha.
    da_score = cna2output(md, cna_res, 'continuous')
    da_cell  = cna2output(md, cna_res, 'label', args.alpha)
    
    # prepare the output dataframe
    DA_df = pd.DataFrame(
        {
            'DA': da_cell,
            'Cond_Corr': da_score
        },
        index=adata.obs.index
    )
    prefix = "benchmark_{}_k={}_alpha={}_seed={}{}".format(args.data_id,
                                                           args.k,
                                                           args.alpha,
                                                           args.seed,
                                                           '_batch={}'.format(args.batch_col) if args.batch_col else '')
    DA_file = osp.join(args.outdir, prefix + '.DAresults.cna.csv')
    DA_df.to_csv(DA_file, index=True)


if __name__ == "__main__":
    args = parse_args()
    runDA(args)