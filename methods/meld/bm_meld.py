import argparse
import time
import anndata
import numpy as np
import pandas as pd
import os.path as osp

from _meld import runMELD, meld2output


def parse_args():
    parser = argparse.ArgumentParser("benchmark the performance of meld method.")

    parser.add_argument("--data_dir", type=str, default="./data", help="output directory")
    parser.add_argument("--data_id", type=str, default="linear", help="the dataset id which specifies topology and reps")
    parser.add_argument("--k", type=int, default=30, help="K parameter to build KNN graph")
    parser.add_argument("--beta", type=float, default=60, help="the beta parameter in Meld")
    parser.add_argument("--pop", type=str, default="M1", help="the population number")
    parser.add_argument("--pop_enr", type=str, default=0.75, help="population enrichment")
    parser.add_argument("--be_sd", type=float, default=0, help="SD of batch effect")
    parser.add_argument("--seed", type=int, default=43, help="random seed")
    parser.add_argument("--out_type", type=str, default="label", help="output data format")
    parser.add_argument("--outdir", type=str, help="dir to save benchmark results")    

    args = parser.parse_args()
    return args


def csv2anndata(data_path: str, obs_path: str):
    # read the csv files
    X_pca = pd.read_csv(data_path, index_col=0)
    obs = pd.read_csv(obs_path, index_col=0)
    assert (X_pca.index == obs.index).all(), \
        "the index of data and meta data is different from each other"
    # create the anndata
    adata = anndata.AnnData(X_pca, obs=obs, dtype=np.float64)
    adata.obs.index.name = "cell"
    return adata
    

def calculate_outcome(bm: pd.DataFrame):
    # calculate the metrics, TP, FP, TN, FN and so on.
    TP = ((bm['true'] == bm['pred']) & (bm['pred'] != 'NotDA')).sum()
    FP = ((bm['true'] != bm['pred']) & (bm['pred'] != 'NotDA')).sum()
    FN = ((bm['true'] != bm['pred']) & (bm['pred'] == 'NotDA')).sum()
    TN = ((bm['true'] == bm['pred']) & (bm['pred'] == 'NotDA')).sum()
    
    metrics_dic = {
        'TP': TP, 'FP': FP,
        'FN': FN, 'TN': TN,
        'TPR': [TP / (TP + FN)],
        'FPR': [FP / (FP + TN)],
        'TNR': [TN / (TN + FP)],
        'FNR': [FN / (FN + TP)],
        'FDR': [FP / (TP + FP)],
        'Precision': [TP / (TP + FP)],
        'Power': [1 - FN / (FN + TP)],
        'Accuracy': [(TP + TN) / (TP + TN + FP + FN)]
    }
    
    return pd.DataFrame(metrics_dic, index=['metric'])


def out2bm(out, adata, runtime=None):
    # prepare the data frame as output
    bm = pd.DataFrame()
    bm['true_prob'] = adata.obs['Condition2_prob']
    bm['true'] = adata.obs['true_labels']
    bm['pred'] = out
    bm['method'] = 'meld'
    if runtime:
        bm['runtime'] = runtime
    
    return bm


def runDA(args):
    # set random seed for reproducibility
    np.random.seed(args.seed)
    # read the data
    data_dir = args.data_dir
    prefix = "benchmark_{}_pop_{}_enr{}_seed{}".format(args.data_id,
                                                        args.pop,
                                                        args.pop_enr,
                                                        args.seed)
    pcaData = osp.join(data_dir, prefix + "_batchEffect{}.pca.csv".format(
        int(args.be_sd) if args.be_sd.is_integer() else args.be_sd))
    colData = osp.join(data_dir, prefix + ".coldata.csv")
    adata = csv2anndata(data_path=pcaData, obs_path=colData)
    
    # run analysis using MELD method
    start_time = time.time()
    meld_res = runMELD(adata, k=args.k, sample_col="synth_samples", label_col="synth_labels", beta=args.beta)
    run_time = time.time() - start_time
    # get da cells with different thresholds
    lower = 1/len(adata.obs['synth_labels'].unique()) + 1e-8
    upper = meld_res.max() - 1e-8
    assert (lower < upper), \
        "lower bound is not less than upper bound"
    thresholds = np.arange(lower, upper, (upper - lower) / 100)
    da_cell = meld2output(meld_res, out_type=args.out_type, thresholds=thresholds) # out: [np.array([str])]
    
    # prepare benchmark data frame
    if args.out_type == "continuous":
        bm_out = out2bm(da_cell, adata, run_time)
    else:
        bm_out = []
        for i, da in enumerate(da_cell):
            bm_out.append(calculate_outcome(out2bm(da, adata)))
            bm_out[i][['method', 'thres', 'runtime']] = ['meld', thresholds[i], run_time]
            
    if args.out_type != "continuous":
        bm_out = pd.concat(bm_out, axis=0)
        # print("AUC: ", metrics.auc(bm_out['FPR'].values, bm_out['TPR'].values))
    
    # save the benchmark result
    bm_resfile = osp.join(args.outdir, prefix + "_batchEffect{}.DAresults.meld".format(
        int(args.be_sd) if args.be_sd.is_integer() else args.be_sd) + ".csv")
    bm_out.to_csv(bm_resfile, index=False)


if __name__ == "__main__":
    args = parse_args()
    runDA(args)