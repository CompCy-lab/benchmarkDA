### Run DA methods in R ###

suppressPackageStartupMessages(
    {
        library(argparse)
        library(tidyverse)
        library(SingleCellExperiment)
        library(scran)
    }
)

source('./benchmark_utils.R')
options(dplyr.summarise.inform = FALSE)

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character", help="Path to RDS storing SingleCellExperiment object")
parser$add_argument("method", type="character", help="DA method to use")
parser$add_argument("seed", type="integer", help="Random seed")
parser$add_argument("population", type="character", help="Cell type of DA")
parser$add_argument("--pop_enrichment", type="double", default=0.85, help="Max condition probability in DA population")
parser$add_argument("--batchEffect_sd", type="double", default=0, help="Standard deviation of batch effects")
parser$add_argument("--k", type="integer", default=20, help = "KNN parameter")
parser$add_argument("--resolution", type="double", default=1, help="Resolution of Louvain clustering")
parser$add_argument("--tol", type="double", default=NULL, help="Scalar proportional to the hypersphere radius in Cydar")
parser$add_argument("--downsample", type="integer", default=3, help="Downsampling ratio of cydar")
parser$add_argument("--step", type="integer", default=50, help="the step size of k in daseq")
parser$add_argument("--data_id", type="character", default="linear", help="ID for the dataset used")
parser$add_argument("--data_dir", type="character", help="path to the input data directory")
parser$add_argument("--outdir", type="character", help="path to the output data directory")
args <- parser$parse_args()


### Assign all the parameters ###
data_path <- args$data_RDS
DA_method <- args$method
seed <- args$seed
pop <- args$population

pop_enr <- args$pop_enrichment
be_sd <- args$batchEffect_sd
k <- args$k
resolution <- args$resolution
tol <- args$tol
downsample <- args$downsample
step <- args$step
data_id <- args$data_id
data_dir <- args$data_dir
bm_outdir <- args$outdir


## Load RDS data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- read_csv(paste0(data_dir, outprefix, ".coldata.csv")) %>% column_to_rownames()
X_pca <- read_csv(str_c(data_dir, outprefix, "_batchEffect", be_sd, ".pca.csv")) %>% column_to_rownames()  

## cydar radius scaler picked w/ heuristic
tol_dataset <- list(cluster=2.8, cluster_balanced=2.85, branch=2.4, linear=2.3, 'covid19-pbmc'=2.1, 'bcr-xl'=0.75,
                    test_scale_4000=1.2, test_scale_10000=1.2, test_scale_15000=1.2, test_scale_30000=1.2,
                    test_scale_50000=1.2, test_scale_100000=1.2, pancreas=2.85, levine32=0.45)
if (is.null(tol)) {
  tol <- 0.5
  if (!is.null(tol_dataset[[data_id]])) {
    tol <- tol_dataset[[data_id]]
  }
}

## set the parameters for all the methods
benchmark_params = list(
  milo = list(k=k),
  milo_batch = list(k=k),
  meld = list(k=k),
  daseq = list(k.vec=seq(k, 500, step)),
  louvain = list(k=k, resolution=resolution),
  louvain_batch = list(k=k, resolution=resolution),
  cydar = list(tol=tol, downsample=downsample),
  cydar_batch = list(tol=tol, downsample=downsample),
  monocle3 = list(k=k, resolution=resolution),
  sclca = list(clust.max=15)
)

## Run DA method ##
out <- runDA(sce, X_pca, coldata=coldata, method=DA_method, params=benchmark_params, d=ncol(X_pca))

## Save results ##
write_csv(out, str_c(bm_outdir, outprefix, "_batchEffect", be_sd, ".DAresults.", DA_method, ".csv"))
