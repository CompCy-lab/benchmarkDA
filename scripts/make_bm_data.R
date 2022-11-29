### Benchmarking methods with synthetic labels and batch effects on real data

suppressPackageStartupMessages(
    {
        library(argparse)
        library(tidyverse)
        library(SingleCellExperiment)
        library(scran)
    }
)

source('./data_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character", help="path to RDS storing SingleCellExperiment object")
parser$add_argument("seed", type="integer", help="random seed")
parser$add_argument("population", type="character", help="Cell type of DA")
parser$add_argument("--pop_enrichment", type="double", default=0.75, help="Max condition probability in DA population")
parser$add_argument("--pop_col", type="character", default="celltype", help="the column name to specify population")
parser$add_argument("--reduced.dim", type="character", default="PCA", help="the name of reduced dim to use")
parser$add_argument("--k", type="integer", default=50, help="KNN parameter")
parser$add_argument("--data_id", type="character", default="linear", help="ID for the dataset used")
parser$add_argument("--make_batch_effect", type="character", default="yes", help="should synthetic batch effects be added? (yes/no)")
parser$add_argument("--outdir", type="character", help = "path to the output directory")
args <- parser$parse_args()

# read the parameters
data_path <- args$data_RDS
seed <- args$seed
pop <- args$population

pop_enr <- args$pop_enrichment
pop_col <- args$pop_col
reduced.dim <- args$reduced.dim
k <- args$k
data_id <- args$data_id
make_batch_effects <- args$make_batch_effect
outdir <- args$outdir

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

# replace "_" with blank
if (str_detect(pop, "_")) {
  pop <- str_replace(pop, "_", " ")
}

# Simulate the dataset
if (data_id == "cluster") {
  sce <- add_synthetic_labels_by_cluster(
    sce, pop=pop, pop_column=pop_col, seed=seed,
    pop_enr=pop_enr, m=2)
} else {
  sce <- add_synthetic_labels_pop(
    sce, pop=pop, pop_column=pop_col, seed=seed,
    pop_enr=pop_enr, redDim=reduced.dim, m=2
  )
}

# Set the threshold to classify NegLFC, PosLFC and NotDA
if (data_id == "cluster") {
  if (pop_enr < 0.5) {
    da_lower <- pop_enr + (pop_enr / 100) * 10
    da_upper <- 1 - da_lower
  } else {
    da_upper <- pop_enr - (pop_enr / 100) * 10
    da_lower <- 1 - da_upper
  }
} else {
  # assign labels using the quantile of #target_pop / #entire_pop
  pop_tbl <- table(sce[[pop_col]])
  if (pop_enr < 0.5) {
    da_lower <- quantile(sce$Condition1_prob, pop_tbl[pop] / sum(pop_tbl), names=FALSE)
    da_upper <- 1 - da_lower
  } else {
    da_upper <- quantile(sce$Condition1_prob, 1 - pop_tbl[pop] / sum(pop_tbl), names=FALSE)
    da_lower <- 1 - da_upper
  }
}

# make sure da_upper > da_lower
stopifnot(da_upper > da_lower)

# assign the synthetic labels
true_labels <- ifelse(sce$Condition2_prob < da_lower, "NegLFC", ifelse(sce$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(sce)[["true_labels"]] <- true_labels
if (str_detect(pop, " ")) {
  pop <- str_replace(pop, " ", "_")
}

## Save coldata ##
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)

message(str_c("Writing ", outprefix, " ..."))
coldata <- data.frame(colData(sce)) %>% rownames_to_column()
write_csv(coldata, str_c(outdir, outprefix, ".coldata.csv"))

## Simulate batch effects of different magnitude if needed
set.seed(seed)
if (make_batch_effects == "yes") {
  message("Simulating batch effects...")
  bm_sce_ls <- lapply(c(0, 0.75, 1, 1.25, 1.5), function(sd){
    sce_be <- add_batch_effect(sce, batch_col="synth_batches", norm_sd=sd)
    X_pca <- reducedDim(sce_be, "pca_batch")
    ## Save reduced dims
    write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(outdir, outprefix, "_batchEffect", sd, ".pca.csv"))
  })
} else {
  X_pca <- reducedDim(sce, reduced.dim)
  ## Save reduced dims
  write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(outdir, outprefix, "_batchEffect0.pca.csv"))
}
