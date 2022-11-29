rm(list=ls())
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DAseq)
  library(miloR)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(igraph)
  library(cydar)
  library(pdist)
  library(reshape2)
})


run_cydar <- function(sce, condition_col="synth_labels",
                      sample_col="synth_samples",
                      reduced.dim="pca.corrected",
                      d=20,
                      batch_col=NULL,
                      alpha=0.1,
                      tol=1.0,
                      downsample=10,
                      returnCd=TRUE){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Make list for each sample
  sample_ls <- split(1:ncol(sce), sce[[sample_col]])
  processed.exprs <- lapply(sample_ls, function(s) reducedDim(sce[,s], reduced.dim)[,1:d])
  cd <- prepareCellData(processed.exprs)
  ## Count cells in hyperspheres
  cd <- cydar::countCells(cd, tol=tol, filter=1, downsample=downsample)
  # do DA testing with edgeR
  cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
  
  sim.design <- model.matrix(design, data=design_df)[colnames(cd.dge),]
  sim.dge <- estimateDisp(cd.dge, sim.design)
  sim.fit <- glmQLFit(sim.dge, sim.design)
  sim.res <- glmQLFTest(sim.fit, coef=2)
  
  # control the spatial FDR
  cydar.res <- sim.res$table
  cydar.res$SpatialFDR <- spatialFDR(intensities(cd), sim.res$table$PValue)
  is.sig <- cydar.res$SpatialFDR <= alpha
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}

cydar2output <- function(sce, cd, da_res, out_type="continuous", alpha=0.1, sample_col="synth_samples", reduced.dim="pca.corrected", d=30){
  nhs <- lapply(cellAssignments(cd), function(hs) as.vector(hs))
  # hs_mat <- sapply(nhs, function(nh) ifelse(1:sum(cd@colData$totals) %in% nh, 1, 0))
  ## Recover cell ids
  ordered.cells <- colnames(cellIntensities(cd))
  hs_mat <- sapply(nhs, function(nh) ifelse(seq_along(ordered.cells) %in% nh, 1, 0))
  rownames(hs_mat) <- ordered.cells
  colnames(hs_mat) <- rownames(da_res)
  if (out_type=="continuous") { 
    da.cell.mat <- hs_mat %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.hs <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.hs.mat <- sapply(unique(da.hs), function(x) as.numeric(da.hs==x))
    da.cell.mat <- hs_mat %*% da.hs.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell[colnames(sce)]
}

calculate_outcome <- function(long_bm){
  long_bm <- long_bm %>%
    mutate(outcome=case_when(true==pred & pred!="NotDA" ~ 'TP',
                             true!=pred & pred!="NotDA" ~ 'FP',
                             true!=pred & pred=="NotDA" ~ 'FN',
                             true==pred & pred=="NotDA" ~ 'TN'
    )) %>%
    group_by(method, outcome) %>%
    summarise(n=n()) %>%
    pivot_wider(id_cols=method, names_from=outcome, values_from=n, values_fill=0)
  
  check_cols <- c("TP","FP","FN","TN") %in% colnames(long_bm)
  if (any(!check_cols)) {
    add_cols <- c("TP","FP","FN","TN")[!check_cols]
    for (col in add_cols) {
      long_bm[[col]] <- rep(0, nrow(long_bm))
    }
  }

  long_bm %>%
    mutate(TPR = TP / (TP + FN),
           FPR = FP / (FP + TN),
           TNR = TN / (TN + FP),
           FNR = FN / (FN + TP),
           FDR = FP / (TP + FP),
           Precision = TP / (TP + FP),
           Power = 1 - FNR,
           Accuracy = (TP + TN) / (TP + TN + FP + FN)
    )
}


# main
sce <- readRDS("/Users/mac/Documents/proj/benchmarkDA/data/synthetic/branch/branch_data_bm.RDS")
X_pca <- read.csv("/Users/mac/Documents/proj/benchmarkDA/data/synthetic/branch/benchmark_branch_pop_M1_enr0.75_seed43_batchEffect0.pca.csv") %>% column_to_rownames()
coldata <- read.csv("/Users/mac/Documents/proj/benchmarkDA/data/synthetic/branch/benchmark_branch_pop_M1_enr0.75_seed43.coldata.csv") %>% column_to_rownames()
## Add reduced dim + coldata to sce
colData(sce) <- DataFrame(coldata)
reducedDim(sce, "pca_batch") <- as.matrix(X_pca)

# run the milo method
cydar_res <- run_cydar(sce, condition_col="synth_labels", sample_col="synth_samples",
                       reduced.dim = "pca_batch", d=50, tol=2, downsample=5)
out <- cydar2output(sce, cydar_res$Cd, cydar_res$DAres, out_type = "label", alpha=0.05)

# post-process outputs
bm <- data.frame(bm_out=out)
bm$true_prob <- sce$Condition2_prob
bm$true <- sce$true_labels
long_bm <- pivot_longer(bm, cols = bm_out, names_to='cydar', values_to="pred")
long_bm[["method"]] <- 'cydar'

# compute metrics
metrics <- calculate_outcome(long_bm)
print(metrics)