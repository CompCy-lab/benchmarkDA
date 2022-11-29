### METHODS ###

suppressPackageStartupMessages(
    {
        library(SingleCellExperiment)
        library(scran)
        library(DAseq)
        library(miloR)
        library(tibble)
        library(dplyr)
        library(igraph)
        library(cydar)
        library(pdist)
        library(reshape2)
    }
)


## Milo method

run_milo <- function(
    sce,
    condition_col,
    sample_col,
    reduced.dim="PCA",
    k=15,
    d=30,
    prop=0.1,
    returnMilo = TRUE,
    batch_col=NULL
){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Build graph neighbourhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim)
  ## Test DA
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  ## Note: Milo uses the standard R design matrix to label the conditions.
  DA_results <- testNhoods(milo, design = design, design.df = design_df, reduced.dim=reduced.dim)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results
  }
}

milo2output <- function(
    milo,
    da_res,
    out_type="continuous",
    alphas=NULL
){
  if (out_type=="continuous") { 
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[, 1]
  } else {
    da.cell <- vector(mode="list", length=length(alphas))
    for (i in seq_len(length(alphas))) {
      da.nhoods <- ifelse(da_res$SpatialFDR < alphas[i], ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
      da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
      da.cell.mat <- milo@nhoods %*% da.nhoods.mat
      da.cell[[i]] <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
    }
  }
  da.cell
}


## DA-seq method

run_daseq <- function(
    sce,
    k.vec,
    condition_col,
    sample_col="synth_samples",
    reduced.dim = "PCA",
    d=30
){
  # prepare the params needed by func: getDAcells
  condition_vec <- colData(sce)[[condition_col]]
  conds <- levels(factor(condition_vec))
  cell.labels <- sce[[sample_col]]
  if (length(conds) != 2){
    stop(str_c("daseq can only handle binary labels, but got ", length(conds)))
  }
  # main routine to get DA cells
  daseq_res <- getDAcells(X=reducedDim(sce, reduced.dim)[,1:d],
                          cell.labels=cell.labels,
                          labels.1=unique(cell.labels[condition_vec == conds[1]]),
                          labels.2=unique(cell.labels[condition_vec == conds[2]]),
                          k.vector=k.vec,
                          size=1,
                          do.plot=FALSE)
  return(daseq_res)
}

daseq2output <- function(
    sce,
    daseq_res,
    sample_col="synth_samples",
    out_type="continuous",
    pred.thres=NULL
){
  if (out_type=="continuous") {
    da.cell <- daseq_res$da.pred
  } else {
    if (is.null(pred.thres)) {
      stop(str_c("when out_type == ", out_type, ", pred.thres must be provided."))
    }
    cell.labels <- sce[[sample_col]]
    da.cell <- vector(mode="list", length=length(pred.thres))
    for (i in seq_len(length(pred.thres))) {
      da.cell[[i]] <- rep("NotDA", length(cell.labels))
      res <- updateDAcells(X=daseq_res, pred.thres=c(-pred.thres[i], pred.thres[i]), force.thres = TRUE)
      da.cell[[i]][res$da.up] <- "PosLFC"
      da.cell[[i]][res$da.down] <- "NegLFC"
    }
  }
  da.cell
}


## Cydar method

run_cydar <- function(
    sce,
    condition_col="synth_labels",
    sample_col="synth_samples",
    reduced.dim="pca.corrected",
    d=30,
    batch_col=NULL,
    tol=0.5,
    downsample=10,
    returnCd=TRUE
){
  ## Make the design matrix as Milo
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
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}

cydar2output <- function(
    sce,
    cd,
    da_res,
    out_type="continuous",
    alphas=NULL
){
  nhs <- lapply(cellAssignments(cd), function(hs) as.vector(hs))
  ## Recover cell ids
  ordered.cells <- colnames(cellIntensities(cd))
  hs_mat <- sapply(nhs, function(nh) ifelse(seq_along(ordered.cells) %in% nh, 1, 0))
  rownames(hs_mat) <- ordered.cells
  colnames(hs_mat) <- rownames(da_res)
  if (out_type=="continuous") { 
    da.cell.mat <- hs_mat %*% da_res$logFC
    da.cell <- da.cell.mat[, 1]
    da.cell <- da.cell[colnames(sce)]
  } else {
    da.cell <- vector(mode="list", length(alphas))
    for (i in seq_len(length(alphas))){
      da.hs <- ifelse(da_res$SpatialFDR < alphas[i], ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
      da.hs.mat <- sapply(unique(da.hs), function(x) as.numeric(da.hs==x))
      da.cell.mat <- hs_mat %*% da.hs.mat
      da.cell[[i]] <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
      da.cell[[i]] <- da.cell[[i]][colnames(sce)]
    }
  }
  da.cell
}

run_louvain_nbglm <- function(
    sce,
    condition_col,
    sample_col,
    reduced.dim="PCA",
    d=30,
    k=15,
    resolution=1,
    batch_col=NULL,
    norm.method="TMM"
){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    dplyr::rename(sample=sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Louvain clustering
  X_red_dim = reducedDim(sce, reduced.dim)[,1:d]
  sce.graph <- buildKNNGraph(t(X_red_dim), k=k)
  louvain.clust <- cluster_louvain(sce.graph, resolution=resolution)
  message(str_c("#cluster: ", length(louvain.clust)))
  louvain.clust.ids <- membership(louvain.clust)
  
  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  clust.df <- data.frame("cell_id"=colnames(sce), "Louvain.Clust"=as.character(louvain.clust.ids))
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec
  
  louvain.count <- table(clust.df$Louvain.Clust, clust.df$Sample)
  attributes(louvain.count)$class <- "matrix"
  
  ## Test with same NB-GLM model as the Milo
  if(norm.method %in% c("TMM")){
    message("Using TMM normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("logMS")){
    message("Using logMS normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
  }
  
  model <- model.matrix(design, data=design_df)
  rownames(model) <- design_df$sample
  model <- model[colnames(louvain.count), ]
  
  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust=TRUE)
  n.coef <- ncol(model)
  louvain.res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
  
  clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, 'logFC']
  clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, 'FDR']
  return(clust.df)
}


louvain2output <- function(
    louvain_res,
    out_type="continuous",
    alphas=NULL,
    lfc_threshold=0
){
  if (out_type=="continuous") {
    da.cell <- louvain_res$logFC
  } else {
    da.cell <- vector(mode="list", length=length(alphas))
    for (i in seq_len(length(alphas))){
      da.cell[[i]] <- ifelse(louvain_res$FDR < alphas[i] & abs(louvain_res$logFC) > lfc_threshold, ifelse(louvain_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    }
  }
  da.cell
}

### RUN BENCHMARK ###

runDA <- function(
    sce,
    coldata,
    X_pca, 
    method, 
    condition_col='synth_labels', 
    sample_col="synth_samples",
    params=list(milo=list(k=30),
                daseq=list(k.vec=c(10,20,30,40)),
                louvain=list(k=30, resolution=1),
                cydar=list(tol=0.5, downsample=3)),
    d=30,
    out_type="label",
    seed=42
){
  ## Check that method name is in params
  set.seed(seed)
  if (!method %in% names(params)) {
    stop(paste("Specify parameters for method", method))
  }
  ## Check valid method
  if (!method %in% c("milo", "milo_batch", "louvain", 'daseq', "louvain_batch", 'cydar', 'cydar_batch')) {
    stop("Unrecognized method")
  }
  
  ## Add reduced dim + coldata to sce
  colData(sce) <- DataFrame(coldata)
  reducedDim(sce, "pca_batch") <- as.matrix(X_pca)
  
  ## Select the method and run DA testing
  if (method=="milo") {
    start.time <- Sys.time()
    milo_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                         reduced.dim = "pca_batch", d=d, k=params$milo$k)
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(milo_res$DAres$SpatialFDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- milo2output(milo_res$Milo, milo_res$DAres, out_type=out_type, alphas=alphas)
    if (out_type == "continuous"){
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length=length(da.cell))
      for (i in seq_len(length(da.cell))){
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "milo_batch") {
    start.time <- Sys.time()
    milo_batch_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                               reduced.dim = "pca_batch", d=d, k=params$milo$k, batch_col="synth_batches")
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(milo_batch_res$DAres$SpatialFDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- milo2output(milo_batch_res$Milo, milo_batch_res$DAres, out_type=out_type, alphas=alphas)
    if (out_type == "continuous"){
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(da.cell))
      for (i in seq_len(length(da.cell))){
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "cydar") {
    start.time <- Sys.time()
    cydar_res <- run_cydar(sce, condition_col=condition_col, sample_col=sample_col,
                           reduced.dim = "pca_batch", d=d, tol=params$cydar$tol,
                           downsample=params$cydar$downsample)
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(cydar_res$DAres$SpatialFDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- cydar2output(sce, cydar_res$Cd, cydar_res$DAres, out_type=out_type, alphas=alphas)
    if (out_type == "continuous") {
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(alphas))
      for (i in seq_len(length(da.cell))) {
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "cydar_batch") {
    start.time <- Sys.time()
    cydar_batch_res <- run_cydar(sce, condition_col=condition_col, sample_col=sample_col,
                                 reduced.dim = "pca_batch", d=d, tol=params$cydar$tol,
                                 downsample=params$cydar$downsample, batch_col="synth_batches")
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(cydar_batch_res$DAres$SpatialFDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- cydar2output(sce, cydar_batch_res$Cd, cydar_batch_res$DAres, out_type=out_type, alphas=alphas)
    if (out_type == "continuous") {
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(alphas))
      for (i in seq_len(length(da.cell))) {
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "daseq") {
    start.time <- Sys.time()
    daseq_res <- run_daseq(sce, k.vec=params$daseq$k.vec, condition_col=condition_col, 
                           reduced.dim = "pca_batch", d=d)
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different thresholds
    pred_thresholds <- quantile(abs(daseq_res$da.pred), seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- daseq2output(sce, daseq_res, out_type=out_type, pred.thres=pred_thresholds)
    if (out_type == "continuous"){
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(da.cell))
      for (i in seq_len(length(da.cell))){
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$thres <- pred_thresholds[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "louvain") {
    ## Run louvain
    start.time <- Sys.time()
    louvain_res <- run_louvain_nbglm(sce, condition_col=condition_col, sample_col=sample_col, reduced.dim = "pca_batch",
                                     d=d, k=params$louvain$k, resolution=params$louvain$resolution)
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(louvain_res$FDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- louvain2output(louvain_res, out_type=out_type, alphas=alphas)
    if (out_type == "continuous"){
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(da.cell))
      for (i in seq_len(length(da.cell))){
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  } else if (method == "louvain_batch") {
    start.time <- Sys.time()
    louvain_batch_res <- run_louvain_nbglm(sce, condition_col=condition_col, sample_col=sample_col, reduced.dim="pca_batch",
                                           d=d, k=params$louvain$k, resolution=params$louvain$resolution, batch_col="synth_batches")
    run_time <- as.numeric(Sys.time() - start.time, units="secs")
    # get da.cells with different alphas
    alphas <- quantile(louvain_batch_res$FDR, seq(1e-8, 1 - 1e-8, 0.01), names=FALSE)
    da.cell <- louvain2output(louvain_batch_res, out_type=out_type, alphas=alphas)
    if (out_type == "continuous"){
      bm.out <- out2bm(da.cell, sce, method, run_time)
    } else {
      bm.out <- vector(mode="list", length = length(da.cell))
      for (i in seq_len(length(da.cell))) {
        bm.out[[i]] <- calculate_outcome(out2bm(da.cell[[i]], sce, method))
        bm.out[[i]]$alpha <- alphas[i]
        bm.out[[i]]$runtime <- run_time
      }
    }
  }
  
  ## Save results
  if (out_type != "continuous") {
    benchmark_res <- bind_rows(bm.out)
  } else {
    benchmark_res <- bm.out
  }
  return(benchmark_res)
}


out2bm <- function(out, sce, method, runtime=NULL) {
  bm <- data.frame(bm_out=out)
  if (length(out) < ncol(sce)) {
    return(out)
  }
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  if (!is.null(sce$true_DA_clust)) {
    bm$true_clust <- sce$true_DA_clust
  }
  long_bm <- pivot_longer(bm, cols = bm_out, names_to='method', values_to="pred")
  long_bm[["method"]] <- method
  if (!is.null(runtime)) {
    long_bm[['runtime']] <- runtime
  }
  return(long_bm)
}


calculate_outcome <- function(long_bm){
  long_bm <- long_bm %>%
    mutate(outcome=case_when(true==pred & pred!="NotDA" ~ 'TP',
                             true!=pred & pred!="NotDA" ~ 'FP',
                             true!=pred & pred=="NotDA" ~ 'FN',
                             true==pred & pred=="NotDA" ~ "TN"
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
    mutate(TPR=TP/(TP+FN), FPR=FP/(FP+TN), TNR=TN/(TN+FP), FNR = FN/(FN+TP),
           FDR = FP/(TP+FP),
           Precision = TP/(TP+FP),
           Power = 1 - FNR,
           Accuracy = (TP + TN)/(TP + TN + FP + FN)
    )
}
