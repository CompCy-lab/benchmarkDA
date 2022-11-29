### SYNTHETIC LABELS ###

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

.find_centroid <- function(X_emb, cluster_membership){
  cl.ixs <- split(1:nrow(X_emb), cluster_membership)  
  centroid_emb <- sapply(cl.ixs, function(x) colMeans(X_emb[x, , drop=FALSE]))
  centroid_emb
}

.member_weight <- function(x, centroid_dist, m=2){
  # centroid_dist <- pdist(t(x), t(centroid_emb))@dist
  w_memb <- sapply(centroid_dist, function(x) 1/sum(x/centroid_dist)^(2/(m-1)))
}

.scale_to_range <- function(x, min=1, max=10){
  ((x - min(x))/(max(x)-min(x)))*(max-min) + min
}


.logit <- function(x, a=1){
  1/(1+exp(- a * x))
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels_pop <- function(sce, # SingleCellExperiment obj
                                     pop, pop_column="celltype",
                                     pop_enr = 0.7,
                                     redDim='pca.corrected', # embedding to use to simulate differential abundance
                                     n_conditions=2, # number of conditions to simulate
                                     n_replicates=3, # number of replicates per condition
                                     n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                     condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                     m=2, # Fuzziness parameter (higher m, more fuzziness)
                                     a_logit=0.5, # logit parameter
                                     cap_enr=NULL,
                                     seed=42){
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  X_emb = reducedDim(sce, redDim)
  
  ## Find cluster center
  cluster_membership = sce[[pop_column]]
  centroid_emb <- .find_centroid(X_emb, cluster_membership)

  ## Assign weight to each cell for each cluster center
  centroid_dist <- pdist(X_emb, t(centroid_emb))
  centroid_dist <- as.matrix(centroid_dist)
  
  w <- sapply(1:ncol(centroid_dist),  function(j) 
    sapply(1:nrow(centroid_dist), function(i) 
      1/sum(centroid_dist[i,j]/centroid_dist[i,])^(2/(m-1))
    ) 
  )
  colnames(w) <- colnames(centroid_emb)
  rownames(w) <- rownames(X_emb)
  w <- apply(scale(w), 2, .logit, a=a_logit)
  ## Normalize weights from enr_score to 0.5
  enr_scores <- rep(0.5, ncol(w)) ## Generate enrichment prob for each cluster
  # enr_scores <- runif(ncol(w)) ## Generate enrichment prob for each cluster
  names(enr_scores) <- colnames(w)
  if(length(pop_enr) == length(pop)){
    enr_scores[pop] <- pop_enr
  } else{
    # assume all pops have the same enrichment
    pop_enr <- rep(pop_enr, length(pop))
    enr_scores[pop] <- pop_enr
  }
  
  
  # altering the baseline probability can induce a skew towards a condition across _all_ cells
  enr_prob <- sapply(1:ncol(w), function(i) .scale_to_range(w[,i], min=0.5*condition_balance,
                                                            max=enr_scores[i]))
  colnames(enr_prob) <- colnames(centroid_emb)
  
  # need to integrate over these to get the condition probabilities
  # need to set relevant pops only, force the others to ~0.5
  prob_matrix <- enr_prob[,pop]
  if(is(prob_matrix, "matrix")){
    cond_probability <- rowMeans(prob_matrix)
    for(x in seq_along(pop)){
      cond_probability[sce[[pop_column]] == pop[x]] <- prob_matrix[sce[[pop_column]] == pop[x], pop[x]]
    }
  } else{
    cond_probability <- prob_matrix
  }
  
  ## Cap probabilities (to have the same number of DA cells w different maximum Fold Change)
  if (!is.null(cap_enr)) {
    cond_probability <- ifelse(cond_probability > cap_enr, cap_enr, cond_probability)
  }
  # sim3$Condition1_prob <- ifelse(sim3$Condition1_prob > 0.8, 0.8, sim3$Condition1_prob)
  
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
   names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

## Simple synthetic condition labelling based on cluster membership
add_synthetic_labels_by_cluster <- function(sce, # SingleCellExperiment obj
                                            pop, pop_column="celltype",
                                            pop_enr = 0.7,
                                            # redDim='pca.corrected', # embedding to use to simulate differential abundance
                                            n_conditions=2, # number of conditions to simulate
                                            n_replicates=3, # number of replicates per condition
                                            n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                            condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                            m=2, # Fuzziness parameter (higher m, more fuzziness)
                                            a_logit=0.5, # logit parameter
                                            cap_enr=NULL,
                                            seed=42){
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  ## Set prop != 0.5 for cells in pop
  cond_probability <- rep(0.5, ncol(sce))
  cond_probability[sce[[pop_column]] == pop] <- pop_enr
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
    names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

## Group true DA cells in clusters (to test coverage of DA methods)
cluster_synthetic_labels <- function(embryo_sce, graph, min_cl_size=5){
  adj <- get.adjacency(graph)
  
  ## Retain DA cells
  da.adj <- adj[embryo_sce$true_labels!='NotDA',embryo_sce$true_labels!='NotDA']
  
  ## REmove edges between cells with discodant LFC sign
  da.cells.mat <- sapply(unique(embryo_sce$true_labels), function(x) as.numeric(embryo_sce$true_labels==x))
  da.cells.mat <- da.cells.mat[embryo_sce$true_labels!='NotDA',c("NegLFC", "PosLFC")]
  concord.sign.adj <- tcrossprod(da.cells.mat[,c("NegLFC", "PosLFC")], da.cells.mat[,c("NegLFC", "PosLFC")])
  concord.sign.adj <- as(concord.sign.adj, 'sparseMatrix')
  da.adj[concord.sign.adj == 0] <- 0
  
  ## Cluster DA cells
  da.graph <- graph_from_adjacency_matrix(da.adj, mode = 'undirected')
  clust <- igraph::cluster_louvain(da.graph)
  embryo_sce$true_DA_clust <- rep(NA, length(embryo_sce$true_labels))
  embryo_sce$true_DA_clust[embryo_sce$true_labels != "NotDA"] <- clust$membership
  
  ## Remove singletons (or less than min_cl_size cells)
  embryo_sce$true_DA_clust[embryo_sce$true_DA_clust %in% which(table(clust$membership) < min_cl_size)] <- NA
  
  embryo_sce
}

### SYNTHETIC BATCH EFFECT ###

add_batch_effect <- function(embryo_sce, batch_col="synth_samples", norm_sd=0.5){
  cellids_sample <- split(embryo_sce$cell_id, embryo_sce[[batch_col]])
  X_pca <- reducedDim(embryo_sce, "PCA")
  X_pca_batch <- X_pca

  for (b in names(cellids_sample)){
    batch_effect <- rnorm(ncol(X_pca), mean=0, sd = norm_sd)
    X_pca_batch[cellids_sample[[b]],] <- t(apply(X_pca_batch[cellids_sample[[b]],], 1, function(x) x + batch_effect))
  }
  
  reducedDim(embryo_sce, "pca_batch") <- X_pca_batch
  embryo_sce
}
