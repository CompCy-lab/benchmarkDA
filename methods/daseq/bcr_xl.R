rm(list=ls())
setwd("~/Documents/proj/benchmark_scCut/daseq/")

library(DAseq)
python2use <- "/usr/local/anaconda3/envs/daseq/bin/python"

# prepare the data for DAseq
X <- read.csv("../data/bcr-csv/X.csv", header = FALSE)
obs <- read.csv("../data/bcr-csv/obs.csv")
embedding <- read.csv("../data/bcr-csv/obsm.csv")

rownames(X) <- obs[,"cell"]
rownames(embedding) <- obs[,"cell"]

# set positive and negative labels
labels_pos <- c("BCR-XL")
labels_neg <- c("Ref")

# Get DA cells
da_cells <- getDAcells(
  X = X,
  cell.labels = obs[,"phenotype"],
  labels.1 = labels_neg,
  labels.2 = labels_pos,
  k.vector = seq(10, 30, 5),
  plot.embedding = embedding)

da_cells$pred.plot

# 


