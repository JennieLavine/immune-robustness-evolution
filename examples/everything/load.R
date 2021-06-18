#!/usr/bin/env Rscript

# This script loads the data from the output of `everything.jl` and reconstructs
# matrices from output data

library(RSQLite)

db <- dbConnect(SQLite(), "output.sqlite")

metaTable <- dbReadTable(db, "meta")
edgesTable <- dbReadTable(db, "edges")
summaryTable <- dbReadTable(db, "summary")

# Initialize S0, C0
n_receptors <- metaTable$n_receptors
n_effectors <- metaTable$n_effectors
S0 <- matrix(FALSE, nrow = n_receptors, ncol = n_effectors)
C0 <- matrix(0.0, nrow = n_receptors, ncol = n_effectors)
for(i in 1:nrow(edgesTable)) {
  stopifnot(edgesTable$id[i] == 0)
  S0[edgesTable$receptor[i], edgesTable$effector[i]] = TRUE
  C0[edgesTable$receptor[i], edgesTable$effector[i]] = edgesTable$weight[i]
}

# Reconstruct a list mapping matrix ID to (matrix, fitness) pairs
all_matrices <- list()
for(i in 1:nrow(summaryTable)) {
  id <- summaryTable$id[i]
  cat(sprintf("Processing id %d...\n", id))
  if(id == 0) {
    all_matrices[["0"]] <- list(S = S0, fitness = summaryTable$fitness[i])
  }
  else {
    parent_id <- summaryTable$parent_id[i]
    S <- all_matrices[[as.character(parent_id)]]$S
    S[summaryTable$removed_receptor[i], summaryTable$removed_effector[i]] <- FALSE
    fitness <- summaryTable$fitness[i]
    all_matrices[[as.character(id)]] <- list(S = S, fitness = fitness)
  }
}

# Do something with all this
