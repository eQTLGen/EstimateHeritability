#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(argparse)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Declare constants
  parser <- ArgumentParser(description='process delete values')
  # Basic LD Score Estimation Flags'
  # Filtering / Data Management for LD Score
  parser$add_argument('--h2', type="character",
                      help='matrix')
  parser$add_argument('--delete-vals', default=NULL, type="character", nargs="+")

  args <- parser$parse_args(argv)

  mean_hsq <- mean(fread(args$h2) %>% pull(h2_obs))

  delete_values <- do.call(cbind, lapply(args$delete_vals, function(path) {
    fread(path, header=F)$V1
  }))

  ensemble_ids <- str_extract(args$delete_vals, "ENSG\\d+")
  colnames(delete_values) <- ensemble_ids

  write.table(delete_values, "delete_values_combined.tsv", quote=F, sep="\t", row.names=F, col.names=T)

  mean_delete_values <- apply(delete_values, 1, mean)

  pseudovalues <- 200 * mean_hsq - 199 * mean_delete_values

  se <- sqrt(var(pseudovalues))
  print(se)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
