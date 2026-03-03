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
  parser$add_argument('--out', type="character",
                      help='output')
  parser$add_argument('--genes', type="character",
                      nargs="+")
  parser$add_argument('--delete-vals', default=NULL, type="character", nargs="+")

  args <- parser$parse_args(argv)

  delete_values <- do.call(cbind, lapply(args$delete_vals, function(path) {
    fread(path, header=F)$V1
  }))

  #ensemble_ids <- str_extract(args$delete_vals, "ENSG\\d+")
  colnames(delete_values) <- args$genes

  write.table(delete_values, args$out, quote=F, sep="\t", row.names=F, col.names=T)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
