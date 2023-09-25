#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)

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

  gene_id <- argv[1]
  ldsc_log <- argv[2]

  # Process input
  input_table <- as_tibble(fread(ldsc_log, sep=":", header=F, nrows = 5))

  # Perform method
  table_processed <- input_table %>%
    separate(V2, c("estimate", "stderr"), sep = " ") %>%
    mutate(stderr = as.double(str_extract(stderr, regex("\\d+\\.?\\d*")))) %>%
    rename(c("name" = "V1")) %>%
    mutate(estimate = as.double(estimate),
           gene = gene_id)

  # Process output
  write.table(table_processed, "ldsc_matrix.txt", col.names = F, row.names = F, sep = "\t", quote = F)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}