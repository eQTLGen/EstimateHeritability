#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)

# Declare constants

# Declare function definitions
read_heritability_table <- function(lines) {
  # Process input
  input_table <- as_tibble(fread(text=lines, sep=":", header=F, nrows = 5, fill = TRUE))

  # Perform method
  table_processed <- input_table %>%
    separate(V2, c("estimate", "stderr"), sep = " ") %>%
    mutate(stderr = as.double(str_extract(stderr, regex("\\d+\\.?\\d*")))) %>%
    rename(c("name" = "V1")) %>%
    mutate(estimate = as.double(estimate))

  return(table_processed)
}

read_covariance_table <- function(lines) {
  # Process input
  input_table <- as_tibble(fread(text=lines, sep=":", header=F, nrows = 3, fill = TRUE))

  # Perform method
  table_processed <- input_table %>%
    separate(V2, c("estimate", "stderr"), sep = " ") %>%
    mutate(stderr = as.double(str_extract(stderr, regex("\\d+\\.?\\d*")))) %>%
    rename(c("name" = "V1")) %>%
    mutate(estimate = as.double(estimate))

  return(table_processed)
}

read_ldsc_logs <- function(filepath) {
  # Open the file
  con <- file(filepath, "r")

  # List of heritability tables
  heritability_tables <- list()
  covariance_tables <- list()

  correlation_table <- NULL
  sumstats <- ""

  while ( TRUE ) {

    # Read first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }

    if (startsWith(line, "Reading summary statistics from")) {
      match <- str_match(line, "Reading summary statistics from (.+) ...")
      sumstats <- match[2]

    } else if (startsWith(line, "Heritability of phenotype")) {
      heritability_tables[[sumstats]] <- read_heritability_table(
        readLines(con, n = 6)[2:6])

    } else if (startsWith(line, "Genetic Covariance")) {
      covariance_tables[[sumstats]] <- read_covariance_table(
        readLines(con, n = 4)[2:4])

    } else if (startsWith(line, "Summary of Genetic Correlation Results")) {
      correlation_table <- as_tibble(fread(text=readLines(con, length(heritability_tables) + 1), header=T))

    }
  }

  close(con)
  return(correlation_table)
}

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
  processed_table <- tibble(read_ldsc_logs(ldsc_log)) %>%
    mutate(gene_id = gene_id)

  # Process output
  write.table(processed_table, "ldsc_matrix_h2.txt", col.names = F, row.names = F, sep = "\t", quote = F)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}