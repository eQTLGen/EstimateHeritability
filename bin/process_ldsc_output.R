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
      break

    } else if (startsWith(line, "Genetic Covariance")) {
      covariance_tables[[sumstats]] <- read_covariance_table(
        readLines(con, n = 4)[2:4])

    } else if (startsWith(line, "Summary of Genetic Correlation Results") & length(heritability_tables) > 0) {
      correlation_table <- as_tibble(fread(text=readLines(con, length(heritability_tables) + 1), header=T))

    }
  }

  close(con)
  return(bind_rows(heritability_tables, .id="sumstats_id"))
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
  annot <- argv[2]
  ldsc_log <- argv[3]

  # Process input
  table_proc <- read_ldsc_logs(ldsc_log)

  if (!is.null(table_proc) & nrow(table_proc)>0) {

    processed_table <- table_proc %>%
      mutate(gene_id = gene_id,
             annot = annot)

    # Process output
    write.table(processed_table, sprintf("ldsc_matrix_%s_%s_h2.txt", gene_id, annot), col.names = F, row.names = F, sep = "\t", quote = F)
  } else {
    write.table(c(""), sprintf("ldsc_matrix_%s_%s_h2.txt", gene_id, annot), col.names = F, row.names = F, sep = "\t", quote = F)
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
