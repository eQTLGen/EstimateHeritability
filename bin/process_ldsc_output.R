#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(argparse)

# Declare constants

# Declare function definitions
read_heritability_table <- function(lines) {
  # Process input
  input_table <- as_tibble(fread(text=lines, sep=":", header=F, nrows = 5, fill = TRUE))

  mapping <- c("Total Observed scale h2" = "h2_obs", "Intercept" = "h2_int")

  # Perform method
  table_processed <- input_table %>%
    filter(V1 %in% names(mapping)) %>%
    mutate(V1 = mapping[V1]) %>%
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
    select(c("h2_obs" = "Total Observed scale h2", "h2_int" = "Intercept")) %>%
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
  sumstats <- c()
  current_sumstats <- ""
  n_variants <- c()

  error <- FALSE

  while ( TRUE ) {

    # Read first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }

    if (startsWith(line, "Reading summary statistics from")) {
      match <- str_match(line, "Reading summary statistics from (.+) ...")
      sumstats <- c(sumstats, match[2])

    } else if (startsWith(line, "Read summary statistics for")) {
      match <- str_match(line, "Read summary statistics for (\\d+) SNPs.")
      n_variants[sumstats[1]] <- match[2]

    } else if (startsWith(line, "Total Observed scale h2")) {
      current_sumstats <- sumstats[1]

      message(current_sumstats)
      table <- read_heritability_table(
        c(line, readLines(con, n = 4)))

      heritability_tables[[sumstats[1]]] <- table

    } else if (startsWith(line, "Heritability of phenotype")) {
      match <- str_match(line, "Heritability of phenotype (\\d+)(\\/\\d+)?")
      current_sumstats <- sumstats[as.numeric(match[2])]

      message(current_sumstats)
      table <- read_heritability_table(
        readLines(con, n = 6)[2:6])

      heritability_tables[[sumstats[as.numeric(match[2])]]] <- table

    } else if (startsWith(line, "Genetic Covariance")) {
      next

      covariance_tables[[sumstats[as.numeric(match[2])]]] <- read_covariance_table(
        readLines(con, n = 4)[2:4])

    } else if (startsWith(line, "Summary of Genetic Correlation Results") & length(heritability_tables) > 0) {
      correlation_table <- as_tibble(fread(text=readLines(con, length(heritability_tables) + 1), header=T))

    } else if (startsWith(line, "ERROR")) {
      error <- TRUE
      break
    }
  }

  close(con)

  if (error) {
      return(NULL)
  }

  heritability_table <- bind_rows(heritability_tables[1], .id="p1") %>%
    rename("se" = "stderr") %>%
    pivot_wider(id_cols = p1, values_from = c("estimate", "se"), names_from = "name", names_glue = "{name}_{.value}") %>%
    rename_with(~str_remove(., '_estimate')) %>%
    mutate(p2 = p1, n_variants = n_variants[p1])

  if (is.null(correlation_table)) {
    return(heritability_table)
  }

  return(bind_rows(correlation_table, heritability_table))
}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  parser <- ArgumentParser(description='process ldsc logs')
  # Basic LD Score Estimation Flags'
  # Filtering / Data Management for LD Score
  parser$add_argument('--out', type="character",
                      help='output')
  parser$add_argument('--genes', type="character",
                      nargs="+")
  parser$add_argument('--ldsc-logs', default=NULL, type="character", nargs="+")
  args <- parser$parse_args(argv)

  # Process input
  processed_tables <- mapply(function(ldsc_log, gene) {
    table_proc <- read_ldsc_logs(ldsc_log) %>%
      mutate(gene_id = gene)
    return(table_proc)
  }, args$ldsc_logs, args$genes, SIMPLIFY=F, USE.NAMES = F)

  # Process output
  write.table(bind_rows(processed_tables), args$out, col.names = T, row.names = F, sep = "\t", quote = F)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
