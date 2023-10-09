#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)

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

  while ( TRUE ) {

    # Read first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }

    if (startsWith(line, "Reading summary statistics from")) {
      match <- str_match(line, "Reading summary statistics from (.+) ...")
      sumstats <- c(sumstats, match[2])

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

    }
  }

  close(con)

  heritability_table <- bind_rows(heritability_tables[1], .id="p1") %>%
    rename("se" = "stderr") %>%
    pivot_wider(id_cols = p1, values_from = c("estimate", "se"), names_from = "name", names_glue = "{name}_{.value}") %>%
    rename_with(~str_remove(., '_estimate')) %>%
    mutate(p2 = p1)

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
    write.table(processed_table, sprintf("ldsc_matrix_%s_%s_h2.txt", gene_id, annot), col.names = T, row.names = F, sep = "\t", quote = F)
  } else {
    write.table(c(""), sprintf("ldsc_matrix_%s_%s_h2.txt", gene_id, annot), col.names = F, row.names = F, sep = "\t", quote = F)
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
