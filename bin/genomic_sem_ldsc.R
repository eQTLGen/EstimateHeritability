#!/usr/bin/env Rscript


# Load libraries
require(GenomicSEM)

library(tidyverse)
library(data.table)
library(argparse)

# Declare constants
parser <- ArgumentParser(description='Calculate LDSC heritability using GenomicSEM')
# Basic LD Score Estimation Flags'
# Filtering / Data Management for LD Score
parser$add_argument('--rg', default=None, type=str,
                    help='List of sumstat files for genetic correlation estimation.')
parser$add_argument('--ref-ld-chr', default=None, type=str,
                    help=paste0(
                      'Same as --ref-ld, but will automatically concatenate .l2.ldscore files split ',
                      'across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz ',
                      'to the filename prefix. If the filename prefix contains the symbol @, LDSC will ',
                      'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome ',
                      'numbers to the end of the filename prefix.',
                      'Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz',
                      'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz'))
parser$add_argument('--w-ld-chr', default=None, type=str,
                    help=paste0(
                      'Same as --w-ld, but will read files split into 22 chromosomes in the same ',
                      'manner as --ref-ld-chr.'))
parser$add_argument('--names', default=None, type=str,
                    help=paste0(
                      'List of names for traits associated to the sumstat files.'))

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Process input
  args <- parser$parse_args(argv)

  traits <- args$rg
  sample.prev <- c(NA,NA)
  population.prev <- c(NA,NA)
  ld <- args$ref_ld_chr
  wld <- args$ref_w_chr
  trait.names <- args$names


  # Perform method
  LDSCoutput <- ldsc(
    traits,
    sample.prev,
    population.prev,
    ld,
    wld,
    chisq.max=10000,
    trait.names)

  print(LDSCoutput)

  save(LDSCoutput, file="LDSCoutput.RData")

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}