#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)

# Declare constants
parser <- ArgumentParser(description='Prepare UKB-PPP data for heritability estimation')

# Basic LD Score Estimation Flags'
# Filtering / Data Management for LD Score
parser$add_argument('--input', default=NULL, type="character", nargs="+")
parser$add_argument('--variant_reference', default=NULL, type="character")
parser$add_argument('--gene', default=NULL, type="character", nargs="4")

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Args
  args <- parser$parse_args(argv)

  # Gene info
  gene_chromosome <- args$gene[1]
  gene_start_putative <- args$gene[2]
  gene_end_putative <- args$gene[3]
  gene_name <- args$gene[4]

  gene_start <- min(gene_start_putative, gene_end_putative)
  gene_end <- max(gene_start_putative, gene_end_putative)

  cis_window_start <- gene_start - 1e6
  cis_window_end <- gene_end + 1e6

  trans_window_start <- gene_start - 5e6
  trans_window_end <- gene_end + 5e6

  fwrite(tibble(chr=gene_chromosome,
                start=cis_window_start,
                end=cis_window_end,
                name=gene_name),
         "cis.bed", sep="\t", col.names=F, row.names=F)

  fwrite(tibble(chr=gene_chromosome,
                start=trans_window_start,
                end=trans_window_end,
                name=gene_name),
         "trans.bed", sep="\t", col.names=F, row.names=F)

  # Variant Reference
  variant_reference <- readRDS(args$variant_reference)

  # Process input
  #significance_threshold <- -log(1.7e-11, 10)

  # Summary stats
  summary_stats <- bind_rows(mapply(fread, args$input, SIMPLIFY=FALSE)) %>%
    inner_join(variant_reference) %>%
    filter((REF == ALLELE0 & ALT == ALLELE1) | (REF == ALLELE1 & ALT == ALLELE0)) %>%
    group_by(RSID) %>% filter(n()==1)

  # Lead effects
  summary_stats_cis <- summary_stats %>% filter(
    gene_chromosome == CHROM,
    between(GENPOS, cis_window_start, cis_window_end)) %>%
    mutate(Z = BETA/SE, P=10^(-LOG10P)) %>%
    select(c("SNP"="RSID", "N"="N", "Z"="Z", "P"="P", "BETA"="BETA", "SE"="SE", "A1"="ALLELE0", "A2"="ALLELE1"))

  summary_stats_trans <- summary_stats %>% filter(
    !(gene_chromosome == CHROM & between(GENPOS, trans_window_start, trans_window_end))) %>%
    mutate(Z = BETA/SE, P=10^(-LOG10P)) %>%
    select(c("SNP"="RSID", "N"="N", "Z"="Z", "P"="P", "BETA"="BETA", "SE"="SE", "A1"="ALLELE0", "A2"="ALLELE1"))

  fwrite(summary_stats_trans, "sumstats_hm3.trans_all.csv.gz", sep="\t", col.names=T, row.names=F)
  fwrite(summary_stats_cis, "sumstats_hm3.cis_all.csv.gz", sep="\t", col.names=T, row.names=F)

  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

# SNP     N       Z       P       BETA    SE      A1      A2
# rs3131972       8228.0  -0.1338368835871334     0.8935315606718679      -0.0010302738555188     0.0076979815123092      G       A
# rs1048488       8888.0  -0.3216649805853721     0.7477065103766418      -0.002615999033942      0.0081326821128658      T       C
# rs3115850       8888.0  -0.2290246308471052     0.8188497722101077      -0.0017742573214724     0.0077470153097066      C       T
# rs2286139       9493.0  -0.4136168743380848     0.679154713230194       -0.0033526268147041     0.0081056335529574      T       C
# rs12562034      10961.0 -1.6195031882107394     0.1053390410708316      -0.0127568745525016     0.007877029601032       A       G