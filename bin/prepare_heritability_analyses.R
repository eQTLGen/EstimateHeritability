#!/usr/bin/env Rscript


# Load libraries
library(arrow)
library(tidyverse)
library(data.table)
library(argparse)


# Declare constants
parser <- ArgumentParser(description='Prepare eQTLGen sumstats for heritability analyses')
trans_window <- 5e6
cis_window <- 1e6
polygenic_window <- 5e6

# Basic LD Score Estimation Flags'
# Filtering / Data Management for LD Score
parser$add_argument('--input', default=NULL, type="character", nargs="+")
parser$add_argument('--lead-variants', default=NULL, type="character")
parser$add_argument('--variant-reference', default=NULL, type="character")
parser$add_argument('--variant-list', default=NULL, type="character")
parser$add_argument('--genes', default=NULL, type="character", nargs="+")
parser$add_argument('--n-min', default=NULL, type="integer")
parser$add_argument('--gene-reference', default=NULL, type="character")

get_variants_in_qtl_windows <- function(lead_trans_effects, hm3_variant_ref) {

  setDT(lead_trans_effects)
  setDT(hm3_variant_ref)

  # prepare
  hm3_variant_ref[, chromosome := variant_chr]
  hm3_variant_ref[, start := variant_pos]
  hm3_variant_ref[, end := variant_pos]

  setkey(hm3_variant_ref, chromosome, start, end)
  setkey(lead_trans_effects, chromosome, lead_start, lead_end)

  # exclude variants
  qtl_variants <- foverlaps(
    hm3_variant_ref,
    lead_trans_effects,
    by.x=c("chromosome","start","end"),
    by.y=c("chromosome","lead_start","lead_end"),
    nomatch=0
  )

  print(head(qtl_variants))

  qtl_variants_summarised <- qtl_variants[
    , .(variants = list(variant_id)),
      by = phenotype
  ]

  return(qtl_variants_summarised)
}

# Declare function definitions
write_bed <- function(df, filename, chrom_col = "chromosome", start_col = "start", end_col = "end", gene_col = "gene_id") {
  bed_df <- df %>%
    select(!!chrom_col, !!start_col, !!end_col, !!gene_col) %>%
    arrange(!!chrom_col, !!start_col)
  fwrite(bed_df, file = filename, sep = "\t", col.names = FALSE)
}

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

  # Gene reference
  gene_ref <- fread(args$gene_reference)

  # eQTL dataset
  eqtl_ds <- arrow::open_dataset(args$input)
  genes <- args$genes

  # Variant list
  variant_dt <- fread(args$variant_list, col.names = c("variant"), header=F)

  # Variant reference
  variant_reference <- arrow::read_parquet(args$variant_reference)

  hm3_variant_ref <- variant_reference %>%
    filter(variant %in% variant_dt$variant) %>%
    rename(variant_chr = chromosome, variant_pos = bp, variant_id = variant)

  # Save gene reference as df
  gene_ref_df <- as.data.frame(gene_ref) %>%
    filter(gene_id %in% genes) %>%
    select(gene_id, start, end, seqid, gene_name) %>%
    filter(seqid %in% c(1:22, "X", "Y", "XY", "MT")) %>%
    mutate(chromosome = as.integer(case_when(
      seqid == "X" ~ "23",
      seqid == "Y" ~ "24",
      seqid == "XY" ~ "25",
      seqid == "MT" ~ "26",
      TRUE ~ as.character(seqid))))

  # For every gene, get the cis-window, trans-window
  gene_windows <- gene_ref_df %>%
    mutate(
      cis_start = pmax(start - cis_window, 1),
      cis_end = end + cis_window,
      trans_start = pmax(start - trans_window, 1),
      trans_end = end + trans_window
    )

  # Write cis-window and trans-window BED files for each gene
  cis_bed_file <- "cis.bed"
  trans_bed_file <- "trans.bed"
  write_bed(gene_windows, cis_bed_file, start_col = "cis_start", end_col = "cis_end")
  write_bed(gene_windows, trans_bed_file, start_col = "trans_start", end_col = "trans_end")

  # For each gene, get lead trans effects and write to BED
  lead_variants <- fread(args$lead_variants)
  lead_trans_effects <- lead_variants %>%
    semi_join(gene_windows, by = c("phenotype" = "gene_id")) %>%
    mutate(
      lead_start = pmax(bp - polygenic_window, 1),
      lead_end = bp + polygenic_window)

  lead_bed_file <- "polygenic.bed"
  write_bed(lead_trans_effects, lead_bed_file, start_col = "lead_start", end_col = "lead_end", gene_col = "phenotype")

  qtl_variants <- get_variants_in_qtl_windows(lead_trans_effects, hm3_variant_ref)

  # For all windows, get the correct variant set from the variant index
  gene_ref_df <- gene_windows %>% rowwise() %>%
    mutate(
      cis_variants = list(
        hm3_variant_ref %>%
          filter(
            variant_chr == chromosome,
            variant_pos >= cis_start,
            variant_pos <= cis_end
          ) %>%
          pull(variant_index)
      ),
      trans_variants = list(
        hm3_variant_ref %>%
          filter(
            variant_chr != chromosome | (variant_pos < trans_start | variant_pos > trans_end)
          ) %>%
          pull(variant_index)
      )
    ) %>%
    ungroup()

  ldsc_selector <- c("SNP"="variant_id", "N"="sample_size", "Z"="z_score", "A1"="eff_allele", "A2"="non_eff_allele")

  # For every gene, extract the variants of interest
  for (gene in gene_ref_df$gene_id) {
    gene_data <- gene_ref_df %>% filter(gene_id == gene)
    qtl_variants_focal_gene <- unlist(qtl_variants %>% filter(phenotype == gene) %>% pull(variants))

    cis_variants <- unlist(gene_data$cis_variants)
    trans_variants <- unlist(gene_data$trans_variants)

    summary_stats <- eqtl_ds %>% filter(
      phenotype == gene_data$gene_id,
      variant_index %in% hm3_variant_ref$variant_index) %>% collect() %>%
      inner_join(hm3_variant_ref, by = "variant_index") %>%
      mutate(z_score = beta / standard_error)

    print(summary_stats)

    mean_sample_size <- mean(summary_stats$sample_size, na.rm=T)
    sd_sample_size <- sd(summary_stats$sample_size, na.rm=T)
    max_sample_size <- max(summary_stats$sample_size, na.rm=T)
    print(mean_sample_size)
    print(max_sample_size)

    summary_stats <- summary_stats %>% filter(between(sample_size, max_sample_size * 0.95, max_sample_size + 1))

    trans_summary_stats <- summary_stats %>% filter(variant_index %in% trans_variants) %>%
      select(all_of(ldsc_selector))

    polygenic_summary_stats <- summary_stats %>%
      filter(
        variant_index %in% trans_variants,
        !variant_id %in% qtl_variants_focal_gene
      ) %>% select(all_of(ldsc_selector))

    print(gene)
    print(length(qtl_variants_focal_gene))
    print(nrow(trans_summary_stats))
    print(nrow(polygenic_summary_stats))
    print(nrow(summary_stats))

    fwrite(trans_summary_stats, sprintf("%s.sumstats_hm3.trans_all.csv.gz", gene), sep=" ", row.names = FALSE)
    fwrite(polygenic_summary_stats, sprintf("%s.sumstats_hm3.trans_polygenic.csv.gz", gene), sep=" ", row.names = FALSE)
    fwrite(summary_stats %>% select(all_of(ldsc_selector)), sprintf("%s.sumstats_hm3.gw_polygenic.csv.gz", gene), sep=" ", row.names = FALSE)
  }

  # Annotate those variants that

  # For all genes, get the correct variant set from the eqtl dataset


  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}