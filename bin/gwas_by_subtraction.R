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
parser$add_argument('--rg', default=NULL, type="character", nargs="+",
                    help='List of sumstat files for correction.')
parser$add_argument('--ref-ld-chr', default=NULL, type="character")
parser$add_argument('--w-ld-chr', default=NULL, type="character")
parser$add_argument('--ref', default=NULL, type="character",
                    help='1kg reference')
parser$add_argument('--names', default=NULL, type="character", nargs="+",
                    help='List of names for traits associated to the sumstat files.')

# Declare function definitions
subtract_ext <- function() {

  latent_factors <- sprintf("%s=~NA*GE_CTC + start(0.4)*%s", latent_factor_names, trait_names)
  ge_latent_factor <- "GE=~NA*GE_CTC +start(0.2)*GE_CTC"

  snp_dependency <- sprintf("%s~SNP", latent_factor_names)
  ge_snp_dependency <- "GE~SNP"

  ge_dependency <- "GE~~1*GE"

  latent_factors <- sprintf("RET~~1*RET\nRET~~0*GE", latent_factor_names)

  model_ext <-'RET=~NA*GE_CTC + start(0.4)*RET_G
          ERY=~NA*GE_CTC + start(0.4)*ERY_G
          PLT=~NA*GE_CTC + start(0.4)*PLT_G
          LYM=~NA*GE_CTC + start(0.4)*LYM_G

          GE=~NA*GE_CTC +start(0.2)*GE_CTC

          RET~SNP
          ERY~SNP
          PLT~SNP
          LYM~SNP

          GE~SNP
          GE~~1*GE

          RET~~1*RET
          RET~~0*GE
          ERY~~1*ERY
          ERY~~0*GE
          PLT~~1*PLT
          PLT~~0*GE
          LYM~~1*LYM
          LYM~~0*GE

          RET_G ~~ 0*GE_CTC
          RET_G~~0*RET_G
          ERY_G ~~ 0*GE_CTC
          ERY_G~~0*ERY_G
          PLT_G ~~ 0*GE_CTC
          PLT_G~~0*PLT_G
          LYM_G ~~ 0*GE_CTC
          LYM_G~~0*LYM_G

          GE_CTC~~0*GE_CTC
          SNP~~SNP'

}

subtract_plt <- function(ldsc_output, p_sumstats) {

  # model with SNP
  model_plt <- 'PLT=~NA*GE_CTC + start(0.4)*PLT_G
          GE=~NA*GE_CTC +start(0.2)*GE_CTC

          PLT~SNP
          GE~SNP

          GE~~1*GE
          PLT~~1*PLT
          PLT~~0*GE

          PLT_G ~~ 0*GE_CTC
          PLT_G~~0*PLT_G
          GE_CTC~~0*GE_CTC
          SNP~~SNP'

  outputGWAS <- userGWAS(
    covstruc=ldsc_output,
    SNPs=p_sumstats,
    estimation="DWLS",
    model=model_plt,
    sub=c("PLT~SNP","GE~SNP"))
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
  print(argv)
  args <- parser$parse_args(argv)

  traits <- args$rg
  sample.prev <- c(NA,NA)
  population.prev <- c(NA,NA)
  ld <- args$ref_ld_chr
  wld <- args$ref_w_chr
  trait.names <- args$names
  n <- args$n

  ldsc_output <- ldsc(traits,
                     sample.prev,
                     population.prev,
                     ld,
                     wld,
                     trait.names)

  ref <- args$ref
  se.logit <- c(F,F)
  betas <- c(NULL, NULL)
  info.filter <- 0.6
  maf.filter <- 0.01

  p_sumstats <- sumstats(traits,
                         ref,
                         trait.names,
                         se.logit,
                         info.filter,
                         maf.filter,
                         OLS=c(T,T),
                         betas=betas,
                         linprob=NULL,
                         N = c(NULL, NULL))

  save(p_sumstats, file="Sumstats.RData")

  subtract_plt(ldsc_output, p_sumstats)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
