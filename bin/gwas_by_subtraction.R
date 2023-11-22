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
                    help='List of munged sumstat files for LDSC calculations.')
parser$add_argument('--sumstats', default=NULL, type="character", nargs="+",
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


subtract_ldsc_extended <- function(ldsc_output, trait_names) {

  latent_factor_names <- paste0("ctc_", trait_names)

  latent_factors <- sprintf("%s=~NA*GE + start(0.4)*%s", latent_factor_names, trait_names)
  ge_latent_factor <- "GeNonCtc=~NA*GE + start(0.2)*GE"

  ge_non_ctc_dependency <- "GeNonCtc~~1*GeNonCtc"

  latent_factor_dependencies <- sprintf("%s~~1*%1$s\\n%1$s~~0*GeNonCtc", latent_factor_names)
  trait_dependencies <- sprintf("%s~~0*%1$s\n%1$s~~0*GE", trait_names)

  ge_dependency <- "GE~~0*GE"

  model <- paste(
    latent_factors, ge_latent_factor, ge_non_ctc_dependency,
    latent_factor_dependencies, trait_dependencies, ge_dependency, sep="\n")

  output<-usermodel(ldsc_output,estimation="DWLS",model=model)
  return(output)
}


subtract_ldsc <- function(ldsc_output) {

  model<-'CTC=~NA*GE + start(0.4)*CC
          GeNonCtc=~NA*GE

          GeNonCtc~~1*GeNonCtc
          CTC~~1*CTC
          CTC~~0*GeNonCtc

          CC ~~ 0*GE
          CC~~0*CC
          GE~~0*GE'

  output<-usermodel(ldsc_output,estimation="DWLS",model=model)
  return(output)
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
  gene_reference_path <- args$gene_reference
  sumstats <- args$sumstats
  sample.prev <- c(NA,NA)
  population.prev <- c(NA,NA)
  ld <- args$ref_ld_chr
  wld <- args$ref_w_chr
  trait.names <- args$names

  ldsc_output <- ldsc(traits,
                     sample.prev,
                     population.prev,
                     ld,
                     wld,
                     trait.names)

  res <- subtract_ldsc_extended(ldsc_output, trait.names)

  print(res)
  #
  # ref <- args$ref
  # se.logit <- c(F,F)
  # betas <- c(NULL, NULL)
  # info.filter <- 0.6
  # maf.filter <- 0.01
  #
  # p_sumstats <- sumstats(sumstats,
  #                        ref,
  #                        trait.names,
  #                        se.logit,
  #                        info.filter,
  #                        maf.filter,
  #                        OLS=c(T,T),
  #                        betas=betas,
  #                        linprob=NULL,
  #                        N = c(NULL, NULL))
  #
  # save(p_sumstats, file="Sumstats.RData")
  #
  # subtract_plt(ldsc_output, p_sumstats)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
