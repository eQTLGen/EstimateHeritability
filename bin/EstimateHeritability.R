#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(arrow)
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)
library(SumTool)

setDTthreads(1)

# Functions
RemoveSigEffects <- function(data, window = 1000000, Pthresh = 5e-8, snp_id_col = "SNP", snp_chr_col = "Chr", snp_pos_col = "Pos", beta_col = "beta", se_col = "se") {

  data <- data.table(SNP = data[[snp_id_col]],
                       snp_chr = data[[snp_chr_col]],
                       snp_pos = data[[snp_pos_col]],
                       beta = data[[beta_col]],
                       se = data[[se_col]])

  data$P <- ZtoP(data$beta/data$se)
  data$Z <- data$beta/data$se

  #data_f <- data[data$P < as.numeric(Pthresh), ]
  data_f <- data
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]

  while (min(data_f$P) <= 5e-8) {
    lead_snp <- data_f[abs(data_f$Z) == max(abs(data_f$Z)), ]
    if (nrow(lead_snp) > 1) {
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$snp_chr == lead_snp$snp_chr & data_f$snp_pos > lead_snp$snp_pos - window & data_f$snp_pos < lead_snp$snp_pos + window), ]
    # message(paste("Added:", lead_snp$snp_chr, lead_snp$snp_pos))
  }
  return(data_f)
}

ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    message("Converting the non-numeric vector to numeric vector.")
    Z <- as.numeric(Z)
  }

  if (largeZ == TRUE) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)

    if (largeZ == TRUE & log10P == TRUE) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)

    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
    }
  }

  return(P)
}


# Paths
gtf <- args[1]
ld_path <- args[2]
snp_ref <- args[3]

assoc_parquet <- args[4]
log_path <- args[5]

# Arguments
sig_pthresh <- as.numeric(args[6])
locus_window <- as.numeric(args[7])

## Read in reference files
# gtf

gtf <- fread(gtf, skip = 5)
gtf <- gtf[gtf$V3 %in% "gene", ]
gtf$gene_id <- gsub(";.*", "", gtf$V9)
gtf$gene_id <- gsub("gene_id ", "", gtf$gene_id)
gtf$gene_id <- gsub("\"", "", gtf$gene_id)

gtf$gene_name <- gsub(".*; gene_name ", "", gtf$V9)
gtf$gene_name <- gsub(";.*", "", gtf$gene_name)
gtf$gene_name <- gsub("gene_id .*", "", gtf$gene_name)
gtf$gene_name <- gsub("\"", "", gtf$gene_name)

colnames(gtf)[c(1, 4, 5, 10, 11)] <- c("seqid", "start", "end", "gene_id", "gene_name")
gtf <- unique(gtf[, colnames(gtf) %in% c("gene_id", "gene_name", "seqid", "start", "end"), with = FALSE])
message("GTF read!")

ld <- fread2(paste0(ld_path, "/", 1:22, ".l2.ldscore.gz"), data.table = TRUE)

message("LD scores read!")

hapmap3 <- fread(paste0(ld_path, "/", "w_hm3.snplist"))

ds <- arrow::open_dataset(snp_ref)
snp_ref <- ds %>% dplyr::filter(ID %in% !!hapmap3$SNP) %>% collect()

setkey(snp_ref, cols = ID)
message("SNP reference data read!")

## Read in log file
log <- fread(log_path)
colnames(log)[1:ncol(log) - 1] <- colnames(log)[2:ncol(log)]
log <- log[, -ncol(log), with = FALSE]

log <- log[order(log$N, decreasing = TRUE), ]

## Read in summary statistics
# get gene name:
gene_i <- gsub(".*phenotype=", "", assoc_parquet)
message(paste("Analysing", gene_i))
gtf_f <- gtf[gtf$gene_id == gene_i, ]

# If gene is not in gtf, stop and write out empty files
if (nrow(gtf_f) == 0){
  
  res <- data.table(gene = NA,
  hgnc = NA,
  int = NA,
  int_se = NA,
  P_int = NA,
  h2 = NA,
  h2_se = NA,
  P_h2 = NA,
  type = NA)
  res <- res[-1, ]
  fwrite(res, paste0(gene_i, "_h2.txt"), sep = "\t", quote = FALSE)

  png(paste0(gene_i, "_h2.png"), height = 5, width = 5, units = "in", res = 100, type = "cairo")
  plot(1,1, col = "white")
  dev.off()

  quit()
  }

and <- read_parquet(list.files(assoc_parquet, full.names = TRUE))

###########
# Analysis#
###########
and <- as.data.table(and)
and <- and[i_squared < 40, ]

chi2 <- (and$beta/and$standard_error)^2
lambda_gc <- median(chi2) / qchisq(0.5,1)
mean_chi2 <- mean(chi2)

and <- and[variant %in% hapmap3$SNP, ]

# If gene is tested in >1 dataset, remove variants tested in only one
if (length(unique(and$N)) > 1){
and <- and[!(and$sample_size %in% log$N), ]
}
and <- setkey(and, cols = "variant")

# If assoc file is empty, stop and write out empty files
if (nrow(and) == 0){
  
  res <- data.table(gene = NA,
  hgnc = NA,
  int = NA,
  int_se = NA,
  P_int = NA,
  h2 = NA,
  h2_se = NA,
  P_h2 = NA,
  lambda = NA,
  mean_chi2 = NA,
  N = NA,
  N_variants = NA,
  type = NA)

  res <- res[-1, ]
  fwrite(res, paste0(gene_i, "_h2.txt"), sep = "\t", quote = FALSE)

  png(paste0(gene_i, "_h2.png"), height = 5, width = 5, units = "in", res = 100, type = "cairo")
  plot(1,1, col = "white")
  dev.off()

  quit()
  }

and <- merge(and, snp_ref, by.x = "variant", by.y = "ID")

message("Ref. added.")

# Check if alleles are same compared to hapmap3 reference
and <- merge(and, hapmap3, by.x = "variant", by.y = "SNP")
and[and$str_allele1 == and$A1 & and$str_allele2 == and$A2, ]$beta <- -and[and$str_allele1 == and$A1 & and$str_allele2 == and$A2, ]$beta

# remove those SNPs which do not suit with HapMap3 alleles
and <- and[((and$str_allele1 == and$A1 | and$str_allele1 == and$A2) & (and$str_allele2 == and$A1 | and$str_allele2 == and$A2)), ]

# Flip the alleles for inverted SNPs
and[(and$str_allele1 == and$A2), ]$str_allele1 <- and[(and$str_allele1 == and$A2), ]$A1
and[(and$str_allele2 == and$A1), ]$str_allele2 <- and[(and$str_allele2 == and$A1), ]$A2

inp1 <- data.table(SNP = and$variant,
                    Chr = and$CHR,
                    Pos = and$bp,
                    A1 = and$str_allele1,
                    A2 = and$str_allele2,
                    BETA = and$beta,
                    SE = and$standard_error,
                    N = and$sample_size)

inp1 <- inp1[inp1$N > quantile(inp1$N, 0.9) / 1.5, ]

sample_size <- max(inp1$N)

ld2 <- merge(hapmap3, ld)

ld2 <- merge(ld2[, -c(6), with = FALSE], inp1[, c(1, 3), with = FALSE], by.x = "SNP", by.y = "SNP")

ld2 <- data.table(SNP = ld2$SNP,
                    Chr = ld2$CHR,
                    Pos = ld2$Pos,
                    A1 = ld2$A1,
                    A2 = ld2$A2,
                    Maf = ld2$MAF,
                    ldscore = ld2$L2)

# Remove HLA region
inp1 <- inp1[!(inp1$Chr == 6 & inp1$Pos > 25726063 & inp1$Pos < 33400644), ]

# How many variants
inp1_SNP2 <- paste(inp1$Chr, inp1$Pos, inp1$A1, inp1$A2)
ld2_SNP2 <- paste(ld2[Maf > 0.01]$Chr, ld2[Maf > 0.01]$Pos, ld2[Maf > 0.01]$A1, ld2[Maf > 0.01]$A2)

SNP_N_all <- length(intersect(inp1_SNP2, ld2_SNP2))

message("Data prepared!")
message("LDSC for full data")

res1 <- LDreg(sumstat = inp1, 
ldscore = ld2,
nblock = 200,
maf = 0.01,
maxz2 = 30)

names(res1)[c(2, 4)] <- c("int_se", "h2_se")

chi2 <- (inp1$BETA/inp1$SE)^2
lambda_gc <- median(chi2) / qchisq(0.5,1)
mean_chi2 <- mean(chi2)
ratio <- (mean_chi2 - 1) / (res1["Intercept"] - 1)

res1 <- data.table(gene = gene_i,
hgnc = gtf_f$gene_name, 
int = res1["Intercept"],
int_se = res1["int_se"],
P_int = pnorm(q = res1["Intercept"] - 1 / res1["int_se"],  lower.tail = FALSE),
h2 = res1["Heritability"],
h2_se = res1["h2_se"],
P_h2 = pnorm(q = res1["Heritability"] / res1["h2_se"],  lower.tail = FALSE),
ratio = ratio,
lambda = lambda_gc,
mean_chi2 = mean_chi2,
N = sample_size,
N_variants = SNP_N_all,
type = "all")

# comb <- merge(inp1[, c(1, 6, 7), with = FALSE], ld2[Maf > 0.01, c(1, 7), with = FALSE], by = "SNP")

# group <- cut(comb$ldscore, 50, labels = F)
# ld_block1 <- tapply(comb$ldscore, group, mean)
# chi2_block1 <- tapply((comb$BETA/comb$SE) * (comb$BETA/comb$SE), group, mean)

message("LDSC for data without cis-eQTL locus")
inp2 <- inp1[!(inp1$Chr == gtf_f$seqid & inp1$Pos > gtf_f$start - locus_window & inp1$Pos < gtf_f$end + locus_window), ]

# How many variants
inp2_SNP2 <- paste(inp2$Chr, inp2$Pos, inp2$A1, inp2$A2)
ld2_SNP2 <- paste(ld2[Maf > 0.01]$Chr, ld2[Maf > 0.01]$Pos, ld2[Maf > 0.01]$A1, ld2[Maf > 0.01]$A2)

SNP_N_cis <- length(intersect(inp2_SNP2, ld2_SNP2))

res2 <- LDreg(sumstat = inp2, 
ldscore = ld2,
nblock = 200,
maf = 0.01,
maxz2 = 30)

names(res2)[c(2, 4)] <- c("int_se", "h2_se")

chi2 <- (inp2$BETA/inp2$SE)^2
lambda_gc <- median(chi2) / qchisq(0.5,1)
mean_chi2 <- mean(chi2)
ratio <- (mean_chi2 - 1) / (res2["Intercept"] - 1)

res2 <- data.table(gene = gene_i,
hgnc = gtf_f$gene_name, 
int = res2["Intercept"],
int_se = res2["int_se"],
P_int = pnorm(q = res2["Intercept"] - 1 / res2["int_se"],  lower.tail = FALSE),
h2 = res2["Heritability"],
h2_se = res2["h2_se"],
P_h2 = pnorm(q = res2["Heritability"] / res2["h2_se"],  lower.tail = FALSE),
ratio = ratio,
lambda = lambda_gc,
mean_chi2 = mean_chi2,
N = sample_size,
N_variants = SNP_N_cis,
type = "cis")

# group <- cut(comb2$ldscore, 50, labels = F)
# ld_block2 <- tapply(comb2$ldscore, group, mean)
# chi2_block2 <- tapply((comb2$BETA/comb2$SE) * (comb2$BETA/comb2$SE), group, mean)

message("LDSC for data without any eQTL locus")

polygenic_background <- RemoveSigEffects(inp1, snp_id_col = "SNP", 
snp_chr_col = "Chr", snp_pos_col = "Pos", 
beta_col = "BETA", se_col = "SE", Pthresh = sig_pthresh, 
window = locus_window)

inp3 <- inp1[inp1$SNP %in% polygenic_background$SNP, ]
message("Data distance pruned")

# How many variants
inp3_SNP2 <- paste(inp3$Chr, inp3$Pos, inp3$A1, inp3$A2)
ld2_SNP2 <- paste(ld2[Maf > 0.01]$Chr, ld2[Maf > 0.01]$Pos, ld2[Maf > 0.01]$A1, ld2[Maf > 0.01]$A2)

SNP_N_omni <- length(intersect(inp3_SNP2, ld2_SNP2))

res3 <- LDreg(sumstat = inp3, 
ldscore = ld2,
nblock = 200,
maf = 0.01,
maxz2 = 30)

names(res3)[c(2, 4)] <- c("int_se", "h2_se")

chi2 <- (inp3$BETA/inp3$SE)^2
lambda_gc <- median(chi2) / qchisq(0.5,1)
mean_chi2 <- mean(chi2)
ratio <- (mean_chi2 - 1) / (res3["Intercept"] - 1)

res3 <- data.table(gene = gene_i,
hgnc = gtf_f$gene_name, 
int = res3["Intercept"],
int_se = res3["int_se"],
P_int = pnorm(q = res3["Intercept"] - 1 / res3["int_se"],  lower.tail = FALSE),
h2 = res3["Heritability"],
h2_se = res3["h2_se"],
P_h2 = pnorm(q = res3["Heritability"] / res3["h2_se"],  lower.tail = FALSE),
ratio = ratio,
lambda = lambda_gc,
mean_chi2 = mean_chi2,
N = sample_size,
N_variants = SNP_N_omni,
type = "omnigenic")

res <- rbind(res1, res2, res3)

fwrite(res, paste0(gene_i, "_h2.txt"), sep = "\t", quote = FALSE)

# png(paste0(gene_i, "_h2.png"), height = 5, width = 15, units = "in", res = 100, type = "cairo")
# par(mfrow = c(1, 3))

# plot(ld_block1, chi2_block1, xlab = "Mean LD score", ylab = "Mean chi2", pch = 19, main = paste0(gene_i, "\n", gtf_f$gene_name, "\nall effects"))
# abline(lm(chi2_block1 ~ ld_block1), col = "blue", lty = 2)
# plot(ld_block2, chi2_block2, xlab = "Mean LD score", ylab = "Mean chi2", pch = 19, main = paste0(gene_i, "\n", gtf_f$gene_name, "\nwithout cis-eQTL locus"))
# abline(lm(chi2_block2 ~ ld_block2), col = "blue", lty = 2)
# plot(ld_block3, chi2_block3, xlab = "Mean LD score", ylab = "Mean chi2", pch = 19, main = paste0(gene_i, "\n", gtf_f$gene_name, "\nwithout any eQTL locus"))
# abline(lm(chi2_block3 ~ ld_block3), col = "blue", lty = 2)

# dev.off()
# message("Analysis done!")
