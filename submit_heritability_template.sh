#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="HeritabilityEstimation"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load jdk/16.0.1
module load openjdk/11.0.2
module load any/singularity
module load squashfs/4.4

# Define paths
nextflow_path=/gpfs/space/GI/eQTLGen/tools/

empirical="/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPCNoExpPC_2023-07-13_PerGene"
genes="/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPCNoExpPC_2023-07-13_PerGene/phenotypes_unique.txt"

variant_reference="/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/MetaAnalysis/bin/hase/data/1000G-30x.ref.gz"
gene_reference="/gpfs/space/GI/eQTLGen/freeze1/InputFilesForPaper/2023-01-28_MetaAnalysis/data/Homo_sapiens.GRCh38.106.gtf.gz"

ld_w_dir="/gpfs/space/GI/eQTLGen/public_data/ldsc/baselineLF_v2.2.UKB"
variants=""

output_folder="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/ldsc_heritability/output_empirical_4GenPCNoExpPC_2023-07-13"

NXF_VER=21.10.6 ${nextflow_path}/nextflow run main.nf \
--input ${empirical}/eqtls \
--genes ${genes} \
--variant_reference ${variant_reference} \
--gene_reference ${gene_reference} \
--output ${output_folder} \
--ld_w_dir ${ld_w_dir} \
-resume \
-profile slurm,singularity


