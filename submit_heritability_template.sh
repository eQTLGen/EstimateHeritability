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
nextflow_path=[path to nextflow executble]
set -f

input_path=[Input folder with eQTL results in .parquet format]
output_folder=[Your output folder]

NXF_VER=21.10.6 ${nextflow_path}/nextflow run main.nf \
--inputdir ${input_path} \
--outdir ${output_folder} \
-resume \
-profile slurm,singularity
