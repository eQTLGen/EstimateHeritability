#!/bin/bash

#SBATCH --time=10:00:00
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
module load squashfs
module load singularity

# Define paths
nextflow_path=/gpfs/space/GI/eQTLGen/tools/

ldsc_sources="/gpfs/space/GI/eQTLGen/freeze2/tools/ldsc"

mastertable="/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/input/MasterTable_newCohorts_2023-07-13_correctedCohortNames.txt"
empirical="/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPCNoExpPC_2023-07-13_PerGene"
#empirical="/gpfs/space/home/warmerda/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPC20ExpPC_2023-05-27_PerGene"
genes="/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPCNoExpPC_2023-07-13_PerGene/phenotypes_unique.txt"
inclusion_step_output="/gpfs/space/GI/eQTLGen/freeze2/InclusionLists/output_2023-10-03_all"

variant_reference="/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/MetaAnalysis/bin/hase/data/1000G-30x.ref.gz"
gene_reference="/gpfs/space/GI/eQTLGen/freeze1/InputFilesForPaper/2023-01-28_MetaAnalysis/data/Homo_sapiens.GRCh38.106.gtf.gz"

ld_w_dir="/gpfs/space/GI/eQTLGen/public_data/ldsc/sldsc_files_hg38/baselineLD_v2.2"
frqfile_dir="/gpfs/space/GI/eQTLGen/public_data/ldsc/sldsc_files_hg38/plink_files"
weights_dir="/gpfs/space/GI/eQTLGen/public_data/ldsc/sldsc_files_hg38/weights"

variants="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/heritability/EstimateHeritability/data/matching_unambiguous_hapmap_variants.txt"
hapmap="/gpfs/space/GI/eQTLGen/public_data/ldsc/w_hm3.snplist"

gwas_map="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/heritability/EstimateHeritability/gwas_studies_map.csv"
maf_table="/gpfs/space/GI/eQTLGen/freeze2/PreMetaQc/output/MafInformation_2023-10-03.txt.gz"

output_folder="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/heritability/output_empirical_4GenPCNoExpPC_2023-07-13/sldsc_2024-01-30"

variants_bed="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/heritability/EstimateHeritability/data/variants.chrALL.hg38.w_hm3.M_5_50.bed"

NXF_VER=21.10.6 ${nextflow_path}/nextflow run main.nf \
--input ${empirical}/eqtls \
--genes ${genes} \
--mastertable ${mastertable} \
--variant_reference ${variant_reference} \
--gene_reference ${gene_reference} \
--output ${output_folder} \
--ld_w_dir ${ld_w_dir} \
--frqfile_dir ${frqfile_dir} \
--weights_dir ${weights_dir} \
--variants ${variants} \
--hapmap ${hapmap} \
--variants_bed ${variants_bed} \
--gwas_map ${gwas_map} \
--ldsc_source ${ldsc_sources} \
--inclusion_step_output ${inclusion_step_output} \
--maf_table ${maf_table} \
-resume \
-profile slurm,singularity


