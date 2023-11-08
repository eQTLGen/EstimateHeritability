#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { EstimateTransHeritabilityLdsc; EstimateTransHeritabilityLdsc as EstimateGwHeritabilityLdsc; EstimateCisHeritabilityLdsc; ProcessLdscOutput as ProcessTransLdscOutput; ProcessLdscOutput as ProcessCisLdscOutput; ProcessLdscOutput as ProcessGwLdscOutput; CountHeritabilitySnps; EstimateHeritabilityLdscAllPairwise; ProcessLdscDeleteVals } from './modules/EstimateHeritability'
include { WriteOutRes } from './modules/WriteOutRes'
include { ExtractResults; ExtractResultsPerCohort; ProcessResults } from './modules/CollectResults.nf'
include { ProcessVuckovicGwasData } from './modules/ProcessGwas.nf'

def helpmessage() {

log.info"""

Estimate heritability ~ v${workflow.manifest.version}"
=================================================
Pipeline for running encoded meta-analysis over multiple studies in the central site.

Pay attention that it does not explicitly overwrite the output folder, so clean it before re-running the pipeline.

Usage:

nextflow run main.nf \
--inputdir [path to parquet folder] \
--lddir [] \
--snpref [] \
--gtf [] \
--metalog [path to dataset QC log file]
--logfile [path to logfile] \
--pthresh [P-value threshold for eQTLs] \
--eqltwindow [window for removing eQTL loci]
--outdir [path to results folder]


Mandatory arguments:
--inputdir    Path to per-gene summary statistics, as written out by HASE meta-analyser.
--lddir       Full path to the directory with HapMap3 LD scores.
--snpref      Full path to SNP reference file in parquet format.
--gtf         Full path to ENSEML .gtf file for hg38.
--pthresh     P-value threshold for removing eQTLs. Defaults to 5e-8.
--eqtlwindow  Genomic window for removing eQTLs. Defaults to (+/-) 5000000(bp).
--logfile     Full path to logfile containing sample sizes, etc.
--outdir      Path to results folder where heritability estimation outputs are written.

""".stripIndent()

}

params.variants = 'NO_FILE'
params.cols = '+z_score'
params.output

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }
Channel.fromPath(params.gwas_map)
    .splitCsv(header:true)
    .map { row-> tuple(row.Name, file(row.Path), row.N) }
    .view()
    .set { gwas_input_ch }

cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.cohort_new_name ]}

//cohorts_ch = Channel
//   .fromList( ["GTEx_2017-06-05_v8_EUR"] )

inclusion_step_output_ch = file(params.inclusion_step_output)
one_kg_bed_ch = file(params.variants_bed)
variants_ch = file(params.variants)
hapmap_ch = file(params.hapmap)
onekg_gwas_by_subtraction_reference = file("data/reference.1000G.maf.0.005.txt.gz").collect()

ld_ch = Channel.fromPath(params.ld_w_dir, type: 'dir').collect()

gene_chunk_size=1

log.info """=================================================
Estimate heritability v${workflow.manifest.version}"
================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Input directory']                          = params.inputdir
summary['LD score directory']                       = params.lddir
summary['SNP ref. file']                            = params.snpref
summary['ENSEMBL gtf']                              = params.gtf
summary['Log file']                                 = params.logfile
summary['eQTL P-thresh']                            = params.pthresh
summary['eQTL window']                              = params.eqtlwindow
summary['Output directory']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "================================================="


workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    if (params.per_cohort) {
        // Extract loci
        loci_extracted_ch = ExtractResultsPerCohort(input_parquet_ch, variant_reference_ch, variants_ch, genes_buffered_ch, params.cols, cohorts_ch.collect())
            .flatten()
            .map { file ->
                   def key = file.name.toString().tokenize('.').get(1)
                   return tuple(key, file) }
            groupTuple()
    } else {
        // Extract loci
        loci_extracted_ch = ExtractResults(input_parquet_ch, variant_reference_ch, variants_ch, genes_buffered_ch, params.cols)
            .flatten()
            .map { file ->
                   def key = file.name.toString().tokenize('.').get(1)
                   return tuple(key, file) }
            groupTuple()
    }

    // Split summary statistics in cis and trans regions
    results_ch = ProcessResults(loci_extracted_ch, variant_reference_ch, gene_reference_ch, inclusion_step_output_ch, cohorts_ch.collect())

    // Count the number of heritability variants for each gene
    heritability_snps_file_ch = CountHeritabilitySnps(gene_reference_ch, one_kg_bed_ch)

    // Heritability SNPs
    heritability_snps_cis = heritability_snps_file_ch.cis.splitCsv(header:false, sep:' ')
        .map { row -> tuple(row[1], "cis", row[0]) }
    heritability_snps_trans = heritability_snps_file_ch.trans.splitCsv(header:false, sep:' ')
        .map { row -> tuple(row[1], "trans", row[0]) }
    heritability_snps_gw = heritability_snps_file_ch.gw.splitCsv(header:false, sep:' ')
        .map { row -> tuple(row[1], "gw", row[0]) }

    ldsc_cis_in_ch = results_ch.cis.join(heritability_snps_cis, by:[0,1], remainder:false)
    ldsc_trans_in_ch = results_ch.trans.join(heritability_snps_trans, by:[0,1], remainder:false)
    ldsc_gw_in_ch = results_ch.gw.join(heritability_snps_gw, by:[0,1], remainder:false)

    // Process GWAS data
    process_gwas_ch = ProcessVuckovicGwasData(gwas_input_ch, hapmap_ch)

    // Run Heritability estimates
    ldsc_cis_output_ch = EstimateCisHeritabilityLdsc(
        ldsc_cis_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)

    ldsc_trans_output_ch = EstimateTransHeritabilityLdsc(
        ldsc_trans_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)

    ldsc_gw_output_ch = EstimateGwHeritabilityLdsc(
        ldsc_gw_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)

    EstimateHeritabilityLdscAllPairwise(
        process_gwas_ch.map { name, gws, file -> name }.collect(),
        process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)

    GwasBySubtraction(
        ldsc_gw_in_ch,
        process_gwas_ch.map { name, gws, file -> name }.collect(),
        process_gwas_ch.map { name, gws, file -> file }.collect(),
        ld_ch,
        onekg_gwas_by_subtraction_reference)

    // Process LDSC logs
    ldsc_cis_matrices_ch = ProcessCisLdscOutput(ldsc_cis_output_ch)
        .collectFile(name:'ldsc_table_cis.txt', skip: 1, keepHeader: true, storeDir: params.output)
    ldsc_trans_matrices_ch = ProcessTransLdscOutput(ldsc_trans_output_ch)
        .collectFile(name:'ldsc_table_trans.txt', skip: 1, keepHeader: true, storeDir: params.output)
    ldsc_gw_matrices_ch = ProcessGwLdscOutput(ldsc_gw_output_ch)
        .collectFile(name:'ldsc_table_gw.txt', skip: 1, keepHeader: true, storeDir: params.output)

    // Process LDSC stuff
    ProcessLdscDeleteVals(
        ldsc_trans_output_ch.map { name, gws, file, del -> name }.collect(),
        ldsc_trans_output_ch.map { name, gws, file, del -> del }.collect(),
        ldsc_trans_matrices_ch)

    // WriteOutRes(heritability_estimates.collectFile(name:'result.txt', sort: true, keepHeader: true))
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
