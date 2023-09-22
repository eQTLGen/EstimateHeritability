#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { EstimateHeritability; EstimateHeritabilityLdsc } from './modules/EstimateHeritability'
include { WriteOutRes } from './modules/WriteOutRes'
include { ExtractResults; ProcessResults } from './modules/CollectResults.nf'

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

params.maf_table = 'NO_FILE'
params.bed = 'NO_FILE'
params.variants = 'NO_FILE'
params.inclusion_step_output = 'NO_FILE'
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

bed_file_ch = file(params.bed)
variants_ch = file(params.variants)

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


// Process input file paths

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }

ld_ch = Channel.fromPath(params.ld_w_dir, type: 'dir').collect()

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    // Extract loci
    loci_extracted_ch = ExtractResults(input_parquet_ch, variant_reference_ch, variants_ch, genes_buffered_ch, params.cols)
        .flatten()
        .map { file ->
               def key = file.name.toString().tokenize('.').get(1)
               return tuple(key, file) }
        groupTuple()
    }

    // Split summary statistics in cis and trans regions
    results_ch = ProcessResults(loci_extracted_ch, variant_reference_ch, gene_reference_ch)

    // Combine results in a single channel
    results_ch_concatenated = results_ch.cis.concat(results_ch.trans)

    // Run Heritability estimates
    heritability_estimates = EstimateHeritability(results_ch_concatenated, ld_ch)

    WriteOutRes(heritability_estimates.out.collectFile(name:'result.txt', sort: true, keepHeader: true))
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}