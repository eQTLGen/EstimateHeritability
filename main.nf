#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { EstimateTransHeritabilityLdsc; EstimateTransHeritabilityLdsc as EstimatePolyHeritabilityLdsc; EstimateCisHeritabilityLdsc; ProcessLdscOutput as ProcessTransLdscOutput; ProcessLdscOutput as ProcessCisLdscOutput; ProcessLdscOutput as ProcessGwLdscOutput; CountHeritabilitySnps; EstimateHeritabilityLdscAllPairwise; ProcessLdscDeleteVals; ProcessLdscDeleteVals as ProcessLdscDeleteValsGw } from './modules/EstimateHeritability'
include { WriteOutRes } from './modules/WriteOutRes'
include { PrepareHeritabilityEstimation } from './modules/CollectResults.nf'
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
params.cols = '+z_score,+p_value'
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

inclusion_step_output_ch = file(params.inclusion_step_output)
one_kg_bed_ch = file(params.variants_bed)
variants_ch = file(params.variants)
hapmap_ch = file(params.hapmap)
i_squared_threshold = 100
onekg_gwas_by_subtraction_reference = Channel.fromPath("data/reference.1000G.maf.0.005.txt.gz").collect()

ld_ch = Channel.fromPath(params.ld_w_dir, type: 'dir').collect()
frqfile_ch = Channel.fromPath(params.frqfile_dir, type: 'dir').collect()
weights_ch = Channel.fromPath(params.weights_dir, type: 'dir').collect()

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

    PrepareHeritabilityEstimation(
        input_parquet_ch, variant_reference_ch, variants_ch, gene_reference_ch, inclusion_step_output_ch,
        genes_buffered_ch, cohorts_ch.collect(), i_squared_threshold, ld_ch, frqfile_ch)

    polygenic_ch = PrepareHeritabilityEstimation.out.sumstats_transpolygenic
        .flatten()
        .map { file ->
               def gene = file.name.toString().tokenize('.').get(0)
               return tuple(gene, file) }

    cis_ch = PrepareHeritabilityEstimation.out.sumstats_cis
        .flatten()
        .map { file ->
               def gene = file.name.toString().tokenize('.').get(0)
               return tuple(gene, file) }

    trans_ch = PrepareHeritabilityEstimation.out.sumstats_trans
        .flatten()
        .map { file ->
               def gene = file.name.toString().tokenize('.').get(0)
               return tuple(gene, file) }

    // List number of variants per gene
    lead_effects_ch = PrepareHeritabilityEstimation.out.leads
        .collectFile(keepHeader: true, skip: 1, name: "lead_effect_variants.txt", storeDir: params.output)

    // Heritability SNPs
    ldsc_cis_in_ch = PrepareHeritabilityEstimation.out.cis_variants.collectFile()
        .splitCsv(header:false, sep:'\t')
        .map { row -> return tuple(row[0], row[1]) }
        .join(cis_ch, by:[0], remainder:false)

    ldsc_trans_in_ch = PrepareHeritabilityEstimation.out.trans_variants.collectFile()
        .splitCsv(header:false, sep:'\t')
        .map { row -> return tuple(row[0], row[1]) }
        .join(trans_ch, by:[0], remainder:false)

    ldsc_polygenic_in_ch = PrepareHeritabilityEstimation.out.transpolygenic_variants.collectFile()
        .splitCsv(header:false, sep:'\t')
        .map { row -> return tuple(row[0], row[1]) }
        .join(polygenic_ch, by:[0], remainder:false)

    // Run Heritability estimates
    ldsc_cis_output_ch = EstimateCisHeritabilityLdsc(
        ldsc_cis_in_ch, ld_ch, frqfile_ch, weights_ch)

    ldsc_trans_output_ch = EstimateTransHeritabilityLdsc(
        ldsc_trans_in_ch, ld_ch, frqfile_ch, weights_ch)

    ldsc_polygenic_output_ch = EstimatePolyHeritabilityLdsc(
        ldsc_polygenic_in_ch, ld_ch, frqfile_ch, weights_ch)

    // Process LDSC logs
    ldsc_cis_matrices_ch = ProcessCisLdscOutput(ldsc_cis_output_ch)
        .collectFile(name:'ldsc_table_cis.txt', skip: 1, keepHeader: true, storeDir: params.output)
    ldsc_trans_matrices_ch = ProcessTransLdscOutput(ldsc_trans_output_ch)
        .collectFile(name:'ldsc_table_trans.txt', skip: 1, keepHeader: true, storeDir: params.output)
    ldsc_polygenic_matrices_ch = ProcessGwLdscOutput(ldsc_polygenic_output_ch)
        .collectFile(name:'ldsc_table_polygenic.txt', skip: 1, keepHeader: true, storeDir: params.output)

    // Process LDSC stuff
    ProcessLdscDeleteVals(
        ldsc_trans_output_ch.map { name, gws, file, del -> name }.collect(),
        ldsc_trans_output_ch.map { name, gws, file, del -> del }.collect(),
        ldsc_trans_matrices_ch, "trans")

    ProcessLdscDeleteValsGw(
        ldsc_polygenic_output_ch.map { name, gws, file, del -> name }.collect(),
        ldsc_polygenic_output_ch.map { name, gws, file, del -> del }.collect(),
        ldsc_polygenic_matrices_ch, "polygenic")

    // WriteOutRes(heritability_estimates.collectFile(name:'result.txt', sort: true, keepHeader: true))
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
