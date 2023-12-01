#!/bin/bash nextflow


/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

// import modules
include { GwasBySubtraction; GwasBySubtraction as GwasBySubtractionTrans; EstimateTransHeritabilityLdsc; EstimateTransHeritabilityLdsc as EstimateGwHeritabilityLdsc; EstimateCisHeritabilityLdsc; ProcessLdscOutput as ProcessTransLdscOutput; ProcessLdscOutput as ProcessCisLdscOutput; ProcessLdscOutput as ProcessGwLdscOutput; CountHeritabilitySnps; EstimateHeritabilityLdscAllPairwise; ProcessLdscDeleteVals; ProcessLdscDeleteVals as ProcessLdscDeleteValsGw } from './modules/EstimateHeritability'
include { WriteOutRes } from './modules/WriteOutRes'
include { ExtractResults; ExtractResultsPerCohort; ProcessResults; LoadResultsAnnotated } from './modules/CollectResults.nf'
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

//cohorts_ch = Channel
//   .fromList( ["GTEx_2017-06-05_v8_EUR"] )

inclusion_step_output_ch = file(params.inclusion_step_output)
one_kg_bed_ch = file(params.variants_bed)
variants_ch = file(params.variants)
hapmap_ch = file(params.hapmap)
i_squared_threshold = 100
onekg_gwas_by_subtraction_reference = Channel.fromPath("data/reference.1000G.maf.0.005.txt.gz").collect()

Channel.fromPath(params.maf_table).collect().set { maf_table_ch }

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

    // if (params.per_cohort) {
    //     // Extract loci
    //     results_ch = ExtractResultsPerCohort(input_parquet_ch, variant_reference_ch, variants_ch, genes_buffered_ch, params.cols, cohorts_ch.collect())
    //         .flatten()
    //         .map { file ->
    //                def key = file.name.toString().tokenize('.').get(1)
    //                return tuple(key, file) }
    //         groupTuple()
    // } else {
    //     // Extract loci
    //     results_ch = ExtractResults(input_parquet_ch, variant_reference_ch, variants_ch, genes_buffered_ch, params.cols)
    //         .flatten()
    //         .map { file ->
    //                def key = file.name.toString().tokenize('.').get(1)
    //                return tuple(key, file) }
    //         groupTuple()
    // }

    LoadResultsAnnotated(
        input_parquet_ch, variant_reference_ch, variants_ch, gene_reference_ch, inclusion_step_output_ch, maf_table_ch,
        genes_buffered_ch, cohorts_ch.collect(), i_squared_threshold)

    results_ch = LoadResultsAnnotated.out
        .flatten()
        .map { file ->
               def gene = file.name.toString().tokenize('.').get(2)
               def annot = file.name.toString().tokenize('.').get(3)
               def cis_trans_gw = annot.tokenize('_').get(0)
               return tuple(gene, cis_trans_gw, annot, file) }
        .groupTuple()
        .view()

    // Split summary statistics in cis and trans regions
    // results_ch = ProcessResults(loci_extracted_ch, variant_reference_ch, gene_reference_ch, inclusion_step_output_ch, cohorts_ch.collect(), i_squared_threshold)

    // List number of variants per gene
    LoadResultsAnnotated.out.variants
        .collectFile(keepHeader: true, skip: 1, name: "variants_per_gene.txt", storeDir: params.output)

    //// Count the number of heritability variants for each gene
    //heritability_snps_file_ch = CountHeritabilitySnps(gene_reference_ch, one_kg_bed_ch)
//
    //// Heritability SNPs
    //heritability_snps = heritability_snps_file_ch.splitCsv(header:false, sep:' ')
    //    .flatten()
    //    .map { file ->
    //           def gene = file.name.toString().tokenize('.').get(2)
    //           def cis_trans_gw = file.name.toString().tokenize('.').get(3)
    //           return tuple(gene, cis_trans_gw, file) }
    //    .groupTuple()
//
    //ldsc_in_ch = results_ch.join(heritability_snps, by:[0,1], remainder:false)
//
    //// Process GWAS data
    //process_gwas_ch = ProcessVuckovicGwasData(gwas_input_ch, hapmap_ch)
//
    //// Run Heritability estimates
    //ldsc_cis_output_ch = EstimateCisHeritabilityLdsc(
    //    ldsc_cis_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)
//
    //ldsc_trans_output_ch = EstimateTransHeritabilityLdsc(
    //    ldsc_trans_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)
//
    //ldsc_gw_output_ch = EstimateGwHeritabilityLdsc(
    //    ldsc_gw_in_ch, process_gwas_ch.map { name, gws, file -> file }.collect(), ld_ch)
//
    //GwasBySubtraction(
    //    ldsc_gw_in_ch,
    //    gwas_input_ch.map { name, gws, n -> gws }.collect(),
    //    process_gwas_ch.map { name, gws, file -> name }.collect(),
    //    process_gwas_ch.map { name, gws, file -> file }.collect(),
    //    ld_ch, onekg_gwas_by_subtraction_reference)
    //    .collectFile(name:'latent_factors_cis_comp.txt', skip: 1, keepHeader: true, storeDir: params.output)
    //GwasBySubtractionTrans(
    //    ldsc_trans_in_ch,
    //    gwas_input_ch.map { name, gws, n -> gws }.collect(),
    //    process_gwas_ch.map { name, gws, file -> name }.collect(),
    //    process_gwas_ch.map { name, gws, file -> file }.collect(),
    //    ld_ch, onekg_gwas_by_subtraction_reference)
    //    .collectFile(name:'latent_factors_trans_comp.txt', skip: 1, keepHeader: true, storeDir: params.output)
//
    //// Process LDSC logs
    //ldsc_cis_matrices_ch = ProcessCisLdscOutput(ldsc_cis_output_ch)
    //    .collectFile(name:'ldsc_table_cis.txt', skip: 1, keepHeader: true, storeDir: params.output)
    //ldsc_trans_matrices_ch = ProcessTransLdscOutput(ldsc_trans_output_ch)
    //    .collectFile(name:'ldsc_table_trans.txt', skip: 1, keepHeader: true, storeDir: params.output)
    //ldsc_gw_matrices_ch = ProcessGwLdscOutput(ldsc_gw_output_ch)
    //    .collectFile(name:'ldsc_table_gw.txt', skip: 1, keepHeader: true, storeDir: params.output)
//
    //// Process LDSC stuff
    //ProcessLdscDeleteVals(
    //    ldsc_trans_output_ch.map { name, gws, file, del -> name }.collect(),
    //    ldsc_trans_output_ch.map { name, gws, file, del -> del }.collect(),
    //    ldsc_trans_matrices_ch, "trans")
//
    //ProcessLdscDeleteValsGw(
    //    ldsc_gw_output_ch.map { name, gws, file, del -> name }.collect(),
    //    ldsc_gw_output_ch.map { name, gws, file, del -> del }.collect(),
    //    ldsc_gw_matrices_ch, "gw")

    // WriteOutRes(heritability_estimates.collectFile(name:'result.txt', sort: true, keepHeader: true))
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
