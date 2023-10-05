#!/bin/bash nextflow

process EstimateHeritability {

    tag "$gene"

    //publishDir "${params.outdir}", mode: 'copy', pattern: '*.png'

    input:
      tuple val(gene), path(input_ch), path(ld_ch), path(snpref), path(gtf), path(logfile)
      val(pthresh)
      val(eqtlwindow)

    output:
      path '*_h2.txt', emit: h2
      //path '*_h2.png', emit: plots

    shell:
    // Should first limit to the trans variants
    '''
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla !{baseDir}/bin/EstimateHeritability.R \
    !{gtf} \
    !{ld_ch} \
    !{snpref} \
    !{input_ch} \
    !{logfile} \
    !{pthresh} \
    !{eqtlwindow}
    '''
}

process CountHeritabilitySnps {

    input:
      path geneReference
      path oneKgBedFiles

    output:
      path "cis_gen_annot_M_5_50.txt", emit: cis
      path "trans_gen_annot_M_5_50.txt", emit: trans

    shell:
    '''
    gene_bed_files.py \
        --gene-ref !{geneReference} \
        --out-prefix genes

    bedtools intersect -a "genes.cis.bed" -b !{oneKgBedFiles} | \
        awk -F '\t' '{print $4}' | sort | uniq -c > "cis_gen_annot_M_5_50.txt"
    bedtools intersect -v -a "genes_trans.bed" -b !{oneKgBedFiles} | \
            awk -F '\t' '{print $4}' | sort | uniq -c > "trans_gen_annot_M_5_50.txt"
    '''
}

process EstimateHeritabilityLdsc {
    container 'quay.io/cawarmerdam/ldsc:v0.3'
    errorStrategy = 'ignore'
    executor = 'local'

    input:
      tuple val(name_a), val(annot_a), path(sumstats_a), val(name_b), val(annot_b), path(sumstats_b)
      path ld_ch

    output:
      tuple val(name_a), val(annot_a), val(name_b), val(annot_b), path('*_rg_*.log')

    shell:
    // Should first limit to the trans variants
    '''
    /ldsc/ldsc.py \
    --rg !{sumstats_a},!{sumstats_b} \
    --ref-ld-chr !{ld_ch}/ \
    --w-ld-chr !{ld_ch}/ \
    --chisq-max 10000 \
    --out !{name_a.replaceAll("\\s","_")}_rg_!{name_b.replaceAll("\\s","_")}
    '''
}

process EstimateHeritabilityGenomicSem {
    errorStrategy = 'ignore'

    input:
      tuple val(gene), val(annot), path(sumstats)
      tuple val(name), path(gwas)
      path ld_ch

    output:
      tuple val(gene), val(annot), path('*_rg.log')

    shell:
    // Should first limit to the trans variants
    '''
    genomic_sem_ldsc.R \
    --rg !{sumstats_a} !{gwas.join(" ")} \
    --ref-ld-chr !{ld_ch}/ \
    --w-ld-chr !{ld_ch}/ \
    --names "!{gene}" "!{name.join("\" \"")}"
    '''
}

process ProcessLdscOutput {
    publishDir "${params.output}", mode: 'copy', pattern: '*_h2.txt'

    input:
      tuple val(name_a), val(annot_a), val(name_b), val(annot_b), path(ldsc_output)

    output:
      path '*_h2.txt'

    shell:
    // Should first limit to the trans variants
    '''
    process_ldsc_output.R !{gene} !{annot} !{ldsc_output}
    '''
}
