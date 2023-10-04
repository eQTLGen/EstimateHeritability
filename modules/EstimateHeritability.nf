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


process ExtractGenomicCoordinates {

    input:
      path geneReference

    output:
      path bedFile

    script:
    """
    gene_bed_files.py \
        --gene-ref ${geneReference}
    """
}


process CalculateHeritabilitySnps {

    input:
      path bedFile
      path oneKgBedFiles

    output:
      path "cis_trans_gen_annot_M_5_50.txt"

    shell:
    // Presume filtering of 1000Kg files to variants with MAF > 5% is already done
    // Per row in the bedFile, count the number of intersecting variants
    '''
    plink2 --bfile
    '''
}

process EstimateHeritabilityLdsc {
    container 'quay.io/cawarmerdam/ldsc:v0.3'
    tag "ldsc_${annot}_${gene}"
    errorStrategy = 'ignore'

    input:
      tuple val(name_a), val(annot_a), path(sumstats_a), val(name_b), val(annot_b), path(sumstats_b), val(m_5_50)
      path gwas
      path ld_ch

    output:
      tuple val(gene), val(annot), path('*_rg.log')

    shell:
    // Should first limit to the trans variants
    '''
    /ldsc/ldsc.py \
    --rg !{sumstats},!{gwas.join(",")} \
    --ref-ld-chr !{ld_ch}/ \
    --w-ld-chr !{ld_ch}/ \
    --chisq-max 10000 \
    --M !{m_5_50}
    --out !{gene}_rg
    '''
}

process ProcessLdscOutput {
    publishDir "${params.output}", mode: 'copy', pattern: '*_h2.txt'

    input:
      tuple val(gene), val(annot), path(ldsc_output)

    output:
      path '*_h2.txt'

    shell:
    // Should first limit to the trans variants
    '''
    process_ldsc_output.R !{gene} !{annot} !{ldsc_output}
    '''
}
