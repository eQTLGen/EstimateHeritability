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


process EstimateHeritabilityLdsc {
    container 'quay.io/cawarmerdam/ldsc:v0.1'
    tag "ldsc_${annot}_${gene}"
    errorStrategy = 'ignore'

    input:
      tuple val(gene), val(annot), path(sumstats)
      path ld_ch
      path ldsc_source

    output:
      tuple val(gene), path('*_h2.log')

    shell:
    // Should first limit to the trans variants
    '''
    !{ldsc_source}/ldsc.py \
    --h2 !{sumstats} \
    --ref-ld-chr !{ld_ch}/ \
    --w-ld-chr !{ld_ch}/ \
    --chisq-max 10000 \
    --out !{gene}_h2
    '''
}

process ProcessLdscOutput {

    input:
      tuple val(gene), path(ldsc_output)

    output:
      path '*_h2.txt'

    shell:
    // Should first limit to the trans variants
    '''
    sed -n '/Total Observed scale h2/,$p' !{ldsc_output} > h2.txt
    Rscript process_ldsc_output.R !{gene} h2.txt
    '''
}