#!/bin/bash nextflow 


process WriteOutRes {

    publishDir "${params.outdir}", mode: 'copy'

    input:
      tuple val(gene_id), path(ldsc_log)

    output:
      path 'result.txt', emit: h2

    script:
        """
        echo "Write out results"
        """
}
