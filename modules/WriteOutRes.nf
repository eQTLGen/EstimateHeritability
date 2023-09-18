#!/bin/bash nextflow 


process WriteOutRes {

    publishDir "${params.outdir}", mode: 'copy'

    input:
      path "result.txt"

    output:
      path 'result.txt', emit: h2

    script:
        """
        echo "Write out results"
        """
}
