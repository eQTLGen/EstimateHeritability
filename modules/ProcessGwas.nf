#!/bin/bash nextflow

process ProcessVuckovicGwasData {
    container 'quay.io/cawarmerdam/ldsc:v0.3'

    input:
        tuple val(name), path(sumstats), val(sample_size)
        path snplist

    output:
        tuple val(name), path("processed.${name}_cis.csv.gz")

    script:
        """
        ldsc/munge_sumstats.py \
        --sumstats ${sumstats} \
        --N ${sample_size} \
        --out ${name}_processed \
        --merge-alleles ${snplist}
        """
}