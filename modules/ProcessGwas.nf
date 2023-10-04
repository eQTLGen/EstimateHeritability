#!/bin/bash nextflow

process ProcessVuckovicGwasData {
    container 'quay.io/cawarmerdam/ldsc:v0.3'

    input:
        tuple val(name), path(sumstats), val(sample_size)
        path snplist

    output:
        tuple val(name), val("ctc"), path("${name.replaceAll("\\s","_")}.processed.sumstats.gz")

    script:
        """
        /ldsc/munge_sumstats.py \
        --sumstats ${sumstats} \
        --N ${sample_size} \
        --out ${name.replaceAll("\\s","_")}.processed \
        --merge-alleles ${snplist} \
        --snp variant_id \
        --chunksize 500000
        """
}
