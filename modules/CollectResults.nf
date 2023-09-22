#!/bin/bash nextflow


process ExtractResults {
    scratch true

    input:
        path input
        path variantReference
        path variants
        val genes
        val cols

    output:
        path "extracted*out.csv"

    shell:
        variants_arg = (variants.name != 'NO_FILE') ? "--variants-file ${variants}" : ""
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        prefix = (genes.size() == 1) ? genes.collect { "extracted.$it" }.join("") : "extracted"
        '''
        mkdir tmp_eqtls
        echo "!{phenotypes_formatted}" > file_matches.txt

        while read gene; do
          cp -r "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            --variant-reference !{variantReference} \
            --genes !{genes.join(' ')} \
            !{variants_arg} \
            --cols '!{cols}' \
            --bed-file !{locus} \
            --output-prefix extracted

        rm -r tmp_eqtls
        '''
}


process ProcessResults {

    input:
        tuple val(locus_string), path(files, stageAs: "locus_*.csv")
        path variantReference
        path geneReference

    output:
        tuple val(locus_string), val("cis"), path("annotated.${locus_string}_cis.csv.gz"), emit: cis
        tuple val(locus_string), val("trans"), path("annotated.${locus_string}_trans.csv.gz"), emit: trans

    script:
        """
        head -n 1 ${files[0]} > concatenated.${locus_string}.csv
        tail -n +2 ${files.join(' ')} >> concatenated.${locus_string}.csv

        annotate_loci.py \
            --input-file concatenated.${locus_string}.csv \
            --variant-reference ${variantReference} \
            --gene-gff ${geneReference} \
            --out-prefix annotated.${locus_string}
        """
}