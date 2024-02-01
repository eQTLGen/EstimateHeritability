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
            --output-prefix extracted.!{genes.join('_')}

        rm -r tmp_eqtls
        '''
}


process ExtractResultsPerCohort {
    scratch true

    input:
        path input
        path variantReference
        path variants
        val genes
        val cols
        val cohort

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
            --cohort !{cohort.join(' ')} \
            --output-prefix extracted.!{genes.join('_')}

        rm -r tmp_eqtls
        '''
}


process ProcessResults {

    input:
        tuple val(locus_string), path(files, stageAs: "locus_*.csv")
        path variantReference
        path geneReference
        path inclusionDir
        val cohorts
        val isqThreshold

    output:
        tuple val(locus_string), val("cis"), path("annotated.${locus_string}_cis.csv.gz"), emit: cis, optional: true
        tuple val(locus_string), val("trans"), path("annotated.${locus_string}_trans.csv.gz"), emit: trans, optional: true
        tuple val(locus_string), val("gw"), path("annotated.${locus_string}_gw.csv.gz"), emit: gw, optional: true
        path "annotated.${locus_string}_passed_variants.csv", emit: variants

    script:
        """
        head -n 1 ${files[0]} > concatenated.${locus_string}.csv
        tail -n +2 ${files.join(' ')} >> concatenated.${locus_string}.csv

        annotate_loci.py \
            --input-file concatenated.${locus_string}.csv \
            --cohorts ${cohorts.join(' ')} \
            --inclusion-path ${inclusionDir} \
            --variant-reference ${variantReference} \
            --gene-ref ${geneReference} \
            --out-prefix annotated.${locus_string} \
            --i2-threshold ${isqThreshold}
        """
}


process PrepareHeritabilityEstimation {
    scratch false

    input:
        path input
        path variantReference
        path variants
        path geneReference
        path inclusionDir
        val genes
        val cohorts
        val isqThreshold
        path ld_ch
        path frqfile_ch

    output:
        path "*.sumstats.*.csv.gz", emit: sumstats, optional:true
        path "*.lead_effects.csv", emit: leads
        path "heritability_snps.*.csv", emit: variants

    shell:
        variants_arg = (variants.name != 'NO_FILE') ? "--variants-file ${variants}" : ""
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
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
            --cols '+z_score,+p_value' \
            --output-prefix extracted

        annotate_loci.py \
            --input-file extracted.out.csv \
            --cohorts !{cohorts.join(' ')} \
            --inclusion-path !{inclusionDir} \
            --variant-reference !{variantReference} \
            !{variants_arg} \
            --gene-ref !{geneReference} \
            --out-prefix annotated.!{genes.join("_")} \
            --i2-threshold !{isqThreshold}

        count_heritability_snps.py \
            --bed-file-incl "cis.bed"
            --bed-file-excl "trans.bed" "polygenic.bed" \
            --annot-chr !{ld_ch}/baselineLD. \
            --frqfile-chr !{frqfile_ch}/1000G.EUR.hg38. \
            --out-prefix counts

        rm -r tmp_eqtls
        '''
}
