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

    n_total=$(cat !{oneKgBedFiles} | wc -l | awk '{print $1}')

    bedtools intersect -a "genes.cis.bed" -b !{oneKgBedFiles} | \
        awk -F '\t' '{print $4}' | sort | uniq -c | sed 's/^ *//' > "cis_gen_annot_M_5_50.txt"
    bedtools intersect -a "genes.trans.bed" -b !{oneKgBedFiles} | \
            awk -F '\t' '{print $4}' | sort | uniq -c | sed 's/^ *//' | awk -v tot=$n_total '{print (tot - $1),$2}' > "trans_gen_annot_M_5_50.txt"
    '''
}

process EstimateHeritabilityLdsc2 {
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
    --print-cov \
    --print-delete-vals \
    --out !{name_a.replaceAll("\\s","_")}_rg_!{name_b.replaceAll("\\s","_")}
    '''
}

process EstimateHeritabilityLdsc {
    container 'quay.io/cawarmerdam/ldsc:v0.3'
    tag "ldsc_${annot}_${gene}"
    errorStrategy = 'ignore'

    input:
      tuple val(gene), val(annot), path(sumstats), val(m_5_50)
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
    --M !{m_5_50} \
    --out !{gene}_rg \
    --print-cov \
    --print-delete-vals \
    '''
}

process EstimateHeritabilityLdscAllPairwise {
    container 'quay.io/cawarmerdam/ldsc:v0.3'
    tag "ldsc_${annot}_${gene}"
    errorStrategy = 'ignore'

    input:
      val name
      path gwas
      path ld_ch

    output:
      val name
      path 'rg_*_*.log'

    shell:
    // Should first limit to the trans variants
    '''
    f=( !{gwas.join(" ")} )
    n=( !{name.join(" ")} )

    function join_by { local IFS="$1"; shift; echo "$*"; }

    for ((i = 0; i+1 < ${#f[@]}; i++)); do

        /ldsc/ldsc.py \
        --rg ${f[i]},$(join_by , "${f[@]:i+1}") \
        --ref-ld-chr !{ld_ch}/ \
        --w-ld-chr !{ld_ch}/ \
        --out rg_${i}_${i+1}
    done
    '''
}


process EstimateHeritabilityGenomicSem {
    errorStrategy = 'ignore'

    input:
      tuple val(gene), val(annot), path(sumstats)
      tuple val(name), val(annot_b), path(gwas)
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
      tuple val(gene), val(annot), path(ldsc_output)

    output:
      path '*_h2.txt'

    shell:
    // Should first limit to the trans variants
    '''
    process_ldsc_output.R !{gene} !{annot} !{ldsc_output}
    '''
}
