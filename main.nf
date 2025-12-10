#!/usr/bin/env nextflow

log.info """\

        ================================================================
        NF-varmit - Variant calling for whole-genome resequence data, iteratively
        https://github.com/MozesBlom/nf-variant
        Authors: Mozes P.K. Blom ; Jon J. Hoffman
        ================================================================
        |indivfile    : ${params.indivs_file}
        |chromos      : ${params.chromos_file}
        |reference    : ${params.ref_file}
        |readsdir     : ${params.readsdir}
        |outputdir    : ${params.outputdir}
        |
        |Sex chromos  : ${params.sex_chromos}
        |
        |Rounds of mapping/consensus : ${params.n_rounds ?: 4}
        |
        |Obtain coverage stats?         : ${params.calc_coverage}
        |Call variants by individual?   : ${params.indiv_var_call}
        |Call consensus sequences?      : ${params.call_consensus}
        |Calcuate missing data?         : ${params.calc_missing_data}
        |
        |Variant calling & filtering
        |---------------------------
        |Filter indels?                 : ${params.filt_indels}
        |Mask heterozygous positions?   : ${params.mask_hets}
        |Mask low/high coverage sites?  : ${params.mask_cov}
        |Min depth per position?        : ${params.mask_min_cov}
        |Max depth per position?        : ${params.mask_max_cov}
        |
        ================================================================
        """
        .stripIndent()

/*
====================================================
~ ~ ~ > *  Input channels and file check  * < ~ ~ ~
====================================================
*/

indivs    = file(params.indivs_file).readLines()
indivs_ch = Channel.fromList(indivs)

chromos    = file(params.chromos_file).readLines()
chromos_ch = Channel.fromList(chromos)

/*
 * Reads channel (paired-end)
 */
params.readsdir      = params.readsdir      ?: params.inputdir
params.reads_suffix1 = params.reads_suffix1 ?: '_R1.fastq.gz'
params.reads_suffix2 = params.reads_suffix2 ?: '_R2.fastq.gz'

reads_ch = indivs_ch.map { indiv ->
    def r1 = file("${params.readsdir}/${indiv}${params.reads_suffix1}")
    def r2 = file("${params.readsdir}/${indiv}${params.reads_suffix2}")
    if( !r1.exists() ) error "Missing R1 reads for ${indiv}: ${r1}"
    if( !r2.exists() ) error "Missing R2 reads for ${indiv}: ${r2}"
    tuple(indiv, r1, r2)
}

/*
 * Reference file
 */
def ref_file_obj = file(params.ref_file)
if( !ref_file_obj.exists() ) {
    error "No reference sequence found from: ${params.ref_file}"
}
ref_ch = Channel.value(ref_file_obj)

/*
 * Number of iterative rounds (1–4)
 */
params.n_rounds = (params.n_rounds ? params.n_rounds.toInteger() : 4)
if( params.n_rounds < 1 || params.n_rounds > 4 ) {
    error "params.n_rounds must be between 1 and 4 (got: ${params.n_rounds})"
}

/*
===============================
~ ~ ~ > *  Processes  * < ~ ~ ~
===============================
*/

process index_bam {

    tag "Index bam file"

    input:
    tuple val(indiv), path(indiv_bam)

    output:
    tuple val(indiv), path(indiv_bam), path("${indiv_bam}.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${indiv_bam}
    """
}

/*
 * BWA-MEM2 mapping processes (per round)
 */

process bwa_map_R1 {

    label 'Mapping'
    tag "BWA-MEM2 mapping (round1)"
    publishDir "${params.outputdir}/00.bam/round1", mode: 'copy'

    input:
    tuple val(individual), path(read1), path(read2), path(reference)

    output:
    tuple val(individual), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    samtools faidx ${reference}
    bwa-mem2 index ${reference}

    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} | \
      samtools sort -@ ${task.cpus} -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}

process bwa_map_R2 {

    label 'Mapping'
    tag "BWA-MEM2 mapping (round2)"
    publishDir "${params.outputdir}/00.bam/round2", mode: 'copy'

    input:
    tuple val(individual), path(read1), path(read2), path(reference)

    output:
    tuple val(individual), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    samtools faidx ${reference}
    bwa-mem2 index ${reference}

    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} | \
      samtools sort -@ ${task.cpus} -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}

process bwa_map_R3 {

    label 'Mapping'
    tag "BWA-MEM2 mapping (round3)"
    publishDir "${params.outputdir}/00.bam/round3", mode: 'copy'

    input:
    tuple val(individual), path(read1), path(read2), path(reference)

    output:
    tuple val(individual), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    samtools faidx ${reference}
    bwa-mem2 index ${reference}

    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} | \
      samtools sort -@ ${task.cpus} -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}

process bwa_map_R4 {

    label 'Mapping'
    tag "BWA-MEM2 mapping (round4)"
    publishDir "${params.outputdir}/00.bam/round4", mode: 'copy'

    input:
    tuple val(individual), path(read1), path(read2), path(reference)

    output:
    tuple val(individual), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    samtools faidx ${reference}
    bwa-mem2 index ${reference}

    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} | \
      samtools sort -@ ${task.cpus} -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}

/*
 * Coverage processes (used on Round 1 only)
 */

process cov_estimate {

    tag "Estimate coverage"

    input:
    tuple val(indiv), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(indiv), path("${chromo}.tsv")

    script:
    def out_fn = "${chromo}.tsv"
    """
    samtools coverage --region ${chromo} --output ${out_fn} ${indiv_bam}
    """
}

process cov_summary_INDIV {

    tag "Coverage summary per individual"
    publishDir "${params.outputdir}/00.coverage", mode:'copy'

    input:
    tuple val(indiv), path(chromo_cov_tsv_list)

    output:
    path("${indiv}_coverage.tsv")

    script:
    """
    00_cov_stats_INDIV.py \
      -c ${chromo_cov_tsv_list} \
      -i ${indiv}
    """
}

process cov_summary_ALL {

    tag "Coverage summary for all indivs"
    publishDir "${params.outputdir}/00.coverage", mode:'copy'

    input:
    path(indivs_cov_tsv_list)

    output:
    file('*')

    script:
    """
    01_cov_stats_SUMMARY.py \
      -c ${indivs_cov_tsv_list} \
      -s ${params.sex_chromos} \
      -p ${params.plots_per_row}
    """
}

/*
 * Variant calling per round
 */

process call_variants_CHROMO_R1 {

    tag "Variant calling per chromosome (round1)"
    publishDir path: { "${params.outputdir}/01.variants/round1/${individual}" }, mode:'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt.vcf.gz")

    script:
    """
    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam} | \
      vcffilter -f "QUAL < 20" -f "( AB > 0 ) & ( AB < 0.2 )" --invert --or | \
      vcfallelicprimitives -k -g | \
      bgzip -c > ${chromo}_vars_filt.vcf.gz
    """
}

process call_variants_CHROMO_R2 {

    tag "Variant calling per chromosome (round2)"
    publishDir path: { "${params.outputdir}/01.variants/round2/${individual}" }, mode:'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt.vcf.gz")

    script:
    """
    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam} | \
      vcffilter -f "QUAL < 20" -f "( AB > 0 ) & ( AB < 0.2 )" --invert --or | \
      vcfallelicprimitives -k -g | \
      bgzip -c > ${chromo}_vars_filt.vcf.gz
    """
}

process call_variants_CHROMO_R3 {

    tag "Variant calling per chromosome (round3)"
    publishDir path: { "${params.outputdir}/01.variants/round3/${individual}" }, mode:'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt.vcf.gz")

    script:
    """
    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam} | \
      vcffilter -f "QUAL < 20" -f "( AB > 0 ) & ( AB < 0.2 )" --invert --or | \
      vcfallelicprimitives -k -g | \
      bgzip -c > ${chromo}_vars_filt.vcf.gz
    """
}

process call_variants_CHROMO_R4 {

    tag "Variant calling per chromosome (round4)"
    publishDir path: { "${params.outputdir}/01.variants/round4/${individual}" }, mode:'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt.vcf.gz")

    script:
    """
    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam} | \
      vcffilter -f "QUAL < 20" -f "( AB > 0 ) & ( AB < 0.2 )" --invert --or | \
      vcfallelicprimitives -k -g | \
      bgzip -c > ${chromo}_vars_filt.vcf.gz
    """
}

/*
 * Remove indels per round
 */

process remove_indels_R1 {

    tag "Remove indels (round1)"
    publishDir path: { "${params.outputdir}/01.variants/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt_indels.vcf.gz")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
      bgzip -c > ${chromo}_vars_filt_indels.vcf.gz
    """
}

process remove_indels_R2 {

    tag "Remove indels (round2)"
    publishDir path: { "${params.outputdir}/01.variants/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt_indels.vcf.gz")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
      bgzip -c > ${chromo}_vars_filt_indels.vcf.gz
    """
}

process remove_indels_R3 {

    tag "Remove indels (round3)"
    publishDir path: { "${params.outputdir}/01.variants/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt_indels.vcf.gz")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
      bgzip -c > ${chromo}_vars_filt_indels.vcf.gz
    """
}

process remove_indels_R4 {

    tag "Remove indels (round4)"
    publishDir path: { "${params.outputdir}/01.variants/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt_indels.vcf.gz")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
      bgzip -c > ${chromo}_vars_filt_indels.vcf.gz
    """
}

/*
 * Het mask per round
 */

process mask_hets_R1 {

    tag "Generate mask for het sites (round1)"
    publishDir path: { "${params.outputdir}/01.variants/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_hets.tsv")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f '( AF < 1 ) & ( AB < 0.8 )' | \
      cut -f 1,2 > ${chromo}_hets.tsv
    """
}

process mask_hets_R2 {

    tag "Generate mask for het sites (round2)"
    publishDir path: { "${params.outputdir}/01.variants/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_hets.tsv")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f '( AF < 1 ) & ( AB < 0.8 )' | \
      cut -f 1,2 > ${chromo}_hets.tsv
    """
}

process mask_hets_R3 {

    tag "Generate mask for het sites (round3)"
    publishDir path: { "${params.outputdir}/01.variants/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_hets.tsv")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f '( AF < 1 ) & ( AB < 0.8 )' | \
      cut -f 1,2 > ${chromo}_hets.tsv
    """
}

process mask_hets_R4 {

    tag "Generate mask for het sites (round4)"
    publishDir path: { "${params.outputdir}/01.variants/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_hets.tsv")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f '( AF < 1 ) & ( AB < 0.8 )' | \
      cut -f 1,2 > ${chromo}_hets.tsv
    """
}

/*
 * Coverage mask per round
 */

process mask_cov_R1 {

    tag "Generate mask for low or excess coverage sites (round1)"
    publishDir path: { "${params.outputdir}/01.variants/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}

process mask_cov_R2 {

    tag "Generate mask for low or excess coverage sites (round2)"
    publishDir path: { "${params.outputdir}/01.variants/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}

process mask_cov_R3 {

    tag "Generate mask for low or excess coverage sites (round3)"
    publishDir path: { "${params.outputdir}/01.variants/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}

process mask_cov_R4 {

    tag "Generate mask for low or excess coverage sites (round4)"
    publishDir path: { "${params.outputdir}/01.variants/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}

/*
 * Merge masks per round
 */

process mask_merge_R1 {

    tag "Merge low cov and het mask files (round1)"
    publishDir path: { "${params.outputdir}/01.variants/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov_hets.tsv")

    script:
    """
    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\\t/g' ${cov_bed}
    sed -i 's/ /\\t/g' ${het_bed}

    cat ${cov_bed} ${het_bed} | \
      sort -Vk1 -Vk2 | \
      uniq > ${chromo}_cov_hets.tsv

    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}

process mask_merge_R2 {

    tag "Merge low cov and het mask files (round2)"
    publishDir path: { "${params.outputdir}/01.variants/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov_hets.tsv")

    script:
    """
    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\\t/g' ${cov_bed}
    sed -i 's/ /\\t/g' ${het_bed}

    cat ${cov_bed} ${het_bed} | \
      sort -Vk1 -Vk2 | \
      uniq > ${chromo}_cov_hets.tsv

    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}

process mask_merge_R3 {

    tag "Merge low cov and het mask files (round3)"
    publishDir path: { "${params.outputdir}/01.variants/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov_hets.tsv")

    script:
    """
    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\\t/g' ${cov_bed}
    sed -i 's/ /\\t/g' ${het_bed}

    cat ${cov_bed} ${het_bed} | \
      sort -Vk1 -Vk2 | \
      uniq > ${chromo}_cov_hets.tsv

    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}

process mask_merge_R4 {

    tag "Merge low cov and het mask files (round4)"
    publishDir path: { "${params.outputdir}/01.variants/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov_hets.tsv")

    script:
    """
    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\\t/g' ${cov_bed}
    sed -i 's/ /\\t/g' ${het_bed}

    cat ${cov_bed} ${het_bed} | \
      sort -Vk1 -Vk2 | \
      uniq > ${chromo}_cov_hets.tsv

    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}

/*
 * Consensus per round
 * R1: standard
 * R2–R4: use bcftools +fixref to ensure REF matches FASTA
 */

process call_consensus_R1 {

    tag "Call consensus without masking (round1)"
    publishDir path: { "${params.outputdir}/02.consensus/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${vcf_fn} -o ${individual}_${chromo}_cons.fa
    """
}

process call_consensus_MASK_R1 {

    tag "Call consensus with masking (round1)"
    publishDir path: { "${params.outputdir}/02.consensus/round1/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${vcf_fn} -m ${mask_fn} -o ${individual}_${chromo}_cons.fa
    """
}

/* R2: with fixref */

process call_consensus_R2 {

    tag "Call consensus without masking (round2)"
    publishDir path: { "${params.outputdir}/02.consensus/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    # Index original VCF
    tabix -p vcf ${vcf_fn}

    # Fix REF alleles to match reference FASTA, output gzipped VCF
    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    # Index fixed VCF
    tabix -p vcf ${chromo}_fixed.vcf.gz

    # Build consensus from fixed VCF
    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -o ${individual}_${chromo}_cons.fa
    """
}

process call_consensus_MASK_R2 {

    tag "Call consensus with masking (round2)"
    publishDir path: { "${params.outputdir}/02.consensus/round2/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    # Index original VCF
    tabix -p vcf ${vcf_fn}

    # Fix REF alleles to match reference FASTA, output gzipped VCF
    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    # Index fixed VCF
    tabix -p vcf ${chromo}_fixed.vcf.gz

    # Build masked consensus from fixed VCF
    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -m ${mask_fn} -o ${individual}_${chromo}_cons.fa
    """
}

/* R3: with fixref */

process call_consensus_R3 {

    tag "Call consensus without masking (round3)"
    publishDir path: { "${params.outputdir}/02.consensus/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    tabix -p vcf ${chromo}_fixed.vcf.gz

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -o ${individual}_${chromo}_cons.fa
    """
}

process call_consensus_MASK_R3 {

    tag "Call consensus with masking (round3)"
    publishDir path: { "${params.outputdir}/02.consensus/round3/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    tabix -p vcf ${chromo}_fixed.vcf.gz

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -m ${mask_fn} -o ${individual}_${chromo}_cons.fa
    """
}

/* R4: with fixref */

process call_consensus_R4 {

    tag "Call consensus without masking (round4)"
    publishDir path: { "${params.outputdir}/02.consensus/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    tabix -p vcf ${chromo}_fixed.vcf.gz

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -o ${individual}_${chromo}_cons.fa
    """
}

process call_consensus_MASK_R4 {

    tag "Call consensus with masking (round4)"
    publishDir path: { "${params.outputdir}/02.consensus/round4/${individual}" }, mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    """
    tabix -p vcf ${vcf_fn}

    bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- -f ${reference} -m flip -d

    tabix -p vcf ${chromo}_fixed.vcf.gz

    samtools faidx ${reference} ${chromo} | \
      bcftools consensus ${chromo}_fixed.vcf.gz -m ${mask_fn} -o ${individual}_${chromo}_cons.fa
    """
}

/*
 * Build per-individual references from consensus (rounds 1–3)
 * NOTE: cleaned to A/C/G/T/N only
 */

process build_ref_from_consensus_R1 {

    tag "Build per-individual reference for next round (round1→2)"
    publishDir "${params.outputdir}/02.consensus_refs/round1", mode:'copy'

    input:
    tuple val(individual), val(round), path(cons_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_ref.fa")

    script:
    """
    # Concatenate per-chromosome consensus FASTAs
    cat ${cons_list} > ${individual}_round${round}_ref.tmp.fa

    # Clean reference: convert any non-ACGT character to N
    awk '
      /^>/ { print; next }
      { gsub(/[^ACGTacgt]/, "N"); print }
    ' ${individual}_round${round}_ref.tmp.fa > ${individual}_round${round}_ref.fa

    rm ${individual}_round${round}_ref.tmp.fa
    """
}

process build_ref_from_consensus_R2 {

    tag "Build per-individual reference for next round (round2→3)"
    publishDir "${params.outputdir}/02.consensus_refs/round2", mode:'copy'

    input:
    tuple val(individual), val(round), path(cons_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_ref.fa")

    script:
    """
    cat ${cons_list} > ${individual}_round${round}_ref.tmp.fa

    awk '
      /^>/ { print; next }
      { gsub(/[^ACGTacgt]/, "N"); print }
    ' ${individual}_round${round}_ref.tmp.fa > ${individual}_round${round}_ref.fa

    rm ${individual}_round${round}_ref.tmp.fa
    """
}

process build_ref_from_consensus_R3 {

    tag "Build per-individual reference for next round (round3→4)"
    publishDir "${params.outputdir}/02.consensus_refs/round3", mode:'copy'

    input:
    tuple val(individual), val(round), path(cons_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_ref.fa")

    script:
    """
    cat ${cons_list} > ${individual}_round${round}_ref.tmp.fa

    awk '
      /^>/ { print; next }
      { gsub(/[^ACGTacgt]/, "N"); print }
    ' ${individual}_round${round}_ref.tmp.fa > ${individual}_round${round}_ref.fa

    rm ${individual}_round${round}_ref.tmp.fa
    """
}

/*
 * Missing data stats on final consensus
 */

process calc_missing_data_INDIV {

    tag "Calculate missing data per individual and chromosome"
    publishDir "${params.outputdir}/02.consensus/final", mode:'copy'

    input:
    tuple val(individual), path(cons_fn_list)

    output:
    file("${individual}_missing_data.tsv")

    script:
    """
    02_cons_stats_INDIV.py \
      -c ${cons_fn_list} \
      -i ${individual}
    """
}

process calc_missing_data_SUMMARY {

    tag "Calculate missing data across all individuals"
    publishDir "${params.outputdir}/02.consensus/final", mode:'copy'

    input:
    path(cons_fn_list)

    output:
    file("all_indivs_missing_data.tsv")
    file("all_indivs_missing_data.pdf")

    script:
    """
    unset DISPLAY

    03_cons_stats_SUMMARY.py \
      -c ${cons_fn_list}
    """
}

/*
============================================
~ ~ ~ > *  Pipeline specification  * < ~ ~ ~
============================================
*/

workflow {

    def rounds = params.n_rounds as int

    /*
     * Round 1
     */

    log.info "===================================="
    log.info "      ITERATIVE ROUND 1"
    log.info "===================================="

    // Mapping R1: attach original reference file to each indiv/read tuple
    def map_input_R1_ch = reads_ch.map { row ->
        def indiv = row[0]
        def r1    = row[1]
        def r2    = row[2]
        tuple(indiv, r1, r2, ref_file_obj)
    }

    def bam_R1_ch         = bwa_map_R1(map_input_R1_ch)
    def bam_R1_chromo_ch  = bam_R1_ch.combine(chromos_ch)
    def bam_R1_chromo_ref_ch = bam_R1_chromo_ch.map { row ->
        def indiv  = row[0]
        def bam    = row[1]
        def bai    = row[2]
        def chromo = row[3]
        tuple(indiv, bam, bai, chromo, ref_file_obj)
    }

    // Coverage (only R1)
    if (params.calc_coverage) {
        cov_estimate(bam_R1_chromo_ref_ch)
        cov_estimate.out.groupTuple() | cov_summary_INDIV
        cov_summary_INDIV.out.collect() | cov_summary_ALL
    }

    // Variants R1
    def vars_filt_R1_ch
    if (params.indiv_var_call && params.filt_indels) {
        vars_filt_R1_ch = call_variants_CHROMO_R1(bam_R1_chromo_ref_ch) | remove_indels_R1
    } else if (params.indiv_var_call && !params.filt_indels) {
        vars_filt_R1_ch = call_variants_CHROMO_R1(bam_R1_chromo_ref_ch)
    } else {
        vars_filt_R1_ch = Channel.empty()
    }

    // Masks R1
    def mask_R1_fn_ch = null
    if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
        def mask_het_R1_ch  = vars_filt_R1_ch | mask_hets_R1
        def mask_cov_R1_ch  = mask_cov_R1(bam_R1_chromo_ref_ch)
        def mask_comb_R1_ch = mask_het_R1_ch.combine(mask_cov_R1_ch, by: [0,1])
        mask_R1_fn_ch = mask_merge_R1(mask_comb_R1_ch)
    } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
        mask_R1_fn_ch = vars_filt_R1_ch | mask_hets_R1
    } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
        mask_R1_fn_ch = mask_cov_R1(bam_R1_chromo_ref_ch)
    }

    // Consensus R1
    def consensus_R1_ch = null
    if (params.indiv_var_call && params.call_consensus) {

        if (params.mask_hets || params.mask_cov) {
            def vars_mask_R1_ch = vars_filt_R1_ch.combine(mask_R1_fn_ch, by: [0,1])
            def vars_mask_ref_R1_ch = vars_mask_R1_ch.map { row ->
                def indiv  = row[0]
                def chromo = row[1]
                def vcf    = row[2]
                def mask   = row[3]
                tuple(indiv, chromo, vcf, mask, ref_file_obj)
            }
            consensus_R1_ch = call_consensus_MASK_R1(vars_mask_ref_R1_ch)
        } else {
            def vcf_ref_R1_ch = vars_filt_R1_ch.map { row ->
                def indiv  = row[0]
                def chromo = row[1]
                def vcf    = row[2]
                tuple(indiv, chromo, vcf, ref_file_obj)
            }
            consensus_R1_ch = call_consensus_R1(vcf_ref_R1_ch)
        }
    }

    // Build per-individual references for Round 2
    def indiv_ref_R2_ch = null
    if (rounds >= 2) {
        def grouped_R1_ch = consensus_R1_ch.groupTuple(by: 0)
        def next_ref_input_R1_ch = grouped_R1_ch.map { row ->
            def indiv     = row[0]
            def cons_list = row[2].flatten()
            tuple(indiv, 2, cons_list)
        }
        def built_ref_R1_ch = build_ref_from_consensus_R1(next_ref_input_R1_ch)
        indiv_ref_R2_ch = built_ref_R1_ch.map { row ->
            def indiv = row[0]
            def ref   = row[2]
            tuple(indiv, ref)
        }
    }

    /*
     * Round 2
     */
    def consensus_R2_ch = null
    def indiv_ref_R3_ch = null

    if (rounds >= 2) {
        log.info "===================================="
        log.info "      ITERATIVE ROUND 2"
        log.info "===================================="

        def map_input_R2_ch = reads_ch
            .combine(indiv_ref_R2_ch, by: 0)
            .map { row ->
                def indiv = row[0]
                def r1    = row[1]
                def r2    = row[2]
                def ref   = row[3]
                tuple(indiv, r1, r2, ref)
            }

        def bam_R2_ch = bwa_map_R2(map_input_R2_ch)

        def bam_R2_chromo_ch = bam_R2_ch.combine(chromos_ch)

        def bam_R2_chromo_ref_ch = bam_R2_chromo_ch
            .combine(indiv_ref_R2_ch, by: 0)
            .map { row ->
                def indiv  = row[0]
                def bam    = row[1]
                def bai    = row[2]
                def chromo = row[3]
                def ref    = row[4]
                tuple(indiv, bam, bai, chromo, ref)
            }

        // Variants R2
        def vars_filt_R2_ch
        if (params.indiv_var_call && params.filt_indels) {
            vars_filt_R2_ch = call_variants_CHROMO_R2(bam_R2_chromo_ref_ch) | remove_indels_R2
        } else if (params.indiv_var_call && !params.filt_indels) {
            vars_filt_R2_ch = call_variants_CHROMO_R2(bam_R2_chromo_ref_ch)
        } else {
            vars_filt_R2_ch = Channel.empty()
        }

        // Masks R2
        def mask_R2_fn_ch = null
        if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
            def mask_het_R2_ch  = vars_filt_R2_ch | mask_hets_R2
            def mask_cov_R2_ch  = mask_cov_R2(bam_R2_chromo_ref_ch)
            def mask_comb_R2_ch = mask_het_R2_ch.combine(mask_cov_R2_ch, by: [0,1])
            mask_R2_fn_ch = mask_merge_R2(mask_comb_R2_ch)
        } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
            mask_R2_fn_ch = vars_filt_R2_ch | mask_hets_R2
        } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
            mask_R2_fn_ch = mask_cov_R2(bam_R2_chromo_ref_ch)
        }

        // Consensus R2
        if (params.indiv_var_call && params.call_consensus) {

            if (params.mask_hets || params.mask_cov) {
                def vars_mask_R2_ch = vars_filt_R2_ch.combine(mask_R2_fn_ch, by: [0,1])
                def vars_mask_ref_R2_ch = vars_mask_R2_ch
                    .combine(indiv_ref_R2_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def mask   = row[3]
                        def ref    = row[4]
                        tuple(indiv, chromo, vcf, mask, ref)
                    }
                consensus_R2_ch = call_consensus_MASK_R2(vars_mask_ref_R2_ch)
            } else {
                def vcf_ref_R2_ch = vars_filt_R2_ch
                    .combine(indiv_ref_R2_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def ref    = row[3]
                        tuple(indiv, chromo, vcf, ref)
                    }
                consensus_R2_ch = call_consensus_R2(vcf_ref_R2_ch)
            }
        }

        // Build per-individual references for Round 3
        if (rounds >= 3) {
            def grouped_R2_ch = consensus_R2_ch.groupTuple(by: 0)
            def next_ref_input_R2_ch = grouped_R2_ch.map { row ->
                def indiv     = row[0]
                def cons_list = row[2].flatten()
                tuple(indiv, 3, cons_list)
            }
            def built_ref_R2_ch = build_ref_from_consensus_R2(next_ref_input_R2_ch)
            indiv_ref_R3_ch = built_ref_R2_ch.map { row ->
                def indiv = row[0]
                def ref   = row[2]
                tuple(indiv, ref)
            }
        }
    }

    /*
     * Round 3
     */
    def consensus_R3_ch = null
    def indiv_ref_R4_ch = null

    if (rounds >= 3) {
        log.info "===================================="
        log.info "      ITERATIVE ROUND 3"
        log.info "===================================="

        def map_input_R3_ch = reads_ch
            .combine(indiv_ref_R3_ch, by: 0)
            .map { row ->
                def indiv = row[0]
                def r1    = row[1]
                def r2    = row[2]
                def ref   = row[3]
                tuple(indiv, r1, r2, ref)
            }

        def bam_R3_ch = bwa_map_R3(map_input_R3_ch)

        def bam_R3_chromo_ch = bam_R3_ch.combine(chromos_ch)

        def bam_R3_chromo_ref_ch = bam_R3_chromo_ch
            .combine(indiv_ref_R3_ch, by: 0)
            .map { row ->
                def indiv  = row[0]
                def bam    = row[1]
                def bai    = row[2]
                def chromo = row[3]
                def ref    = row[4]
                tuple(indiv, bam, bai, chromo, ref)
            }

        // Variants R3
        def vars_filt_R3_ch
        if (params.indiv_var_call && params.filt_indels) {
            vars_filt_R3_ch = call_variants_CHROMO_R3(bam_R3_chromo_ref_ch) | remove_indels_R3
        } else if (params.indiv_var_call && !params.filt_indels) {
            vars_filt_R3_ch = call_variants_CHROMO_R3(bam_R3_chromo_ref_ch)
        } else {
            vars_filt_R3_ch = Channel.empty()
        }

        // Masks R3
        def mask_R3_fn_ch = null
        if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
            def mask_het_R3_ch  = vars_filt_R3_ch | mask_hets_R3
            def mask_cov_R3_ch  = mask_cov_R3(bam_R3_chromo_ref_ch)
            def mask_comb_R3_ch = mask_het_R3_ch.combine(mask_cov_R3_ch, by: [0,1])
            mask_R3_fn_ch = mask_merge_R3(mask_comb_R3_ch)
        } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
            mask_R3_fn_ch = vars_filt_R3_ch | mask_hets_R3
        } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
            mask_R3_fn_ch = mask_cov_R3(bam_R3_chromo_ref_ch)
        }

        // Consensus R3
        if (params.indiv_var_call && params.call_consensus) {

            if (params.mask_hets || params.mask_cov) {
                def vars_mask_R3_ch = vars_filt_R3_ch.combine(mask_R3_fn_ch, by: [0,1])
                def vars_mask_ref_R3_ch = vars_mask_R3_ch
                    .combine(indiv_ref_R3_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def mask   = row[3]
                        def ref    = row[4]
                        tuple(indiv, chromo, vcf, mask, ref)
                    }
                consensus_R3_ch = call_consensus_MASK_R3(vars_mask_ref_R3_ch)
            } else {
                def vcf_ref_R3_ch = vars_filt_R3_ch
                    .combine(indiv_ref_R3_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def ref    = row[3]
                        tuple(indiv, chromo, vcf, ref)
                    }
                consensus_R3_ch = call_consensus_R3(vcf_ref_R3_ch)
            }
        }

        // Build per-individual references for Round 4
        if (rounds >= 4) {
            def grouped_R3_ch = consensus_R3_ch.groupTuple(by: 0)
            def next_ref_input_R3_ch = grouped_R3_ch.map { row ->
                def indiv     = row[0]
                def cons_list = row[2].flatten()
                tuple(indiv, 4, cons_list)
            }
            def built_ref_R3_ch = build_ref_from_consensus_R3(next_ref_input_R3_ch)
            indiv_ref_R4_ch = built_ref_R3_ch.map { row ->
                def indiv = row[0]
                def ref   = row[2]
                tuple(indiv, ref)
            }
        }
    }

    /*
     * Round 4 (final, optional)
     */
    def consensus_R4_ch = null

    if (rounds >= 4) {
        log.info "===================================="
        log.info "      ITERATIVE ROUND 4 (final)"
        log.info "===================================="

        def map_input_R4_ch = reads_ch
            .combine(indiv_ref_R4_ch, by: 0)
            .map { row ->
                def indiv = row[0]
                def r1    = row[1]
                def r2    = row[2]
                def ref   = row[3]
                tuple(indiv, r1, r2, ref)
            }

        def bam_R4_ch = bwa_map_R4(map_input_R4_ch)

        def bam_R4_chromo_ch = bam_R4_ch.combine(chromos_ch)

        def bam_R4_chromo_ref_ch = bam_R4_chromo_ch
            .combine(indiv_ref_R4_ch, by: 0)
            .map { row ->
                def indiv  = row[0]
                def bam    = row[1]
                def bai    = row[2]
                def chromo = row[3]
                def ref    = row[4]
                tuple(indiv, bam, bai, chromo, ref)
            }

        // Variants R4
        def vars_filt_R4_ch
        if (params.indiv_var_call && params.filt_indels) {
            vars_filt_R4_ch = call_variants_CHROMO_R4(bam_R4_chromo_ref_ch) | remove_indels_R4
        } else if (params.indiv_var_call && !params.filt_indels) {
            vars_filt_R4_ch = call_variants_CHROMO_R4(bam_R4_chromo_ref_ch)
        } else {
            vars_filt_R4_ch = Channel.empty()
        }

        // Masks R4
        def mask_R4_fn_ch = null
        if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
            def mask_het_R4_ch  = vars_filt_R4_ch | mask_hets_R4
            def mask_cov_R4_ch  = mask_cov_R4(bam_R4_chromo_ref_ch)
            def mask_comb_R4_ch = mask_het_R4_ch.combine(mask_cov_R4_ch, by: [0,1])
            mask_R4_fn_ch = mask_merge_R4(mask_comb_R4_ch)
        } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
            mask_R4_fn_ch = vars_filt_R4_ch | mask_hets_R4
        } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
            mask_R4_fn_ch = mask_cov_R4(bam_R4_chromo_ref_ch)
        }

        // Consensus R4 (final)
        if (params.indiv_var_call && params.call_consensus) {

            if (params.mask_hets || params.mask_cov) {
                def vars_mask_R4_ch = vars_filt_R4_ch.combine(mask_R4_fn_ch, by: [0,1])
                def vars_mask_ref_R4_ch = vars_mask_R4_ch
                    .combine(indiv_ref_R4_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def mask   = row[3]
                        def ref    = row[4]
                        tuple(indiv, chromo, vcf, mask, ref)
                    }
                consensus_R4_ch = call_consensus_MASK_R4(vars_mask_ref_R4_ch)
            } else {
                def vcf_ref_R4_ch = vars_filt_R4_ch
                    .combine(indiv_ref_R4_ch, by: 0)
                    .map { row ->
                        def indiv  = row[0]
                        def chromo = row[1]
                        def vcf    = row[2]
                        def ref    = row[3]
                        tuple(indiv, chromo, vcf, ref)
                    }
                consensus_R4_ch = call_consensus_R4(vcf_ref_R4_ch)
            }
        }
    }

    /*
     * Missing data on final consensus
     */
    if (params.calc_missing_data && params.indiv_var_call && params.call_consensus) {

        def final_consensus_ch
        if (rounds == 1) {
            final_consensus_ch = consensus_R1_ch
        } else if (rounds == 2) {
            final_consensus_ch = consensus_R2_ch
        } else if (rounds == 3) {
            final_consensus_ch = consensus_R3_ch
        } else {
            final_consensus_ch = consensus_R4_ch
        }

        def grouped_final_ch = final_consensus_ch.groupTuple(by: 0)
        def md_input_ch = grouped_final_ch.map { row ->
            def indiv   = row[0]
            def cons_fn = row[2].flatten()
            [indiv, cons_fn]
        }

        calc_missing_data_INDIV(md_input_ch)
        calc_missing_data_INDIV.out.collect() | calc_missing_data_SUMMARY
    }
}
