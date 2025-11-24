#!/usr/bin/env nextflow

log.info """\

        ================================================
        NF-varmit - Variant calling for whole-genome resequence data, iteratively
        https://github.com/MozesBlom/nf-varmit   
        Author: Mozes P.K. Blom ; Jon J. Hoffman
        ================================================
        |indivfile    : ${params.indivs_file}
        |chromos      : ${params.chromos_file}
        |reference    : ${params.ref_file}
        |readsdir     : ${params.readsdir}
        |outputdir    : ${params.outputdir}
        |
        |Sex chromos  : ${params.sex_chromos}
        |
        |Rounds of mapping/consensus : ${params.n_rounds}
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
        ================================================
        """
        .stripIndent()

/*
====================================================
~ ~ ~ > *  Input channels and file check  * < ~ ~ ~ 
====================================================
*/

/*
 * Create two lists:
 * - All individuals to include (expects that the read files are called: indiv_R1.fastq.gz / indiv_R2.fastq.gz by default)
 * - All chromosomes to include (must match exactly with reference genome fa ids)
 */

indivs      = file(params.indivs_file).readLines()
indivs_ch   = Channel.fromList(indivs)

chromos     = file(params.chromos_file).readLines()
chromos_ch  = Channel.fromList(chromos)

/*
 * Reads channel:
 * - Paired-end reads per individual:
 *   (individual, R1, R2)
 */
reads_ch = indivs_ch.map { indiv ->
    def r1 = file("${params.readsdir}/${indiv}${params.reads_suffix1}")
    def r2 = file("${params.readsdir}/${indiv}${params.reads_suffix2}")
    if( !r1.exists() ) error "Missing R1 reads for ${indiv}: ${r1}"
    if( !r2.exists() ) error "Missing R2 reads for ${indiv}: ${r2}"
    tuple(indiv, r1, r2)
}

/*
 * Set data channels:
 * - Reference genome for round 1
 */
ref_ch = Channel.fromPath(params.ref_file)
                .ifEmpty { error "No reference sequence found from: ${params.ref_file}" }

/*
 * Default number of iterative rounds (if not set in config/CLI)
 */
params.n_rounds = params.n_rounds ?: 1


/*
===============================
~ ~ ~ > *  Processes  * < ~ ~ ~ 
===============================
*/

/*
 * NOTE: index_bam process is kept for backwards compatibility,
 * but is not used in the iterative mapping workflow below.
 */

process index_bam  {

/*
 * If need be index each of the bam files before proceeding.
 */

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


process bwa_map {

/*
 * Map reads with bwa-mem to a given reference and produce a sorted, indexed BAM.
 */

    tag "BWA-MEM mapping"
    publishDir "${params.outputdir}/00.bam/${individual}", mode: 'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(read1), path(read2), path(reference)

    output:
    tuple val(individual), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    # Index reference if needed
    samtools faidx ${reference}
    bwa index ${reference}

    # Map and sort
    bwa mem -t ${task.cpus} ${reference} ${read1} ${read2} | \
      samtools sort -@ ${task.cpus} -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}


process cov_estimate {

/*
 * Estimate the mean coverage for each individual and each of the focal chromosomes/scaffolds.
 * NOTE: I tried to lump all chromos per individual into one job but it didn't work.
 * so for now it will result in a larger number of short jobs (which is not optimal for a cluster)
 */

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

/*
 * Accurately summarise the coverage across all target chromos per individual
 */

    tag "Coverage summary per individual"
    publishDir "${params.outputdir}/00.coverage/", mode:'copy'

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

/*
 * Accurately summarise the coverage across all target chromos per individual
 */

    tag "Coverage summary for all indivs"
    publishDir "${params.outputdir}/00.coverage/", mode:'copy'

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


process call_variants_CHROMO {

/*
 * Call variants for each individual by chromosome and an initial round of filtering takes place.
 *
 * We use vcfallelicprimitives to deconstruct mnps to snps.
 * My original script has the -g flag after vcf ap. Double check if this is still supported, since it's not listed in the -h list
 */

    tag "Variant calling per chromosome"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'
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


process remove_indels {

/*
 * Optional: Remove indel variation from vcf file
 */

    tag "Remove indels"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

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


process mask_hets {

/*
 * Optional: Create mask file for heterozygous sites
 */

    tag "Generate mask for het sites"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

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


process mask_cov {

/*
 * Optional: Create mask file for low coverage or excess coverage sites
 */

    tag "Generate mask for low or excess coverage sites"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), file(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}


process mask_merge {

/*
 * Optional: Merge two mask files
 */

    tag "Merge low cov and het mask files"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

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

    # Merge two mask files and sort by position. Remove duplicate entries
    cat ${cov_bed} ${het_bed} | \
      sort -Vk1 -Vk2 | \
      uniq > ${chromo}_cov_hets.tsv

    # Add a header for each column to make it a proper bed-like tsv file
    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}


process call_consensus {

/*
 * Call a consensus sequence for each individual by chromosome WITHOUT MASKING
 * (reference is provided per round and per individual/chromosome)
 */

    tag "Call consensus without masking"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

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


process call_consensus_MASK {

/*
 * Call a consensus sequence for each individual by chromosome WITH MASKING
 */

    tag "Call consensus with masking"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

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


process build_ref_from_consensus {

/*
 * Concatenate per-chromosome consensus sequences into a per-individual
 * reference fasta for the next mapping round.
 */

    tag "Build per-individual reference for next round"
    publishDir "${params.outputdir}/02.consensus_refs/${individual}", mode:'copy'

    input:
    tuple val(individual), val(round), path(cons_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_ref.fa")

    script:
    """
    cat ${cons_list} > ${individual}_round${round}_ref.fa
    """
}


process calc_missing_data_INDIV {

/*
 * Optional: Estimate the amount of missing data for each individual and each chromosome
 * NOTE: Within the framework of the current pipeline this can be optimised by using the mask files.
 * However, I already had existing scripts borrowed from nf-phylo to make this work.
 */

    tag "Calculate missing data per individual and chromosome"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

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

/*
 * Optional: Summarise the amount of missing data across all individuals
 * NOTE: Within the framework of the current pipeline this can be optimised by using the mask files.
 * However, I already had existing scripts borrowed from nf-phylo to make this work.
 */

    tag "Calculate missing data across all individual"
    publishDir "${params.outputdir}/02.consensus/", mode:'copy'

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
=========================================================
~ ~ ~ > *  Helper closures for iterative workflow * < ~ ~ ~ 
=========================================================
*/

/*
 * Extract (indiv, chromo, ref) from (indiv, bam, bai, chromo, ref) tuples
 * for combining with variant/consensus channels.
 */
def ref_for_round = { ch ->
    ch.map { row ->
        def indiv  = row[0]
        def chromo = row[3]
        def ref    = row[4]
        tuple(indiv, chromo, ref)
    }
}


/*
============================================
~ ~ ~ > *  Pipeline specification  * < ~ ~ ~ 
============================================
*/

workflow {

    /*
     * Per-individual reference channel for round 1:
     * every individual starts from the same params.ref_file.
     * indiv_ref_ch: (indiv, ref)
     */
    def indiv_ref_ch = indivs_ch.combine(ref_ch)

    // Keep track of last round’s consensus (for stats) and mapping channel
    def consensus_last       = null
    def bam_chromo_ref_last  = null

    /*
     * Main iterative loop:
     *  - map reads to per-individual reference
     *  - call variants + masks
     *  - build consensus
     *  - build new reference for next round
     */
    (1..params.n_rounds).each { round ->

        log.info "===================================="
        log.info "      ITERATIVE ROUND ${round}"
        log.info "===================================="

        //
        // 1) Mapping with bwa-mem: (indiv, R1, R2) + (indiv, ref) → (indiv, bam, bai)
        //
        def map_input_ch = reads_ch
            .combine(indiv_ref_ch, by: 0)
            .map { left, right ->
                /*
                 * left  = [indiv, r1, r2]
                 * right = [indiv, ref]
                 */
                def indiv = left[0]
                def r1    = left[1]
                def r2    = left[2]
                def ref   = right[1]
                tuple(indiv, r1, r2, ref)
            }

        def bam_round_ch = bwa_map(map_input_ch)   // (indiv, bam, bai)

        //
        // 2) Build (indiv, bam, bai, chromo, ref) for this round
        //
        def bam_chromo_ch = bam_round_ch
            .combine(chromos_ch)                  // (indiv, bam, bai, chromo)

        def bam_chromo_ref_ch = bam_chromo_ch
            .combine(indiv_ref_ch, by: 0)         // join on individual
            .map { left, right ->
                /*
                 * left  = [indiv, bam, bai, chromo]
                 * right = [indiv, ref]
                 */
                def indiv  = left[0]
                def bam    = left[1]
                def bai    = left[2]
                def chromo = left[3]
                def ref    = right[1]
                tuple(indiv, bam, bai, chromo, ref)
            }

        bam_chromo_ref_last = bam_chromo_ref_ch   // remember last round’s mapping

        //
        // 3) Coverage estimation (optional; per round)
        //
        if (params.calc_coverage) {
            cov_estimate(bam_chromo_ref_ch)
            cov_estimate.out.groupTuple() | cov_summary_INDIV
            cov_summary_INDIV.out.collect() | cov_summary_ALL
        }

        //
        // 4) Variant calling + filtering
        //
        def vars_filt_fn_ch
        if (params.indiv_var_call && params.filt_indels) {
            vars_filt_fn_ch = call_variants_CHROMO(bam_chromo_ref_ch) | remove_indels
        } else if (params.indiv_var_call && !params.filt_indels) {
            vars_filt_fn_ch = call_variants_CHROMO(bam_chromo_ref_ch)
        } else {
            vars_filt_fn_ch = Channel.empty()
        }

        //
        // 5) Masking (hets / cov)
        //
        def mask_fn_ch = null
        if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
            def mask_het_fn_ch   = vars_filt_fn_ch | mask_hets
            def mask_cov_fn_ch   = mask_cov(bam_chromo_ref_ch)
            def mask_combined_ch = mask_het_fn_ch.combine(mask_cov_fn_ch, by: [0,1])
            mask_fn_ch = mask_merge(mask_combined_ch)
        } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
            mask_fn_ch = vars_filt_fn_ch | mask_hets
        } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
            mask_fn_ch = mask_cov(bam_chromo_ref_ch)
        }

        //
        // 6) Consensus for this round, using round-specific reference
        //
        def consensus_fn = null

        if (params.indiv_var_call && params.call_consensus) {

            // (indiv, chromo, ref) for this round
            def ref_for_this_round = ref_for_round(bam_chromo_ref_ch)

            if (params.mask_hets || params.mask_cov) {
                // vcf + mask + ref → call_consensus_MASK
                def vars_mask_ch = vars_filt_fn_ch.combine(mask_fn_ch, by: [0,1])   // (indiv, chromo, vcf, mask)
                def vars_mask_ref_ch = vars_mask_ch.combine(ref_for_this_round, by: [0,1]) // (indiv, chromo, vcf, mask, ref)

                consensus_fn = call_consensus_MASK(vars_mask_ref_ch)
            } else {
                // vcf + ref → call_consensus
                def vcf_with_ref_ch = vars_filt_fn_ch.combine(ref_for_this_round, by: [0,1]) // (indiv, chromo, vcf, ref)
                consensus_fn = call_consensus(vcf_with_ref_ch)
            }
        }

        consensus_last = consensus_fn  // store final consensus of this round

        //
        // 7) Build per-individual reference for the next round
        //
        if (round < params.n_rounds && params.indiv_var_call && params.call_consensus) {

            def new_ref_ch = consensus_fn
                .groupTuple(by: 0)     // (indiv, [chromo_list], [cons_list])
                .map { row ->
                    def indiv     = row[0]
                    def cons_list = row[2].flatten()
                    tuple(indiv, round + 1, cons_list)
                } \
                | build_ref_from_consensus

            // For the next round we only need (indiv, ref)
            indiv_ref_ch = new_ref_ch.map { row ->
                def indiv = row[0]
                def ref   = row[2]
                tuple(indiv, ref)
            }
        }

    } // end rounds loop

    //
    // 8) After last round: missing data on final consensus (optional)
    //
    if (params.calc_missing_data && consensus_last) {
        consensus_last
            .groupTuple(by: 0)
            .map { it ->
                def indiv   = it[0]
                def cons_fn = it[2].flatten()
                [indiv, cons_fn]
            } | calc_missing_data_INDIV

        calc_missing_data_INDIV.out.collect() | calc_missing_data_SUMMARY
    }
}
