#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-rnasnv: RNA-Seq variant calling nextflow pipeline
 *
 * Copyright (C) 2019 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


import nextflow.splitter.CsvSplitter
import nextflow.config.ConfigParser

nextflow_config = file("${baseDir}/nextflow.config").text
parsed_config = new ConfigParser().setIgnoreIncludes(true).parse(nextflow_config)
defaults = parsed_config.params

check_params()


/*
 * Log workflow options
 */

// Nextflow exectution options
log.info("Config profile: ${workflow.profile}")
log.info("Workflow revision: ${workflow.revision}")
log.info("Workflow work-dir: ${workflow.workDir}")

// Input options
log.info("Input csv: ${params.csv}")
log.info("Paired-end readgroups: ${String.valueOf(params.paired_end)}")

// Reference genome options
log.info("Reference genome: ${params.genome}")
log.info("Reference FASTA file: ${params.fasta}")
log.info("Reference GTF file: ${params.gtf}")

// Sequencing adapter options
log.info("Adapters to trim: ${params.adapters}")
log.info("R1 adapter: ${params.r1_adapter}")
log.info("R2 adapter: ${params.r2_adapter}")

// Output options
log.info("Reference directory: ${params.refdir}")
log.info("Output directory: ${params.outdir}")


/*
 * Log advanced options
 */

// STAR genome generate options
log.info("STAR size of the bins for genome storage: ${params.star_genome_chr_bin_n_bits}")
log.info("STAR length of the SA pre-indexing string: ${params.star_genome_sa_index_n_bases}")
log.info("STAR SJDB overhang: ${params.star_sjdb_overhang}")

// Cutadapt options
log.info("Cutadapt base quality cutoff: ${params.cutadapt_base_qual_cutoff}")
log.info("Cutadapt min read length: ${params.cutadapt_min_read_length}")

// BCFtools mpileup options
log.info("BCFtools mpileup num regions per process: ${defaults.mpileup_num_regions}")
log.info("BCFtools mpileup max depth: ${defaults.mpileup_max_depth}")
log.info("BCFtools mpileup min base quality: ${defaults.mpileup_min_bq}")

// Ensembl VEP options
log.info("Ensembl VEP cache type: ${params.vep_cache_type}")
log.info("Ensembl VEP indexed cache file: ${params.vep_indexed_cache_file}")
log.info("Ensembl VEP species name: ${params.vep_species_name}")
log.info("Ensembl VEP assembly name: ${params.vep_assembly_name}")

// MultiQC options
log.info("MultiQC config: ${params.multiqc_config}")

// Reports options
log.info("Execution report: ${params.execution_report}")
log.info("Trace report: ${params.trace_report}")
log.info("Timeline report: ${params.timeline_report}")
log.info("Flowchart: ${params.flowchart}")

// AWS Batch options
log.info("AWS Batch JobQueue: ${params.aws_queue}")
log.info("AWS Region: ${params.aws_region}")


/*
 * Validate input readgroups
 */

validate_input_csv()


/*
 * File placeholders
 */

csv_file = file(params.csv)
fasta_file = file(params.fasta)
gtf_file = file(params.gtf)

vep_indexed_cache_file = file(params.vep_indexed_cache_file)

multiqc_cfg = file(params.multiqc_config)




/*
 * PREPROCESSING - Gunzip reference FASTA file
 */
process gunzip_fasta {
    storeDir params.refdir

    input:
    file fasta from fasta_file

    output:
    file "${fasta.getBaseName()}" into gunzipped_fasta

    when:
    fasta.getExtension() == "gz"

    """
    gzip -dc "${fasta}" > "${fasta.getBaseName()}"
    """
}


/*
 * PREPROCESSING - Gunzip reference GTF file
 */
process gunzip_gtf {
    storeDir params.refdir

    input:
    file gtf from gtf_file

    output:
    file "${gtf.getBaseName()}" into gunzipped_gtf

    when:
    gtf.getExtension() == "gz"

    """
    gzip -dc "${gtf}" > "${gtf.getBaseName()}"
    """
}


/*
 * PREPROCESSING - Create a value channel for the reference FASTA file
 */
gunzipped_fasta
    .ifEmpty { fasta_file }
    .set { ref_fasta }


/*
 * PREPROCESSING - Create a value channel for the reference GTF file
 */
gunzipped_gtf
    .ifEmpty { gtf_file }
    .set { ref_gtf }


/*
 * PREPROCESSING - STAR index
 */
process star_index {
    storeDir "${params.refdir}/STAR"

    label 'star'

    input:
    file ref_fasta
    file ref_gtf

    output:
    file "${params.genome}" into star_index

    """
    mkdir "${params.genome}"
    STAR \\
        --runThreadN "${task.cpus}" \\
        --runMode genomeGenerate \\
        --genomeDir "${params.genome}" \\
        --genomeFastaFiles "${ref_fasta}" \\
        --genomeChrBinNbits "${params.star_genome_chr_bin_n_bits}" \\
        --genomeSAindexNbases "${params.star_genome_sa_index_n_bases}" \\
        --sjdbGTFfile "${ref_gtf}" \\
        --sjdbOverhang "${params.star_sjdb_overhang}"
    """
}


/*
 * PREPROCESSING - SAMtools faidx
 */
process samtools_faidx {
    storeDir params.refdir

    label 'samtools'

    input:
    file ref_fasta

    output:
    file "${ref_fasta}.fai" into ref_faidx

    """
    samtools faidx "${ref_fasta}"
    """
}


/*
 * PREPROCESSING - Unpack the indexed VEP cache files
 */
process unpack_cache {
    storeDir "${params.refdir}/VEP"

    input:
    file vep_indexed_cache_file

    output:
    file "cache/*" into indexed_vep_cache

    """
    mkdir cache
    tar -xf "${vep_indexed_cache_file}" -C cache
    """
}


/*
 * PREPROCESSING - Parse in the readgroups and create the inputs for FastQC and Cutadapt
 */
(rgids, r1_inputs, r2_inputs) = Channel
    .fromPath( params.csv )
    .splitCsv( header:true )
    .ifEmpty { exit 1, "No readgroups found post-validation. Exiting." }
    .separate(3) { row ->
        def rgid = row.readgroup ? [row.sample, row.readgroup].join(params.rgid_sep) : row.sample

        // files to stage
        def fq1 = params.paired_end ? file(row.fastq1) : file(row.fastq)
        def fq2 = params.paired_end ? file(row.fastq2) : 'NO_FILE'

        // filenames to stage with
        def fq1_fn = params.paired_end ? "${rgid}.1.${get_fastq_extn(fq1)}" : "${rgid}.${get_fastq_extn(fq1)}"
        def fq2_fn = params.paired_end ? "${rgid}.2.${get_fastq_extn(fq2)}" : 'NO_FILE'

        tuple(tuple(row.sample, rgid), tuple(fq1_fn, fq1), tuple(fq2_fn, fq2))
    }

rgids.into { samples; fastqc_raw_rgids; cutadapt_rgids }
r1_inputs.into { fastqc_raw_r1_inputs; cutadapt_r1_inputs }
r2_inputs.into { fastqc_raw_r2_inputs; cutadapt_r2_inputs }


/*
 * STEP 1 - Run FastQC on the raw reads
 */
process fastqc_raw {
    tag { rgid }

    label 'fastqc'

    publishDir "${params.outdir}/FastQC/${rgid}/raw", mode: 'copy'

    input:
    set sample, rgid from fastqc_raw_rgids
    set fq1, file("${fq1}") from fastqc_raw_r1_inputs
    set fq2, file("${fq2}") from fastqc_raw_r2_inputs

    output:
    file "*_fastqc.{zip,html}" into fastqc_raw_results

    script:
    if (params.paired_end) {

        """
        fastqc -q "${fq1}" "${fq2}"
        """

    } else {

        """
        fastqc -q "${fq1}"
        """
    }
}


/*
 * STEP 2 - Trim adapters using Cutadapt
 */
process cutadapt {
    tag { rgid }

    label 'cutadapt'

    publishDir "${params.outdir}/Cutadapt/${rgid}", mode: 'copy'

    input:
    set sample, rgid from cutadapt_rgids
    set fq1, file("${rgid}/${fq1}") from cutadapt_r1_inputs
    set fq2, file("${rgid}/${fq2}") from cutadapt_r2_inputs

    output:
    set rgid, file("*.fastq.gz") into fastqc_trimmed_inputs, trimmed_readgroups
    file "*.log" into cutadapt_logs

    script:
    def r1_adapter = params.r1_adapter != 'NO_R1_ADAPTER' ? "-a ${params.r1_adapter}" : ""
    def r2_adapter = params.r2_adapter != 'NO_R2_ADAPTER' ? "-A ${params.r2_adapter}" : ""

    if (params.paired_end) {

        """
        cutadapt \\
            "${r1_adapter}" \\
            "${r2_adapter}" \\
            -q "${params.cutadapt_base_qual_cutoff}" \\
            -m "${params.cutadapt_min_read_length}" \\
            --trim-n \\
            -o "${rgid}.1.fastq.gz" \\
            -p "${rgid}.2.fastq.gz" \\
            "${rgid}/${fq1}" \\
            "${rgid}/${fq2}" \\
            > "${rgid}.log"
        """

    } else {

        """
        cutadapt \\
            "${r1_adapter}" \\
            -q "${params.cutadapt_base_qual_cutoff}" \\
            -m "${params.cutadapt_min_read_length}" \\
            --trim-n \\
            -o "${rgid}.fastq.gz" \\
            "${rgid}/${fq1}" \\
            > "${rgid}.log"
        """
    }
}


/*
 * STEP 3 - Run FastQC on the trimmed reads
 */
process fastqc_trimmed {
    tag { rgid }

    label 'fastqc'

    publishDir "${params.outdir}/FastQC/${rgid}/trimmed", mode: 'copy'

    input:
    set rgid, file(fastqs) from fastqc_trimmed_inputs

    output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed_results

    script:
    if (params.paired_end) {

        """
        fastqc -q "${rgid}.1.fastq.gz" "${rgid}.2.fastq.gz"
        """

    } else {

        """
        fastqc -q "${rgid}.fastq.gz"
        """
    }
}


/*
 * STEP 4 - Group the trimmed readgroups by sample
 */
samples
    .groupTuple()
    .map { sample, rgids ->
        tuple( groupKey(sample, rgids.size()), rgids )
    }
    .transpose()
    .map { sample, readgroup ->
        tuple(readgroup, sample)
    }
    .join( trimmed_readgroups )
    .groupTuple( by:1 )
    .map { rgids, sample, fastqs ->
        tuple(sample.toString(), rgids, fastqs.flatten())
    }
    .set { star_inputs }


/*
 * STEP 5 - Align the trimmed reads using STAR
 */
process star {
    tag { sample }

    label 'star'

    publishDir "${params.outdir}/STAR/${sample}", mode: 'copy'

    input:
    set sample, rgids, file(fastqs) from star_inputs
    file index from star_index.collect()

    output:
    set sample, file("*.Aligned.out.bam") into star_alignments
    set sample, file("*.Aligned.sortedByCoord.out.bam") into star_csorted_bam_files, star_sorted_alignments
    set sample, file("*.out.bg") into star_signal_output
    set sample, file("*.SJ.out.tab") into star_splice_junctions
    file "*.out" into star_logs

    script:
    def rgs = rgids.collect { /"ID:${it}" "SM:${sample}"/ }.join(" , ")

    if (params.paired_end) {

        def fq1 = rgids.collect { "${it}.1.fastq.gz" }.join(',')
        def fq2 = rgids.collect { "${it}.2.fastq.gz" }.join(',')

        """
        STAR \\
            --genomeDir "${index}" \\
            --readFilesCommand zcat \\
            --readFilesIn "${fq1}" "${fq2}" \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outSAMattributes All \\
            --outSAMattrRGline ${rgs} \\
            --outSAMtype BAM Unsorted SortedByCoordinate \\
            --outSAMunmapped Within KeepPairs \\
            --outWigType bedGraph \\
            --outFileNamePrefix "${sample}." \\
            --outFilterMultimapNmax 1
        """

    } else {

        def fq = rgids.collect { "${it}.fastq.gz" }.join(',')

        """
        STAR \\
            --genomeDir "${index}" \\
            --readFilesCommand zcat \\
            --readFilesIn "${fq}" \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outSAMattributes All \\
            --outSAMattrRGline ${rgs} \\
            --outSAMtype BAM Unsorted SortedByCoordinate \\
            --outSAMunmapped Within KeepPairs \\
            --outWigType bedGraph \\
            --outFileNamePrefix "${sample}." \\
            --outFilterMultimapNmax 1
        """
    }
}


/*
 * STEP 5.1 - Index the coordinate-sorted alignments
 */
process samtools_index {
    tag { sample }

    label 'samtools'

    publishDir "${params.outdir}/STAR/${sample}", mode: 'copy'

    input:
    set sample, file(bam) from star_csorted_bam_files

    output:
    file "*.bai"

    """
    samtools index "${bam}"
    """
}


/*
 * STEP 6 - Mark duplicates using Picard
 */
process mark_duplicates {
    tag { sample }

    label 'picard'

    publishDir "${params.outdir}/MarkDuplicates/${sample}", mode: 'copy', saveAs: { fn ->
        fn.endsWith(".bai") ? "${sample}.bam.bai" : "${fn}"
    }

    input:
    set sample, file(bam) from star_sorted_alignments

    output:
    set sample, file("*.bam"), file("*.bai") into duplicate_marked_alignments, mpileup_alignments
    file "*.metrics.txt" into mark_duplicates_metrics

    """
    picard \\
        -Xmx${task.memory.toGiga()}g \\
        -XX:+UseSerialGC \\
    MarkDuplicates \\
        INPUT="${bam}" \\
        OUTPUT="${sample}.bam" \\
        METRICS_FILE="${sample}.metrics.txt" \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}


/*
 * STEP 7 - Call variants using Strelka
 */
process strelka {
    tag { sample }

    label 'strelka'

    publishDir "${params.outdir}/Strelka/${sample}", mode: 'copy'

    input:
    set sample, file(bam), file("${bam}.bai") from duplicate_marked_alignments
    file ref_fasta
    file ref_faidx

    output:
    set sample, file("${sample}.vcf.gz{,.tbi}") into strelka_variants

    script:
    def variants = "StrelkaGermlineWorkflow/results/variants/variants.vcf.gz"

    """
    configureStrelkaGermlineWorkflow.py \\
        --bam "${bam}" \\
        --rna \\
        --referenceFasta "${ref_fasta}"
    ./StrelkaGermlineWorkflow/runWorkflow.py \\
        -m local \\
        -j ${task.cpus}
    mv -v "${variants}" "${sample}.vcf.gz"
    mv -v "${variants}.tbi" "${sample}.vcf.gz.tbi"
    """
}


/*
 * STEP 8 - Subset sites passing all filters
 */
process subset_pass_variants {
    tag { sample }

    label 'bcftools'

    publishDir "${params.outdir}/Strelka/${sample}", mode: 'copy'

    input:
    set sample, file(indexed_vcf) from strelka_variants

    output:
    set sample, file("${sample}.pass.vcf.gz{,.tbi}") into pass_variants

    script:
    def (vcf, tbi) = indexed_vcf

    """
    bcftools view \\
        --no-version \\
        -Oz \\
        -o "${sample}.pass.vcf.gz" \\
        -f PASS \\
        "${vcf}"
    bcftools index \\
        -t \\
        "${sample}.pass.vcf.gz"
    """
}


/*
 * STEP 9 - Convert the PASS variants to BED format
 */
process convert2bed {
    tag { sample }

    label 'bedops'

    input:
    set sample, file(indexed_vcf) from pass_variants

    output:
    set sample, file("${sample}.bed.gz") into bed_regions

    script:
    def (vcf, tbi) = indexed_vcf

    """
    zcat "${vcf}" |
        convert2bed -i vcf -d - |
        cut -f -3 |
        gzip > "${sample}.bed.gz"
    """
}


/*
 * STEP 10 - Shuffle and split the BED files into chunks of size `params.mpileup_num_regions`
 */
process split_regions {
    tag { sample }

    label 'coreutils'

    input:
    set sample, file(bed) from bed_regions

    output:
    set sample, file("*.bed") into regions_files

    shell:
    '''
    zcat "!{bed}" | shuf | split \\
        -a "!{params.mpileup_suffix_length}" \\
        -d \\
        -l "!{params.mpileup_num_regions}" \\
        --filter='LC_ALL=C sort -k1,1V -k2,2n -k3,3n > ${FILE}.bed' \\
        - \\
        "!{sample}."
    '''
}


/*
 * STEP 11 - Pileup each regions file, re-call variants and apply soft-filters
 */
process mpileup {
    tag { jobname }

    label 'bcftools'

    input:
    set sample, file(bed), file(bam), file("*") from regions_files
        .map { sample, regions -> tuple( groupKey(sample, [regions].flatten().size()), regions ) }
        .transpose()
        .combine(mpileup_alignments, by: 0)
    file ref_fasta
    file ref_faidx

    output:
    set sample, file("${jobname}.vcf.gz{,.tbi}") into mpileup_region_variants

    script:
    jobname = bed.getBaseName()

    def info_fields = ['AD', 'ADF', 'ADR'].collect { "INFO/$it" }
    def format_fields = ['SP', 'DP', 'AD', 'ADF', 'ADR'].collect { "FORMAT/$it" }

    def mpileup_exclude_filters = params.mpileup_exclude_filters
        .collect { name, expr -> "bcftools filter --no-version -Ou -m+ -s '${name}' -e '${expr}' - |" }
        .join(' ')
    def mpileup_include_filters = params.mpileup_include_filters
        .collect { name, expr -> "bcftools filter --no-version -Ou -m+ -s '${name}' -i '${expr}' - |" }
        .join(' ')

    """
    bcftools mpileup \\
        --no-version \\
        -Ou \\
        -d "${params.mpileup_max_depth}" \\
        -f "${ref_fasta}" \\
        -Q "${params.mpileup_min_bq}" \\
        -R "${bed}" \\
        -a "${info_fields.join(',')},${format_fields.join(',')}" \\
        "${bam}" |
    bcftools norm \\
        --no-version \\
        -Ou \\
        -f "${ref_fasta}" \\
        -m +any \\
        - |
    bcftools call \\
        --no-version \\
        -Ou \\
        -m \\
        -v \\
        -f GQ,GP \\
        - |
    bcftools norm \\
        --no-version \\
        -Ou \\
        -f "${ref_fasta}" \\
        -m -any \\
        - |
    ${mpileup_exclude_filters} \\
    ${mpileup_include_filters} \\
    bcftools view \\
        --no-version \\
        -Oz \\
        -o "${jobname}.vcf.gz" \\
        -
    bcftools index \\
        -t \\
        "${jobname}.vcf.gz"
    """
}


/*
 * STEP 11.1 - Concatenate the chunks of soft-filtered variants
 */
process concat {
    tag { sample }

    label 'bcftools'

    publishDir "${params.outdir}/BCFtools/${sample}", mode: 'copy', saveAs: { fn ->
        fn.endsWith(".stats.txt") ? null : "${fn}"
    }

    input:
    set sample, file(vcf_files), file("*") from mpileup_region_variants
        .groupTuple()
        .map { sample, indexed_vcf_files ->

            def vcf_files = indexed_vcf_files*.head()
            def tbi_files = indexed_vcf_files*.tail().flatten()

            tuple(sample.toString(), vcf_files, tbi_files)
        }

    output:
    set sample, file("${sample}.vcf.gz{,.tbi}") into soft_filtered_variants
    file "*.stats.txt" into raw_variant_stats

    script:
    def quoted_vcf_files = vcf_files.collect { /"${it}"/ }.join(' ')

    """
    bcftools concat \\
        --no-version \\
        -a \\
        -D \\
        -Oz \\
        -o "${sample}.vcf.gz" \\
        ${quoted_vcf_files}
    bcftools index \\
        -t \\
        "${sample}.vcf.gz"
    bcftools stats \\
        "${sample}.vcf.gz" \\
        > "${sample}.stats.txt"
    """
}


/*
 * STEP 12 - Hard filter for PASS variants
 */
process filter {
    tag { sample }

    label 'bcftools'

    publishDir "${params.outdir}/BCFtools/${sample}", mode: 'copy', saveAs: { fn ->
        fn.endsWith(".filtered.stats.txt") ? null : "${fn}"
    }

    input:
    set sample, file(indexed_vcf) from soft_filtered_variants

    output:
    set sample, file("${sample}.filtered.vcf.gz{,.tbi}") into filtered_variants
    file "*.filtered.stats.txt" into filtered_variant_stats

    script:
    def (vcf, tbi) = indexed_vcf

    """
    bcftools view \\
        --no-version \\
        -f PASS \\
        -o "${sample}.filtered.vcf.gz" \\
        -Oz \\
        "${vcf}"
    bcftools index \\
        -t \\
        "${sample}.filtered.vcf.gz"
    bcftools stats \\
        "${sample}.filtered.vcf.gz" \\
        > "${sample}.filtered.stats.txt"
    """
}


/*
 * STEP 13 - Annotate the filtered variants using the Ensembl Variant Effect Predictor
 */
process vep {
    tag { sample }

    label 'ensembl_vep'

    publishDir "${params.outdir}/VEP/${sample}", mode: 'copy'

    input:
    set sample, file(indexed_vcf) from filtered_variants
    file "cache/*" from indexed_vep_cache
    file ref_fasta
    file ref_faidx

    output:
    set sample, file("${sample}.filtered.vep.vcf.gz{,.tbi}") into vep_results
    file "*.html" into vep_stats

    script:
    def (vcf, tbi) = indexed_vcf
    def vep_cache_type = params.vep_cache_type ? '--' + params.vep_cache_type : ''

    """
    vep \\
        --everything \\
        --species "${params.vep_species_name}" \\
        --assembly "${params.vep_assembly_name}" \\
        --input_file "${vcf}" \\
        --output_file "${sample}.filtered.vep.vcf.gz" \\
        --stats_file "${sample}.filtered.stats.html" \\
        --fork ${task.cpus - 1} \\
        --dir cache \\
        --offline \\
        --fasta "${ref_fasta}" \\
        ${vep_cache_type} \\
        --vcf \\
        --compress_output bgzip \\
        --dont_skip \\
        --allow_non_variant
    tabix "${sample}.filtered.vep.vcf.gz"
    """
}


/*
 * STEP 13.1 - Parse the Ensembl VEP VCF into tab-separated values
 */
process vepvcf2csv {
    tag { sample }

    label 'ensembl_vep'

    publishDir "${params.outdir}/VEP/${sample}", mode: 'copy'

    input:
    set sample, file(indexed_vcf) from vep_results

    output:
    set sample, file("${sample}.filtered.vep.csv.gz")

    script:
    def (vcf, tbi) = indexed_vcf

    """
    vepvcf2csv.pl \\
        -i "${vcf}" \\
        -o "${sample}.filtered.vep.csv.gz"
    """
}


/*
 * STEP 14 - Run MultiQC
 */
process multiqc {

    label 'multiqc'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file config from multiqc_cfg
    file 'fastqc-raw/*' from fastqc_raw_results.collect()
    file 'cutadapt/*' from cutadapt_logs.collect()
    file 'fastqc-trimmed/*' from fastqc_trimmed_results.collect()
    file 'star/*' from star_logs.collect()
    file 'picard/markduplicates/*' from mark_duplicates_metrics.collect()
    file 'bcftools-raw/*' from raw_variant_stats.collect()
    file 'bcftools-filtered/*' from filtered_variant_stats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    """
    multiqc \\
        --config "${config}" \\
        -m fastqc \\
        -m cutadapt \\
        -m star \\
        -m picard \\
        -m bcftools \\
        .
    """
}


/* Workflow */

workflow.onComplete {

    log.info "Workflow completed at: ${workflow.complete}"
    log.info "Time taken: ${workflow.duration}"
    log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
    log.info "Output directory: ${params.outdir}"
}

workflow.onError {

    log.info "Execution halted: ${workflow.errorMessage}"
}


/* Functions */

def usage() {

    log.info"""
    Usage:
        nextflow run -profile <profile> -revision <revision> ampatchlab/nf-rnasnv [options]


    Nextflow execution options:

        -profile STR
            Nextflow configuration profile to use. Available profiles include:
            'awsbatch', 'conda', 'docker' and 'singularity'

        -revision STR
            Git branch/tag (version) of this workflow to use

        -work-dir DIR
            Directory where intermediate result files are stored

        -help
            Show additional execution options and exit


    Input options:

        --csv FILE
            Comma-separated list of sample and readgroup inputs

        --paired_end
            Expect entries for 'fastq1' and 'fastq2' in the input CSV


    Reference genome options:

        --genome STR
            Reference genome name [Either: ${defaults.genomes.keySet().join(", ")}; Default: ${defaults.genome}]

        --fasta FILE
            Override the reference genome FASTA with FILE [Default: ${defaults.fasta ?: null}]

        --gtf FILE
            Override the reference genome GTF with FILE [Default: ${defaults.gtf ?: null}]


    Sequencing adapter options:

        --adapters STR
            The adapters to trim [Either: ${defaults.seq_adapters.keySet().join(", ")}; Default: ${defaults.adapters}]

        --r1_adapter STR
            Override the sequence of the R1 adapter with STR [Default: ${defaults.r1_adapter ?: null}]

        --r2_adapter STR
            Override the sequence of the R2 adapter with STR [Default: ${defaults.r2_adapter ?: null}]


    Output options:

        --refdir DIR
            Path where the reference index files will be saved [Default: ${defaults.refdir}]

        --outdir DIR
            Path where the results will be saved [Default: ${defaults.outdir}]


    Standard options:

        --advanced
            Show advanced usage and exit

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}

def advanced() {

    log.info"""
    STAR genome generate options:

        --star_genome_chr_bin_n_bits INT
            Size of the bins for genome storage [Default: ${defaults.star_genome_chr_bin_n_bits}]

        --star_genome_sa_index_n_bases INT
            Length (bases) of the SA pre-indexing string [Default: ${defaults.star_genome_sa_index_n_bases}]

        --star_sjdb_overhang INT
            Length of the donor/acceptor sequence on each side of the junctions [Default: ${defaults.star_sjdb_overhang}]


    Cutadapt options:

        --cutadapt_base_qual_cutoff [INT,]INT
            Trim low-quality bases from each read [Default: ${defaults.cutadapt_base_qual_cutoff}]

        --cutadapt_min_read_length INT[:INT]
            Discard reads shorter than INT [Default: ${defaults.cutadapt_base_qual_cutoff}]


    BCFtools mpileup options:

        --mpileup_num_regions INT
            Number of regions per mpileup process [Default: ${defaults.mpileup_num_regions}]

        --mpileup_max_depth INT
            Maximum number of reads to pileup at a given position [Default: ${defaults.mpileup_max_depth}]

        --mpileup_min_bq INT
            Minimum base quality for a base to be considered [Default: ${defaults.mpileup_min_bq}]


    Ensembl VEP options:

        --vep_cache_type STR
            Alternate cache to use [Either: ${defaults.vep_cache_types.join(", ")}; Default: ${defaults.vep_cache_type}]

        --vep_indexed_cache_file FILE
            Override the Ensembl VEP indexed cache file with FILE [Default: ${defaults.vep_indexed_cache_file ?: null}]

        --vep_species_name STR
            Override the Ensembl VEP species name with STR [Default: ${defaults.vep_species_name ?: null}]

        --vep_assembly_name STR
            Override the Ensembl VEP assembly name with STR [Default: ${defaults.vep_assembly_name ?: null}]


    MultiQC options:

        --multiqc_config FILE
            MultiQC YAML config file [Default: ${defaults.multiqc_config}]


    Report options

        --execution_report STR
            Name of the Nextflow execution report to generate [Default: ${defaults.execution_report}]

        --trace_report STR
            Name of the Nextflow trace report to generate [Default: ${defaults.trace_report}]

        --timeline_report STR
            Name of the Nextflow timeline report to generate [Default: ${defaults.timeline_report}]

        --flowchart STR
            Name of the Nextflow flowchart to generate [Default: ${defaults.flowchart}]


    AWS Batch options

        --aws_queue STR
            AWS Batch JobQueue definition [Default: ${defaults.aws_queue}]

        --aws_region STR
            AWS Region definition [Default: ${defaults.aws_region}]
    """.stripIndent()
}

def die() {
    usage()
    exit 1
}

def check_params() {

    // Standard options

    if (params.advanced) {
        advanced()
        exit 0
    }

    if (params.help) {
        usage()
        exit 0
    }

    if (params.version) {
        log.info(workflow.manifest.version)
        exit 0
    }


    // Required options

    if (!params.csv) {
        log.error("A list of samples and readgroups is required. Please use the `--csv` option.")
        die()
    }

    if (file(params.csv).getExtension() != "csv") {
        log.error("Readgroup input file `${params.csv}` must be a CSV file with the '.csv' extension.")
        die()
    }


    // Reference genome options

    if (!params.genome) {
        log.error("Please specify a value for `--genome`; can be one of ${params.genomes.keySet().join(", ")}")
        die()
    }

    params.fasta = params.genome ? params.genomes[ params.genome ].fasta : null
    params.gtf = params.genome ? params.genomes[ params.genome ].gtf : null

    if (!params.fasta) {
        log.error("A reference FASTA file is required. Please use the `--fasta` option.")
        die()
    }

    if (!params.gtf) {
        log.error("A reference GTF file is required. Please use the `--gtf` option.")
        die()
    }


    // Sequencing adapter options

    if (!params.adapters) {
        log.error("Please specify a value for `--adapters`; can be one of ${params.seq_adapters.keySet().join(", ")}")
        die()
    }

    params.r1_adapter = params.adapters ? params.seq_adapters[ params.adapters ].r1 : null
    params.r2_adapter = params.adapters ? params.seq_adapters[ params.adapters ].r2 : null

    if (!params.r1_adapter) {
        log.error("An R1 adapter sequence is required. Please use the `--r1_adapter` option.")
        die()
    }

    if (!params.r2_adapter) {
        log.error("An R2 adapter sequence is required. Please use the `--r2_adapter` option.")
        die()
    }


    // STAR genome generate options

    if (!(params.star_genome_chr_bin_n_bits.toString().isInteger())) {
        log.error("Unknown `--star_genome_chr_bin_n_bits` entry: `${params.star_genome_chr_bin_n_bits}`")
        die()
    }

    if (!(params.star_genome_sa_index_n_bases.toString().isInteger())) {
        log.error("Unknown `--star_genome_sa_index_n_bases` entry: `${params.star_genome_sa_index_n_bases}`")
        die()
    }

    if (!(params.star_sjdb_overhang.toString().isInteger())) {
        log.error("Unknown `--star_sjdb_overhang` entry: `${params.star_sjdb_overhang}`")
        die()
    }


    // Cutadapt options

    if (!(params.cutadapt_base_qual_cutoff.toString().isInteger())) {

        if (params.cutadapt_base_qual_cutoff.toString().contains(',')) {
            def (five_prime_cutoff, three_prime_cutoff) = params.cutadapt_base_qual_cutoff.split(',', 2)

            if (!five_prime_cutoff.isInteger() || !three_prime_cutoff.isInteger()) {
                log.error("Unknown `--cutadapt_base_qual_cutoff` entry: `${params.cutadapt_base_qual_cutoff}`")
                die()
            }
        }
        else {
            log.error("Unknown `--cutadapt_base_qual_cutoff` entry: `${params.cutadapt_base_qual_cutoff}`")
            die()
        }
    }

    if (!(params.cutadapt_min_read_length.toString().isInteger())) {

        if (params.cutadapt_min_read_length.toString().contains(':')) {
            def (r1_min_length, r2_min_length) = params.cutadapt_min_read_length.split(':', 2)

            if (!r1_min_length.isInteger() || !r2_min_length.isInteger()) {
                log.error("Unknown `--cutadapt_min_read_length` entry: `${params.cutadapt_min_read_length}`")
                die()
            }
        }
        else {
            log.error("Unknown `--cutadapt_min_read_length` entry: `${params.cutadapt_min_read_length}`")
            die()
        }
    }


    // BCFtools mpileup options

    if (!(params.mpileup_num_regions.toString().isInteger())) {
        log.error("Unknown `--mpileup_num_regions` entry: `${params.mpileup_num_regions}`")
        die()
    }

    if (!(params.mpileup_max_depth.toString().isInteger())) {
        log.error("Unknown `--mpileup_max_depth` entry: `${params.mpileup_max_depth}`")
        die()
    }

    if (!(params.mpileup_min_bq.toString().isInteger())) {
        log.error("Unknown `--mpileup_min_bq` entry: `${params.mpileup_min_bq}`")
        die()
    }


    // Ensembl VEP options

    if (params.vep_cache_type && !(params.vep_cache_type in params.vep_cache_types)) {
        log.error("Unknown `--vep_cache_type` entry: `${params.vep_cache_type}`; can be one of ${params.vep_cache_types.join(", ")}")
        die()
    }

    params.vep_indexed_cache_file = params.vep_indexed_cache_files[params.genome][params.vep_cache_type ?: 'ensembl'].cache ?: null

    if (!params.vep_indexed_cache_file) {
        log.error("An Ensembl VEP indexed cache file is required. Please use the `--vep_indexed_cache_file` option.")
        die()
    }

    params.vep_species_name = params.vep_cache_info[params.genome].species ?: null

    if (!params.vep_species_name) {
        log.error("An Ensembl VEP species name is required. Please use the `--vep_species_name` option.")
        die()
    }

    params.vep_assembly_name = params.genome ?: null

    if (!params.vep_assembly_name) {
        log.error("An Ensembl VEP assembly name is required. Please use the `--vep_assembly_name` option.")
        die()
    }


    // MultiQC options

    if (!params.multiqc_config) {
        log.error("A configuration file for MultiQC is required. Please use the `--multiqc_config` option.")
        die()
    }

    if (file(params.multiqc_config).getExtension() != "yaml") {
        log.error("MultiQC config file `${params.multiqc_config}` must be a YAML file with the '.yaml' extension.")
        die()
    }


    // Report options

    if (!params.execution_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--execution_report` must end with '.html'")
        die()
    }

    if (!params.trace_report.toString().endsWith('.txt')) {
        log.error("The filename specified using `--trace_report` must end with '.txt'")
        die()
    }

    if (!params.timeline_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--timeline_report` must end with '.html'")
        die()
    }

    def flowchart_extns = ['.dot', '.html', '.pdf', '.png', '.svg']

    if (!(flowchart_extns.any { params.flowchart.toString().endsWith(it) })) {
        log.error("The filename specified using `--flowchart` must end with one of ${flowchart_extns.join(", ")}")
        die()
    }
}

def validate_input_csv() {

    def csv = new CsvSplitter().target(file(params.csv)).options(header:true)
    def rows = csv.list()

    def fastq_columns = params.paired_end ? ["fastq1", "fastq2"] : ["fastq"]
    def required_columns = ["sample"] + fastq_columns

    required_columns.each { col ->

        if (!csv.columnsHeader.contains(col)) {
            log.error("Readgroup input file `${params.csv}` does not contain a '${col}' column. Exiting.")
            exit 1
        }
    }

    def readgroups = [:].withDefault { [] }
    def fastq_files = [:].withDefault { [] }

    def valid_file_extns = ["fastq", "fq", "fastq.gz", "fq.gz"]

    log.info("Validating ${rows.size()} entries...")

    rows.indexed(1).each { idx, row ->

        log.info("Validating entry ${idx}...")

        if (!row.sample) {
            log.error("Entry ${idx}: Invalid 'sample' value. Exiting.")
            exit 1
        }

        if (row.readgroup in readgroups[row.sample]) {
            log.error("Entry ${idx}: Sample `${row.sample}` must have unique readgroup entries. Exiting.")
            exit 1
        }

        readgroups[row.sample].add(row.readgroup)

        fastq_columns.each { fastq ->

            if (!row."${fastq}") {
                log.error("Entry ${idx}: Invalid '${fastq}' value. Exiting.")
                exit 1
            }

            fq = file(row."${fastq}")
            fq_extn = get_fastq_extn(fq)

            if (!valid_file_extns.contains(fq_extn)) {
                log.error("Entry ${idx}: ${fq} does not have a valid file extension. Must be one of: ${valid_file_extns.join(", ")}")
                exit 1
            }

            if (fq.getName() in fastq_files[row.sample]) {
                log.error("Entry ${idx}: File '${fq.getName()}' cannot be used more than once. Exiting.")
                exit 1
            }

            fastq_files[row.sample].add(fq.getName())
        }
    }

    log.info("Done")
}

def get_fastq_extn(fq) {
    def extn
    switch (fq) {
        case { it.name.endsWith('.fastq') }:
            extn = 'fastq'
            break
        case { it.name.endsWith('.fq') }:
            extn = 'fq'
            break
        case { it.name.endsWith('.fastq.gz') }:
            extn = 'fastq.gz'
            break
        case { it.name.endsWith('.fq.gz') }:
            extn = 'fq.gz'
            break
        default:
            extn = ''
            break
    }
    extn
}
