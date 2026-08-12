#!/usr/bin/env nextflow
include { competitive_mapping } from './process/competitive_mapping.nf'
include { build_manifest } from './process/competitive_mapping.nf'
include { has_enough_reads } from './process/competitive_mapping.nf'
include { sylph } from './process/sylph.nf'

// input parameters
params.input_dir = ''
params.manifest = ''
params.species_list = ''
params.seq_platform = ''
params.ref_for_fastqs = ''
params.workflow = 'comp_mapping'

// default thresholds
params.illumina_threshold = 100000
params.ont_threshold = 1000
params.query_ani_threshold = 90
params.profile_ani_threshold = 95

// input defaults
params.input_paired_suffix = "*_{1,2}.fastq.gz"
params.input_single_suffix = "*.fastq.gz"

workflow {

    // Get input fastqs
    if (params.seq_platform == 'ont') {
        input_files = channel.fromPath("${params.input_dir}/${params.input_single_suffix}", checkIfExists: true)
            .ifEmpty { error("cannot find any reads matching ${params.input_single_suffix} in ${params.input_dir}") }
            .map { it -> tuple(it.simpleName, it) }
    }
    else if (params.seq_platform == 'illumina') {
        input_files = channel.fromFilePairs(
                "${params.input_dir}/${params.input_paired_suffix}",
                flat: false,
                checkIfExists: true,
                size: -1,
            )
            .map { it -> tuple(it[0].replace("_mycobacterial_reads", ""), it[1]) }
            .ifEmpty { error("cannot find any reads matching ${params.input_paired_suffix} in ${params.input_dir}") }
    }

    // Show first 3 in channel so user can check if they are correct
    input_files.take(3).view()

    if (params.workflow == 'comp_mapping') {
        // Using fromPath means they can be provided as relative paths
        manifest = channel.fromPath(params.manifest, checkIfExists: true).first()
        species_list = channel.fromPath(params.species_list, checkIfExists: true).first()

        competitive_mapping_wf(
            input_files,
            manifest,
            species_list,
            params.seq_platform,
            params.refs_for_fastqs,
        )
    }
    else if (params.workflow == 'dynamic') {

        // Some params could be provided as comma-separated strings or lists, so we handle both cases here
        ref_genome_dirs_list = params.ref_genome_dirs instanceof List ? params.ref_genome_dirs : params.ref_genome_dirs.split(',').collect { it -> it.trim() }
        ref_genome_dirs_ch = channel.fromList(ref_genome_dirs_list).map { it -> file(it) }.collect()

        sylph_dbs_list = params.sylph_dbs instanceof List ? params.sylph_dbs : params.sylph_dbs.split(',').collect { it -> it.trim() }
        sylph_dbs_ch = channel.fromList(sylph_dbs_list).map { it -> file(it) }.collect()

        taxonomy_files_list = params.taxonomy_files instanceof List ? params.taxonomy_files : params.taxonomy_files.split(',').collect { it -> it.trim() }
        taxonomy_files_ch = channel.fromList(taxonomy_files_list).map { it -> file(it) }.collect()

        dynamic_competitive_mapping_wf(
            input_files,
            ref_genome_dirs_ch,
            sylph_dbs_ch,
            taxonomy_files_ch,
            params.fixed_refs,
            params.refs_for_fastqs,
            params.seq_platform,
        )
    }
    else {
        exit(1, "error: --workflow must be one of 'comp_mapping' or 'dynamic'")
    }
}


workflow competitive_mapping_wf {
    take:
    input_files // tuple of sample name and list of fastq paths
    manifest
    species_list
    seq_platform
    refs_for_fastqs

    main:

    check_seq_platform(seq_platform)

    // converting strings/paths to channels if needed (allows for more flexible input)
    if (manifest instanceof Path || manifest instanceof String) {
        print("converting manifest to channel\n")
        manifest = channel.fromPath(manifest, checkIfExists: true).first()
    }
    if (species_list instanceof Path || species_list instanceof String) {
        print("converting species_list to channel\n")
        species_list = channel.fromPath(species_list, checkIfExists: true).first()
    }

    input_files_with_manifest = input_files
        .combine(manifest)
        .combine(species_list)

    competitive_mapping(input_files_with_manifest, seq_platform, refs_for_fastqs, "false")

    // The following will group reads by reference producing a channel of tuples (sample_name, ref_name, fastqs)
    // where fastqs is a list for illumina and a single file for ont
    ref_reads = competitive_mapping.out.ref_reads.flatMap { sample_name, reads ->
        def read_list = (reads instanceof List) ? reads : [reads]
        read_list
            .groupBy { it -> it.name.replaceFirst(/_reads_for_assembly.*$/, '') }
            .collect { ref, ref_reads ->
                def output_reads = ref_reads.size() == 1 ? ref_reads[0] : ref_reads
                tuple(sample_name, ref, output_reads)
            }
    }

    // Then apply the has_enough_reads process to each tuple
    // result is a channel of tuples (sample_name, ref_name, fastqs, has_enough_reads)
    threshold = seq_platform == 'illumina' ? params.illumina_threshold : params.ont_threshold
    has_enough_reads(
        competitive_mapping.out.report_json.combine(
            ref_reads.map { sample_name, ref_name, _reads -> tuple(sample_name, ref_name) },
            by: 0
        ),
        threshold,
    )
    ref_reads = ref_reads.join(has_enough_reads.out, by: [0,1])


    emit:
    report_json = competitive_mapping.out.report_json
    report_csv = competitive_mapping.out.report_csv
    ref_reads = ref_reads
}


workflow dynamic_competitive_mapping_wf {
    take:
    input_files // tuple of sample name and list of fastq paths
    ref_genome_dirs // Each directory must contain a genome_paths.tsv file
    sylph_dbs // Paths .syldb files
    taxonomy_files // Paths to taxonomy tsv filse
    fixed_refs // optional comma-separated list of accessions to always include as references
    refs_for_fastqs // optional reference to extract reads for from comp mapping
    seq_platform

    main:
    check_seq_platform(seq_platform)

    sylph(
        input_files,
        sylph_dbs,
        taxonomy_files,
        seq_platform,
        params.query_ani_threshold,
        params.profile_ani_threshold,
    )

    build_manifest(
        sylph.out.sylph_report,
        ref_genome_dirs,
        taxonomy_files,
        fixed_refs,
    )

    input_files_with_manifest = input_files
        .join(build_manifest.out.manifest)
        .join(build_manifest.out.contigs)

    competitive_mapping(input_files_with_manifest, seq_platform, refs_for_fastqs, "true")

    // The following will group reads by reference producing a channel of tuples (sample_name, ref_name, fastqs)
    // where fastqs is a list for illumina and a single file for ont
    ref_reads = competitive_mapping.out.ref_reads.flatMap { sample_name, reads ->
        def read_list = (reads instanceof List) ? reads : [reads]
        read_list
            .groupBy { it -> it.name.replaceFirst(/_reads_for_assembly.*$/, '') }
            .collect { ref, ref_reads ->
                def output_reads = ref_reads.size() == 1 ? ref_reads[0] : ref_reads
                tuple(sample_name, ref, output_reads)
            }
    }

    // Then apply the has_enough_reads process to each tuple
    // result is a channel of tuples (sample_name, ref_name, fastqs, has_enough_reads)
    threshold = seq_platform == 'illumina' ? params.illumina_threshold : params.ont_threshold
    has_enough_reads(
        competitive_mapping.out.report_json.combine(
            ref_reads.map { sample_name, ref_name, _reads -> tuple(sample_name, ref_name) },
            by: 0
        ),
        threshold,
    )
    ref_reads = ref_reads.join(has_enough_reads.out, by: [0,1])

    emit:
    report_csv = competitive_mapping.out.report_csv
    report_json = competitive_mapping.out.report_json
    ref_reads = ref_reads
    sylph_report = sylph.out.sylph_report
    sylph_query = sylph.out.sylph_query
    sylph_taxonomy_report = sylph.out.taxonomy_report
    depth_plot = competitive_mapping.out.depth_plot
}

def check_seq_platform(seq_platform) {
    def seq_platforms = ['ont', 'illumina']

    if (seq_platform.getClass() != java.lang.String) {
        throw new Exception("seq_platform should be a string, not a ${seq_platform.getClass()}")
    }

    if (!(seq_platform in seq_platforms)) {
        throw new Exception("seq platform invalid. Should be one of ${seq_platforms}!")
    }
}
