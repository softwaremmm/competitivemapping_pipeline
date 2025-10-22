#!/usr/bin/env nextflow
include { competitiveMapping } from './process/competitive_mapping.nf'
include { tie_break } from './process/competitive_mapping.nf'
include { tie_break_multi } from './process/competitive_mapping.nf'
include { dynamic_tie_break } from './process/competitive_mapping.nf'
include { has_enough_reads } from './process/competitive_mapping.nf'
include { filter_by_depth } from './process/competitive_mapping.nf'
include { extract_reads } from './process/competitive_mapping.nf'

// input parameters
params.input_dir = ''
params.manifest = ''
params.seq_platform = ''
params.reference_name = ''
params.workflow = 'comp_mapping'

// default thresholds
params.illumina_threshold = 100000
params.ont_threshold = 1000
params.use_whole_genera_in_dynamic_cm = true

// input defaults
params.input_paired_suffix = "*_{1,2}.fastq.gz"
params.input_single_suffix = "*.fastq.gz"

workflow {
    //Define ANSI colours for ease
    ANSI_GREEN = '\033[1;32m'
    ANSI_RESET = '\033[0m'

    if (params.help) {
        log.info(
            '''
            ========================================================================
            Competitive Mapping

            Determination of species by Competitive Mapping using minimap2.

            Parameters:
            ------------------------------------------------------------------------

            --input_dir  Directory holding the fastq files *_{1,2}.fastq.gz
            --manifest
            --species_list
            --seq_platform
            --reference_name: Reference name to use for the fastq files. Default: ''. If you want to use a list of names, input a string separated by commas.
            --workflow: Workflow to run. Options are 'comp_mapping' (default) or 'tie_break'
            '''.stripIndent()
        )
        exit(0)
    }

    if (params.input_dir == '') {
        exit(1, 'error: --input_dir is mandatory')
    }
    if (params.manifest == '') {
        exit(1, 'error: --manifest is mandatory')
    }
    if (params.species_list == '') {
        exit(1, 'error: --species_list is mandatory')
    }
    if (params.seq_platform == '') {
        exit(1, 'error: --seq_platform is mandatory')
    }


    log.info(
        """
        ========================================================================

        Competitive Mapping

        Determination of species by Competitive Mapping using minimap2.

        Parameters:
        ------------------------------------------------------------------------

        --input_dir    ${params.input_dir}
        --manifest     ${params.manifest}
        --species_list ${params.species_list}
        --seq_platform ${params.seq_platform}
        --reference_name ${params.reference_name}
        --workflow     ${params.workflow}

        Runtime data:
        ------------------------------------------------------------------------

        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        Project directory     ${ANSI_GREEN}${projectDir}${ANSI_RESET}
        """.stripIndent()
    )

    if (params.seq_platform == 'ont') {
        input_files = Channel.fromPath("${params.input_dir}/${params.input_single_suffix}", checkIfExists: true)
            .ifEmpty { error("cannot find any reads matching ${params.input_single_suffix} in ${params.input_dir}") }
            .map { it -> tuple(it.simpleName, it) }
    }
    else if (params.seq_platform == 'illumina') {
        input_files = Channel.fromFilePairs(
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

    // Using fromPath means they can be provided as relative paths
    manifest = Channel.fromPath(params.manifest, checkIfExists: true).first()
    species_list = Channel.fromPath(params.species_list, checkIfExists: true).first()

    if (params.workflow == 'comp_mapping') 
        competitive_mapping(input_files, manifest, species_list, params.seq_platform, params.reference_name)
    else if (params.workflow == 'tie_break')
        tie_break_workflow(input_files, manifest, species_list, params.seq_platform, params.reference_name)
    else
        exit(1, "error: --workflow must be one of 'comp_mapping' or 'tie_break'")
}


workflow competitive_mapping {
    take:
    input_files
    manifest
    species_list
    seq_platform
    reference_name

    main:

    check_seq_platform(seq_platform)

    competitive_mapping_output = competitiveMapping(input_files, manifest, species_list, seq_platform, reference_name)
    threshold = seq_platform == 'illumina' ? params.illumina_threshold : params.ont_threshold
    has_enough_reads(competitive_mapping_output.report_json, threshold)

    emit:
    ref_reads = competitive_mapping_output.ref_reads
    report_json = competitive_mapping_output.report_json
    report_csv = competitive_mapping_output.report_csv
    cm_enough_reads = has_enough_reads.out
}

// WARNING: Experimental process
workflow tie_break_workflow {
    take:
    input_files
    manifest
    species_list
    seq_platform
    reference_name

    main:
    check_seq_platform(seq_platform)

    analyzer_params = Channel.fromPath("${moduleDir}/process/params_${seq_platform}.yml").first()

    tie_break(
        input_files,
        manifest,
        species_list,
        seq_platform,
        analyzer_params,
        reference_name,
    )

    emit:
    report_csv = tie_break.out.report_csv
    stats = tie_break.out.stats
    ref_reads = tie_break.out.ref_reads
}

// WARNING: Experimental process
workflow tie_break_multi_workflow {
    take:
    input_files
    manifest
    species_list
    seq_platform
    reference_name

    main:
    check_seq_platform(seq_platform)

    analyzer_params = Channel.fromPath("${moduleDir}/process/params_${seq_platform}.yml").first()

    tie_break_multi(
        input_files,
        manifest,
        species_list,
        seq_platform,
        analyzer_params,
        reference_name,
    )

    tie_break.out.report_csv.view()

    filter_by_depth(tie_break.out.report_csv, species_list, 5)

    filter_by_depth.out.high_depth_list.view()

    high_depth_refs_ch = filter_by_depth.out.high_depth_list.flatMap { sample_name, file ->
        file.text
            .readLines()
            .collect { reference ->
                tuple(sample_name, reference)
            }
    }

    high_depth_refs_ch.view { "High depth references: ${it}" }

    high_depth_accessions_ch = filter_by_depth.out.high_depth_accessions.flatMap { sample_name, file ->
        file.text
            .readLines()
            .collect { accession ->
                tuple(sample_name, accession)
            }
    }

    high_depth_accessions_ch.view { "High depth accessions: ${it}" }

    high_depth_refs_ch = input_files
        .combine(high_depth_refs_ch, by: 0)
        .combine(tie_break.out.round_two_alignments, by: 0)
        .combine(tie_break.out.references, by: 0)
        .combine(high_depth_accessions_ch, by: 0)

    high_depth_refs_ch.view { "Input to extract_reads process: ${it}" }

    extract_reads(high_depth_refs_ch)

    extract_reads.out.ref_reads.view { "Extracted reads: ${it}" }

    emit:
    report_csv = tie_break.out.report_csv
    stats = tie_break.out.stats
    mapped_reads = extract_reads.out.ref_reads
}

// WARNING: Experimental process
// currently also runs standard workflow for comparison
workflow dynamic_tie_break_workflow {
    take:
    input_files
    genome_dirs // Each directory must contain a genome_paths.tsv file
    assembly_metadata
    seq_platform

    main:
    check_seq_platform(seq_platform)

    analyzer_params = Channel.fromPath("${moduleDir}/process/params_${seq_platform}.yml").first()

    dynamic_tie_break(
        input_files,
        genome_dirs,
        assembly_metadata,
        seq_platform,
        params.use_whole_genera_in_dynamic_cm,
        analyzer_params,
    )

    emit:
    report_csv = dynamic_tie_break.out.report_csv
    stats = dynamic_tie_break.out.stats
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
