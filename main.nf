#!/usr/bin/env nextflow
include { competitiveMapping } from './process/competitive_mapping.nf'
include { dynamicCompetitiveMapping } from './process/competitive_mapping.nf'
include { has_enough_reads } from './process/competitive_mapping.nf'
include { contigMapping } from './process/contig_mapping.nf'

//Define parameters
params.help = ''
params.input_dir = ''
params.manifest = ''
params.seq_platform = ''

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

        Runtime data:
        ------------------------------------------------------------------------

        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        Project directory     ${ANSI_GREEN}${projectDir}${ANSI_RESET}
        """.stripIndent()
    )

    if (params.seq_platform == 'ont') {
        input_files = Channel
            .fromPath("${params.input_dir}/${params.input_single_suffix}", checkIfExists: true)
            .ifEmpty { error("cannot find any reads matching ${params.input_single_suffix} in ${params.input_dir}") }
            .map { it -> tuple(it.simpleName, it) }
            .first()
    }
    else if (params.seq_platform == 'illumina') {
        input_files = Channel
            .fromFilePairs(
                "${params.input_dir}/${params.input_paired_suffix}",
                flat: false,
                checkIfExists: true,
                size: -1
            )
            .ifEmpty { error("cannot find any reads matching ${params.input_paired_suffix} in ${params.input_dir}") }
            .first()
    }

    input_files.view()

    manifest = Channel.fromPath(params.manifest, checkIfExists: true)
    species_list = Channel.fromPath(params.species_list, checkIfExists: true)
    competitive_mapping(input_files, manifest, species_list, params.seq_platform)
}


workflow competitive_mapping {
    take:
    input_files
    manifest
    species_list
    seq_platform

    main:

    check_seq_platform(seq_platform)

    competitive_mapping_output = competitiveMapping(input_files, manifest, species_list, seq_platform)
    threshold = seq_platform == 'illumina' ? params.illumina_threshold : params.ont_threshold
    has_enough_reads(competitive_mapping_output.cm_report, threshold)

    emit:
    cm_tb_reads = competitive_mapping_output.cm_tb_reads
    cm_report = competitive_mapping_output.cm_report
    cm_csv = competitive_mapping_output.cm_csv
    cm_enough_reads = has_enough_reads.out
}

workflow dynamic_competitive_mapping {
    take:
    input_files
    gtdb_genomes_path
    assembly_metadata
    seq_platform

    main:

    check_seq_platform(seq_platform)

    competitive_mapping_output = dynamicCompetitiveMapping(
        input_files,
        gtdb_genomes_path,
        assembly_metadata,
        seq_platform,
        params.use_whole_genera_in_dynamic_cm
    )

    emit:
    cm_report = competitive_mapping_output.cm_report
    cm_csv = competitive_mapping_output.cm_csv
    manifest = competitive_mapping_output.manifest
}

workflow dynamic_contig_mapping {
    take:
    contig_files
    manifest_files

    main:
    input_files = contig_files.join(manifest_files)
    contigMapping(input_files)


    emit:
    species_comparison = contigMapping.species_comparison
    blast_mapping = contigMapping.blast_mapping
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
