#!/usr/bin/env nextflow

//Define ANSI colours for ease

ANSI_GREEN = '\033[1;32m'
ANSI_RESET = '\033[0m'

params.help = ''
params.input_dir = ''
params.manifest = ''
params.seq_platform = ''

params.illumina_threshold = 100000
params.ont_threshold = 1000

include { competitiveMapping } from './process/competitive_mapping.nf'
include { has_enough_reads } from './process/competitive_mapping.nf'

//Constants
input_paired_suffix = "*_{1,2}.fastq.gz"
input_single_suffix = "*.fastq.gz"
seq_platforms = ['ont', 'illumina']


workflow competitive_mapping {
    take:
        input_files
        manifest
        species_list
        seq_platform

    main:

    // seq_platform should be a String, not a channel
    if (seq_platform.getClass() != java.lang.String) {
        throw new Exception("seq_platform should be a string, not a ${seq_platform.getClass()}")
    }

    // Should be supported by this workflow
    if (! (seq_platform in seq_platforms)) {
        throw new Exception("seq platform invalid. Should be one of $seq_platforms!")
    }

    // WARNING: Previous version used params for manifest and species_list, 
    // instead of using input channels
    competitive_mapping_output = competitiveMapping(input_files, manifest, species_list, seq_platform)
    threshold = seq_platform == 'illumina' ? params.illumina_threshold : params.ont_threshold
    has_enough_reads(competitive_mapping_output.cm_report, threshold)

    emit:
        cm_sample_paths = competitive_mapping_output.cm_sample
        cm_report = competitive_mapping_output.cm_report
        cm_error = competitive_mapping_output.cm_error
        cm_enough_reads =  has_enough_reads.out
}

workflow.onComplete {
    if (workflow.success) {
        log.info '''
        ===========================================
        Competitive Mapping Workflow completed successfully
        '''
        .stripIndent()
    }
    else {
        log.info '''
        ===========================================
        Competitive Mapping finished with errors
        '''
        .stripIndent()
    }
}

workflow {
    main:
        if (params.help) {
        log.info '''
                ========================================================================
                Competitive Mapping

                Determination of species by Competitive Mapping using minimap2.

                Parameters:
                ------------------------------------------------------------------------

                --input_dir  Directory holding the fastq files *_{1,2}.fastq.gz
                --manifest
                --species_list
                --seq_platform
                '''
                .stripIndent()
        exit(0)
        }

        if (params.input_dir == '') {
            exit 1, 'error: --input_dir is mandatory'
        }
        if (params.manifest == '') {
            exit 1, 'error: --manifest is mandatory'
        }
        if (params.species_list == '') {
            exit 1, 'error: --species_list is mandatory'
        }
        if (params.seq_platform == '') {
            exit 1, 'error: --seq_platform is mandatory'
        }


        log.info """
        ========================================================================

        Competitive Mapping

        Determination of species by Competitive Mapping using minimap2.

        Parameters:
        ------------------------------------------------------------------------

        --input_dir    $params.input_dir
        --manifest     $params.manifest
        --species_list $params.species_list
        --seq_platform $params.seq_platform

        Runtime data:
        ------------------------------------------------------------------------

        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        Project directory     ${ANSI_GREEN}${projectDir}${ANSI_RESET}
        """
        .stripIndent()

        if (params.seq_platform == 'ont') {
            input_files = Channel.fromPath("${params.input_dir}/${input_single_suffix}", checkIfExists: true)
                .ifEmpty { error "cannot find any reads matching ${input_single_suffix} in ${params.input_dir}" }
                .map { it -> tuple(it.simpleName, it)}
                .first()
        }
        else if (params.seq_platform == 'illumina') {
            input_files = Channel
                .fromFilePairs("${params.input_dir}/${input_paired_suffix}",
                    flat: false,
                    checkIfExists: true,
                    size: -1)
                .ifEmpty { error "cannot find any reads matching ${input_paired_suffix} in ${params.input_dir}" }
                .first()
        }

        input_files.view()

        manifest = Channel.fromPath(params.manifest, checkIfExists: true)
        species_list = Channel.fromPath(params.species_list, checkIfExists: true)
        competitive_mapping(input_files, manifest, species_list, params.seq_platform)
}
