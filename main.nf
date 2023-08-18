#!/usr/bin/env nextflow

//Define ANSI colours for ease

ANSI_GREEN = '\033[1;32m'
ANSI_RESET = '\033[0m'

params.help = ''
params.input_dir = ''
params.manifest_dir = ''
params.manifest_filename = 'manifest_v2.fasta'

include { competitiveMapping } from './process/competitive_mapping.nf'
include { has_enough_reads } from './process/competitive_mapping.nf'

//Constants
fastq_pattern = "*{1,2}.f*q.gz"




workflow competitive_mapping {
    take:
        input_files
        manifest

    main:
  
    Channel.fromPath(params.manifest)
        .set{ manifest_ch }

    competitive_mapping_output = competitiveMapping(input_files, manifest_ch)

    emit:
        cm_sample_paths = competitive_mapping_output.cm_sample
        cm_report = competitive_mapping_output.cm_report
        cm_enough_reads =  has_enough_reads(cm_report, 100000)
        cm_enough_reads.view{it}
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

                '''
                .stripIndent()
        exit(0)
        }

        if (params.input_dir == '') {
            exit 1, 'error: --input_dir is mandatory'
        }
        if (params.manifest_dir == '') {
            exit 1, 'error: --manifest_dir is mandatory'
        }
        if (params.manifest_file == '') {
            exit 1, 'error: --manifest_file is mandatory'
        }

        log.info """
        ========================================================================

        Competitive Mapping

        Determination of species by Competitive Mapping using minimap2.

        Parameters:
        ------------------------------------------------------------------------

        --input_dir          $params.input_dir
        --manifest_dir       $params.manifest_dir
        --manifest_filename  $params.manifest_filename

        Runtime data:
        ------------------------------------------------------------------------

        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        Project directory     ${ANSI_GREEN}${projectDir}${ANSI_RESET}
        """
        .stripIndent()

        Channel.fromFilePairs("${params.input_dir}/${fastq_pattern}", flat: true, checkIfExists: true, size: -1)
                .ifEmpty { error "cannot find any reads matching ${fastq_pattern} in ${params.input_dir}" }
                .set { input_files }
        input_files.view { it }

        competitive_mapping(input_files,params.manifest_dir,params.manifest_filename)
}
