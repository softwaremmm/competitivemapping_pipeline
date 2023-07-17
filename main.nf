#!/usr/bin/env nextflow

//Set DSL2 syntax

nextflow.enable.dsl = 2

//Define ANSI colours for ease

ANSI_GREEN = '\033[1;32m'
ANSI_RESET = '\033[0m'

params.help = ''
params.input_dir = ''
params.output_dir = ''

if (workflow.profile != 'kubernetes') {
    params.knowledge_bucket = "$projectDir/data/relatedness/knowledge"
} else {
    params.knowledge_bucket = '/data/relatedness/knowledge'
}

include { competitiveMapping } from './process/competitive_mapping.nf'

params.help = ''

//Constants
fastq_pattern = '*_*{1,2}.f*q*'

workflow competitive_mapping {
    take:
    input_dir
    output_dir

    main:

        if (params.input_dir == '') {
        exit 1, 'error: --input_dir is mandatory'
        }

        if (params.output_dir == '') {
        exit 1, 'error: --output_dir is mandatory'
        }

        inputdir_amended = "${params.input_dir}".replaceFirst(/$/, '/')
        indir = "${inputdir_amended}"
        reads = indir + fastq_pattern

        Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
                .ifEmpty { error "cannot find any reads matching ${fastq_pattern} in ${indir}" }
                .set { input_files }
        input_files.view { it }

        competitiveMapping(input_files)
}

workflow.onComplete {
    if (workflow.success) {
        log.info '''
        ===========================================
        Workflow completed successfully
        '''
        .stripIndent()
    }
    else {
        log.info '''
        ===========================================
        Finished with errors
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
                --output_dir Directory for output

                '''
                .stripIndent()
        exit(0)
        }

        log.info """
        ========================================================================

        Competitive Mapping

        Determination of species by Competitive Mapping using minimap2.

        Parameters:
        ------------------------------------------------------------------------

        --input_dir    $params.input_dir
        --output_dir   $params.output_dir

        Runtime data:
        ------------------------------------------------------------------------

        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        Project directory     ${ANSI_GREEN}${projectDir}${ANSI_RESET}
        """
        .stripIndent()

        competitive_mapping(params.input_dir, params.output_dir)
}
