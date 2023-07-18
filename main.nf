#!/usr/bin/env nextflow

include {competitiveMapping} from './process/competitive_mapping.nf'

params.help = ""
params.input_dir = ""

//Constants
fastq_pattern = '*_*{1,2}.f*q*'

workflow competitive_mapping{
    take:
    input_dir
    output_dir

    main:
        // Setup so --help triggers the help message
        if (params.help) {
            log.info """
            ========================================================================
            Competitive Mapping Workflow

            Parameters:
            ------------------------------------------------------------------------
            --input_dir    Path to the sample's directory
            --output_dir   Path to the workflow's output directory

            """
            .stripIndent()
            exit(0)
        }

        if (params.input_dir == '') {
            exit 1, 'error: --input_dir is mandatory'
        }


        inputdir_amended = "${params.input_dir}".replaceFirst(/$/, '/')
        indir = "${inputdir_amended}"
        reads = indir + fastq_pattern

        Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
                .ifEmpty { error "cannot find any reads matching ${fastq_pattern} in ${indir}" }
                .set { input_files }
        input_files.view { it }

        competitive_mapping_output = competitiveMapping(input_files)

    emit:
        cm_sample_paths = competitive_mapping_output.cm_sample
        cm_report = competitive_mapping_output.cm_report
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


workflow{

    main:
        competitive_mapping(params.input_dir, params.output_dir)
}



