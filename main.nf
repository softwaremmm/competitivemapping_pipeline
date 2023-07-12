
//Constants
fastq_pattern = '*_*{1,2}.f*q*'

if (params.input_dir == '') {
  exit 1, 'error: --input_dir is mandatory'
}

workflow {
  inputdir_amended = "${params.input_dir}".replaceFirst(/$/, '/')
  indir = "${inputdir_amended}"
  reads = indir + fastq_pattern

  Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
        .ifEmpty { error "cannot find any reads matching ${fastq_pattern} in ${indir}" }
        .set { input_files }
  input_files.view { it } // print channel contents to console
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
