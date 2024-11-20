
process competitiveMapping {
    container = {
        params.test_container=="" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:1.3.1' : params.test_container
    }

    cpus = 4
    memory = {
        params.testing=="" ? {12.GB + (36.GB * (task.attempts - 1))} : "16GB"
    }

    debug true
    pod label: "name", value: "competitive_mapping_pipeline:competitiveMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs)
    path (manifest)
    path (species_list)
    val(seq_platform)

    output:
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: cm_tb_reads
    tuple val(sample_name), path(competitive_mapping_report), emit: cm_report
    tuple val(sample_name), path(competitive_mapping_csv), emit: cm_csv

    script:
    tb_reads = "reads_for_assembly.fastq.gz"
    tb_reads_1 = "reads_for_assembly_1.fastq.gz"
    tb_reads_2 = "reads_for_assembly_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
    h37rv_ref="M.tuberculosis"
    """
    mkdir outputs
    competitive_mapping manifest --seq_platform ${seq_platform} \
        --reads ${fqs} \
        --ref_for_fastq ${h37rv_ref} \
        --manifest ${manifest} \
        --contigs ${species_list} \
        --cpus ${task.cpus} \
        --output_root "out."
    
    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}

    if [ ${seq_platform} == 'ont' ]
    then
        mv out.reads.fastq.gz ${tb_reads}
    elif [ ${seq_platform} == 'illumina' ]
    then
        mv out.reads_1.fastq.gz ${tb_reads_1}
        mv out.reads_2.fastq.gz ${tb_reads_2}
    fi
    """
}

process dynamicCompetitiveMapping {
    container = {
        params.test_container=="" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:1.3.1' : params.test_container
    }

    cpus 8
    memory { 8.GB * task.attempt}

    pod label: "name", value: "competitive_mapping_pipeline:dynamicCompetitiveMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs), path(sylph_report)
    path (gtdb_genomes)
    val(seq_platform)

    output:
    tuple val(sample_name), path(competitive_mapping_report), emit: cm_report
    tuple val(sample_name), path(competitive_mapping_csv), emit: cm_csv

    script:
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
    h37rv_ref="M.tuberculosis"
    """
    mkdir outputs
    competitive_mapping sylph --seq_platform ${seq_platform} \
        --reads ${fqs} \
        --sylph_report ${sylph_report} \
        --genomes ${gtdb_genomes} \
        --cpus ${task.cpus} \
        --output_root "out."
    
    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}
    """
}

process has_enough_reads {
    container = {
        params.test_container=="" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:1.3.1' : params.test_container
    }

    cpus = 1
    memory = "128MB"

    debug true
    pod label: "name", value: "competitive_mapping_pipeline:has_enough_reads"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path (json)
    val (threshold)

    output:
    tuple val(sample_name), stdout

    script:
    """
    check_read_count --json_file_path ${json} --read_threshold ${threshold}
    """
}
