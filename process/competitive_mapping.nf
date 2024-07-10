
process competitiveMapping {
    container = {
        params.test_container=="" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:1.2.6' : params.test_container
    }

    cpus = 8
    memory = {
        params.testing=="" ? "48GB" : "16GB"
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
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: cm_sample
    path("species_comparison_report.json"), emit: cm_report
    path("species_comparison_error.json"), emit: cm_error

    script:
    competitive_mapping_file = "reads_for_assembly_0.fastq.gz"
    competitive_mapping_file_1 = "reads_for_assembly_0_1.fastq.gz"
    competitive_mapping_file_2 = "reads_for_assembly_0_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_error = "species_comparison_error.json"
    h37rv_rname="AL123456.3"

    template "run_competitive_mapping.sh"

    stub:
    competitive_mapping_file = "reads_for_assembly_0.fastq.gz"
    competitive_mapping_file_1 = "reads_for_assembly_0_1.fastq.gz"
    competitive_mapping_file_2 = "reads_for_assembly_0_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_error = "species_comparison_error.json"

    """
    if [ $seq_platform == 'ont' ]
    then
        touch ${competitive_mapping_file}
    elif [ $seq_platform == 'illumina' ]
    then
        touch ${competitive_mapping_file_1}
        touch ${competitive_mapping_file_2}
    fi
    touch ${competitive_mapping_report}
    touch ${competitive_mapping_error}
    """
}

process has_enough_reads {
    container 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:1.2.6'

    cpus = 1
    memory = "128MB"

    debug true
    pod label: "name", value: "competitive_mapping_pipeline:has_enough_reads"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    path (json)
    val (threshold)

    output:
    stdout

    script:
    """
    num_reads=\$(jq '.references[] | select(.genome_name == "M.tuberculosis") | .numreads' ${json})

    if  [ \$num_reads -ge ${threshold} ]
    then
        echo "true" | tr -d '\n'
    else
        echo "false" | tr -d '\n'
    fi
    """
}
