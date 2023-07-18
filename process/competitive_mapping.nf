project_dir = projectDir

process competitiveMapping {
    container 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:latest'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: cm_paths
    tuple val(sample_name), path("${sample_name}_cm_1.fastq.gz"), path("${sample_name}_cm_2.fastq.gz"), emit: cm_sample
    path("${sample_name}_competitive_mapping.json"), emit: cm_report

    script:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq.gz"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq.gz"

    """
    cp ${fq1} ${competitive_mapping_file_1}
    cp ${fq2} ${competitive_mapping_file_2}
    cp /app/cov_sorted_h37rv_100k.json ${competitive_mapping_json}
    """

    stub:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq.gz"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq.gz"

    """
    touch ${competitive_mapping_json}
    touch "${competitive_mapping_file_1}"
    touch "${competitive_mapping_file_2}"
    """
}
