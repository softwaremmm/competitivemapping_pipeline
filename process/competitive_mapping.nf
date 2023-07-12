project_dir = projectDir

process competitiveMapping {

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_competitive_mapping.json}'
    publishDir "${params.output_dir}/$sample_name/competitive_mapping_out", mode: 'copy', pattern: '*_cm_{1,2}.fastq', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: cm_paths
    tuple val(sample_name), path("${sample_name}_cm_1.fastq"), path("${sample_name}_cm_2.fastq"), emit: cm_sample
    path ("${sample_name}_competitive_mapping.json"), emit: cm_report

    script:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq"

    """
    cp ${fq1} ${competitive_mapping_file_1}
    cp ${fq2} ${competitive_mapping_file_2}
    cp ${baseDir}/lib/test_data/cov_sorted_h37rv_100k.json ${competitive_mapping_json}
    """

    stub:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq"

    """
    touch ${competitive_mapping_json}
    touch "${competitive_mapping_file_1}"
    touch "${competitive_mapping_file_2}"
    """
}
