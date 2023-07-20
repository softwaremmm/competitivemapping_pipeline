project_dir = projectDir

process competitiveMappingMock {
    container 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:latest'

    input:
    tuple val(sample_name), path(fq1), path(fq2)
    path (manifest)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: cm_paths
    tuple val(sample_name), path("${sample_name}_cm_1.fastq.gz"), path("${sample_name}_cm_2.fastq.gz"), emit: cm_sample
    path("${sample_name}_competitive_mapping.json"), emit: cm_report
    path("manifest_summary.tsv"), emit: cm_manifest_summary

    script:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq.gz"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq.gz"
    manifest_summary = "manifest_summary.tsv"

    """
    bash ${baseDir}/lib/manifest_summary.sh ${manifest} ${manifest_summary}
    cp ${fq1} ${competitive_mapping_file_1}
    cp ${fq2} ${competitive_mapping_file_2}
    cp /app/cov_sorted_h37rv_100k.json ${competitive_mapping_json}
    """

    stub:
    competitive_mapping_json = "${sample_name}_competitive_mapping.json"
    competitive_mapping_file_1 = "${sample_name}_cm_1.fastq.gz"
    competitive_mapping_file_2 = "${sample_name}_cm_2.fastq.gz"
    manifest_summary = "manifest_summary.tsv"

    """
    touch ${competitive_mapping_json}
    touch "${competitive_mapping_file_1}"
    touch "${competitive_mapping_file_2}"
    touch "${manifest_summary}
    """
}


process competitiveMapping{
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    path (manifest)

    output:

    tuple val(sample_name), path("h37rv_1.fastq.gz"), path("h37rv_2.fastq.gz"), emit: cm_sample
    path("competitive_mapping.json"), emit: cm_report


    script:
    competitive_mapping_json = "competitive_mapping.json"
    competitive_mapping_file_1 = "h37rv_1.fastq.gz"
    competitive_mapping_file_2 = "h37rv_2.fastq.gz"
    manifest_summary = "manifest_summary.tsv"
    cov = "cov_${sample_name}.tsv"
    h37rv_rname="AL123456.3"

    """
    bash ${baseDir}/lib/manifest_summary.sh ${manifest} ${manifest_summary}

    minimap2 -ax sr -t12 ${manifest} ${fq1} ${fq2} | samtools sort -@ 2 -o ${sample_name}_sorted_alignments.bam

    samtools coverage ${sample_name}_sorted_alignments.bam > ${cov}

    csvjoin --columns "#rname" ${cov} ${manifest_summary} > joined.csv

    cat joined.csv |
    csvcut --delimiter="," -c "#rname","genome_name","endpos","coverage","numreads","meandepth" |
    csvsort -r -c "coverage" |
    csvjson -i 4 > ${competitive_mapping_json}

    # Index alignments
    samtools index ${sample_name}_sorted_alignments.bam index.bai

    # Extract reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome
    samtools view -X ${sample_name}_sorted_alignments.bam index.bai index.bai ${h37rv_rname} -o h37rv.bam


    # Sort reads
    samtools sort -n h37rv.bam -o h37rv_sorted.bam

    # Convert BAM output to FASTQ
    samtools fastq -@ 2 -1 ${competitive_mapping_file_1} -2 ${competitive_mapping_file_2}  -0 /dev/null -s /dev/null h37rv_sorted.bam

    """

    stub:
    competitive_mapping_json = "competitive_mapping.json"
    competitive_mapping_file_1 = "h37rv_1.fastq.gz"
    competitive_mapping_file_2 = "h37rv_2.fastq.gz"

    """
    touch ${competitive_mapping_json}
    touch "${competitive_mapping_file_1}"
    touch "${competitive_mapping_file_2}"
    """
}
