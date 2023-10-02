project_dir = projectDir

process competitiveMapping{

    container 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:v0.0.1-rc1'

    cpus = 4
    memory = "8GB"

    input:
    tuple val(sample_name), path(fq1), path(fq2)
    path (manifest)

    output:
    tuple val(sample_name), path("h37rv_1.fastq.gz"), path("h37rv_2.fastq.gz"), emit: cm_sample
    path("species_comparison_report.json"), emit: cm_report
    path("species_comparison_error.json"), emit: cm_error

    script:
    competitive_mapping_file_1 = "h37rv_1.fastq.gz"
    competitive_mapping_file_2 = "h37rv_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_error = "species_comparison_error.json"
    manifest_summary = "manifest_summary.tsv"
    cov = "cov_${sample_name}.tsv"
    h37rv_rname="AL123456.3"

    """
    if [ ${workflow.profile} == 'kubernetes' ]
    then
        echo "Running with kubernetes"
        /bin/bash ${projectDir}/lib/s3fs_setup.sh $WORKSPACE
    fi

    set +e

    touch ${competitive_mapping_error}

    # Create manifest summary
    bash ${moduleDir}/../lib/manifest_summary.sh ${manifest} ${manifest_summary}

    # Perform competitive mapping
    minimap2 -ax sr -t12 ${manifest} ${fq1} ${fq2} |

    # Sort competitive mapping output
    samtools sort -@ 2 -o ${sample_name}_sorted_alignments.bam 2>>${competitive_mapping_error}

    # Generate a TSV file summarising coverage against each reference genome
    samtools coverage ${sample_name}_sorted_alignments.bam > ${cov} 2>>${competitive_mapping_error}

    # Index alignments
    samtools index ${sample_name}_sorted_alignments.bam index.bai 2>>${competitive_mapping_error}

    # Extract reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome
    samtools view -X ${sample_name}_sorted_alignments.bam index.bai index.bai ${h37rv_rname} -o h37rv.bam 2>>${competitive_mapping_error}

    # Extract unmapped reads
    samtools view -X ${sample_name}_sorted_alignments.bam index.bai index.bai "*" -o unmapped.bam 2>>${competitive_mapping_error}

    # Merge reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome and unmapped reads
    samtools merge -o merged.bam h37rv.bam unmapped.bam

    # Convert BAM output to FASTQ
    samtools fastq -@ 2 -1 ${competitive_mapping_file_1} -2 ${competitive_mapping_file_2}  -0 /dev/null -s /dev/null merged.bam

    # Generate competitive mapping json
    bash ${moduleDir}/../lib/generate_competitive_mapping_json.sh --cov ${cov}  --manifest-summary ${manifest_summary} --competitive-mapping-json ${competitive_mapping_report}
    """

    stub:
    competitive_mapping_file_1 = "h37rv_1.fastq.gz"
    competitive_mapping_file_2 = "h37rv_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_error = "species_comparison_error.json"

    """
    touch ${competitive_mapping_file_1}
    touch ${competitive_mapping_file_2}
    touch ${competitive_mapping_report}
    touch ${competitive_mapping_error}
    """
}

process has_enough_reads {

    container 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:latest'

    cpus = 1
    memory = "128MB"


    input:
    path (json)
    val (threshold)

    output:
    stdout

    script:

    """
    if [ ${workflow.profile} == 'kubernetes' ]
    then
        /bin/bash ${projectDir}/lib/s3fs_setup.sh $WORKSPACE
    fi

    num_reads=\$(jq '.[] | select(.genome_name == "Mycobacterium tuberculosis H37Rv complete genome") | .numreads' ${json})

    if  [ \$num_reads -ge ${threshold} ]
    then
        echo "true" | tr -d '\n'
    else
        echo "false" | tr -d '\n'
    fi
    """
}

