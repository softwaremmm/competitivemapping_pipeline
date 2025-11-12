process competitiveMapping {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:f363973' : params.test_container_cm
    }

    cpus 4
    memory { 12.GB + (4.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:competitiveMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs)
    path manifest
    path species_list
    val seq_platform
    val reference_name

    output:
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: ref_reads, optional: true
    tuple val(sample_name), path(competitive_mapping_report), emit: report_json
    tuple val(sample_name), path(competitive_mapping_csv), emit: report_csv

    script:
    tb_reads = "reads_for_assembly.fastq.gz"
    ref_reads_1 = "reads_for_assembly_1.fastq.gz"
    ref_reads_2 = "reads_for_assembly_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
    ref_for_fastq = reference_name == "" ? "" : "--ref_for_fastq " + reference_name
    """
    manifest_mapper \
        --seq_platform ${seq_platform} \
        --manifest ${manifest} \
        --reads ${fqs} \
        --cpus ${task.cpus} \
        -o aln.bam

    competitive_mapping \
        --input_bam aln.bam \
        --seq_platform ${seq_platform} \
        --contigs ${species_list} \
        ${ref_for_fastq} \
        --cpus ${task.cpus} \
        --output_root "out."

    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}

    # Rename filtered fastqs if we're filtering reads
    if [ ${reference_name} != '' ]
    then
        if [ ${seq_platform} == 'ont' ]
        then
            mv out.reads.fastq.gz ${tb_reads}
        elif [ ${seq_platform} == 'illumina' ]
        then
            mv out.reads_1.fastq.gz ${ref_reads_1}
            mv out.reads_2.fastq.gz ${ref_reads_2}
        fi
    fi

    # clean up large intermediate files
    rm aln.bam
    """
}

// WARNING: Experimental process
process tie_break {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:f363973' : params.test_container_cm
    }

    cpus 8
    memory { 8.GB + (8.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:tie_break"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs)
    path manifest
    path species_list
    val seq_platform
    path tie_break_params
    val reference_name

    output:
    tuple val(sample_name), path(tie_break_report), emit: report_csv
    tuple val(sample_name), path(tie_break_stats), emit: stats
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: ref_reads, optional: true

    script:
    tie_break_report = "species_comparison.csv"
    tie_break_stats = "tie_break_stats.yaml"
    ref_reads_root = "reads_for_assembly"
    """
    manifest_mapper \
        --seq_platform ${seq_platform} \
        --manifest ${manifest} \
        --reads ${fqs} \
        --cpus ${task.cpus} \
        --sort_by_name \
        -o aln.bam

    tie_break \
        --input-bam aln.bam \
        --contigs ${species_list} \
        --threads ${task.cpus} \
        --parameters ${tie_break_params} \
        --output-root "out."

    mv out.alignment_summary.csv ${tie_break_report}
    mv out.stats.yaml ${tie_break_stats}

    # If reference_name is provided, filter reads
    if [ "${reference_name}" != "" ]
    then
        extract_reads -f ${fqs} -a out.alns_round_2.csv \
            -c out.references.csv \
            -r ${reference_name} \
            -o ${ref_reads_root}

        gzip ${ref_reads_root}*
    fi

    # clean up large intermediate files if present
    find . -type f -name "aln.bam" -delete
    find . -type f -name "out.alns_round_*.csv" -delete
    """
}

// WARNING: Experimental process
process dynamic_tie_break {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:f363973' : params.test_container_cm
    }

    cpus 4
    memory { 8.GB + (12.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:dynamic_tie_break"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs), path(sylph_report)
    path genome_dirs
    // Pattern used to avoid name conflicts
    path "metadata?/*"
    val seq_platform
    val include_whole_genus
    path tie_break_params

    output:
    tuple val(sample_name), path(tie_break_report), emit: report_csv
    tuple val(sample_name), path(tie_break_stats), emit: stats

    script:
    whole_genera_arg = include_whole_genus ? "--include_whole_genus" : ""
    tie_break_report = "species_comparison.csv"
    tie_break_stats = "tie_break_stats.yaml"
    """
    manifest_builder --sylph_report ${sylph_report} \
        --genome_dirs ${genome_dirs} \
        --metadata_files metadata*/* \
        ${whole_genera_arg} \
        --cpus ${task.cpus} \
        --output_root "out."

    manifest_mapper \
        --seq_platform ${seq_platform} \
        --manifest out.manifest.fasta.gz \
        --reads ${fqs} \
        --cpus ${task.cpus} \
        --sort_by_name \
        -o aln.bam

    tie_break \
        --input-bam aln.bam \
        --contigs out.contigs.csv \
        --threads ${task.cpus} \
        --parameters ${tie_break_params} \
        --output-root "out."

    mv out.alignment_summary.csv ${tie_break_report}
    mv out.stats.yaml ${tie_break_stats}

    # clean up large intermediate files if present
    find . -type f -name "aln.bam" -delete
    find . -type f -name "out.manifest.fasta.gz" -delete
    """
}

process has_enough_reads {
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:f363973' : params.test_container_cm
    }

    cpus 1
    memory "128MB"

    debug true
    pod label: "name", value: "competitive_mapping_pipeline:has_enough_reads"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(json)
    val threshold

    output:
    tuple val(sample_name), stdout

    script:
    """
    check_read_count --json_file_path ${json} --read_threshold ${threshold}
    """
}
