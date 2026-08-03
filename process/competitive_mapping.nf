process competitive_mapping {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:ccd624c' : params.test_container_cm
    }

    cpus 4
    memory { 12.GB + (4.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:competitive_mapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs), path(manifest), path(species_list)
    val seq_platform
    val ref_for_fastq

    output:
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: ref_reads, optional: true
    tuple val(sample_name), path(competitive_mapping_report), emit: report_json
    tuple val(sample_name), path(competitive_mapping_csv), emit: report_csv

    script:
    ref_reads = "reads_for_assembly.fastq.gz"
    ref_reads_1 = "reads_for_assembly_1.fastq.gz"
    ref_reads_2 = "reads_for_assembly_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
    ref_for_fastq_arg = ref_for_fastq == "" ? "" : "--ref_for_fastq " + ref_for_fastq
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
        ${ref_for_fastq_arg} \
        --cpus ${task.cpus} \
        --output_root "out."

    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}

    # Rename filtered fastqs if we're filtering reads
    if [ ${ref_for_fastq} != '' ]
    then
        if [ ${seq_platform} == 'ont' ]
        then
            mv out.reads.fastq.gz ${ref_reads}
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

process build_manifest {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:ccd624c' : params.test_container_cm
    }

    cpus 4
    memory { 8.GB * task.attempt }

    pod label: "name", value: "competitive_mapping_pipeline:build_manifest"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(sylph_report)
    path ref_genome_dirs
    // Pattern used to avoid name conflicts
    path "taxonomy?/*"
    // optional comma-separated list of accessions to always include as references
    val fixed_refs

    output:
    tuple val(sample_name), path("out.manifest.fasta.gz"), emit: manifest
    tuple val(sample_name), path("out.contigs.csv"), emit: contigs

    script:
    fixed_refs_arg = fixed_refs ? "--fixed_refs " + fixed_refs : ""
    """
    manifest_builder --sylph_report ${sylph_report} \
        ${fixed_refs_arg} \
        --genome_dirs ${ref_genome_dirs} \
        --taxonomy_files taxonomy*/* \
        --cpus ${task.cpus} \
        --output_root "out."
    """
}

process has_enough_reads {
    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:ccd624c' : params.test_container_cm
    }

    cpus 1
    memory "128MB"

    debug true
    pod label: "name", value: "competitive_mapping_pipeline:has_enough_reads"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(json)
    val ref_name
    val threshold

    output:
    tuple val(sample_name), stdout

    script:
    ref_name_arg = ref_name == "" ? "--genome_name no_ref" : "--genome_name " + ref_name
    """
    check_read_count --json_file_path ${json} --read_threshold ${threshold} ${ref_name_arg}
    """
}
