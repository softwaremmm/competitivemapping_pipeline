process competitiveMapping {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:bca3398' : params.test_container_cm
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
    tuple val(sample_name), path("reads_for_assembly*fastq.gz"), emit: cm_tb_reads
    tuple val(sample_name), path(competitive_mapping_report), emit: cm_report
    tuple val(sample_name), path(competitive_mapping_csv), emit: cm_csv

    script:
    tb_reads = "reads_for_assembly.fastq.gz"
    tb_reads_1 = "reads_for_assembly_1.fastq.gz"
    tb_reads_2 = "reads_for_assembly_2.fastq.gz"
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
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
        --ref_for_fastq ${reference_name} \
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

    # clean up large intermediate files
    rm aln.bam
    """
}

process dynamicCompetitiveMapping {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:bca3398' : params.test_container_cm
    }

    cpus 4
    memory { 8.GB + (12.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:dynamicCompetitiveMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs), path(sylph_report)
    path gtdb_genomes_dir
    path assembly_metadata
    // Used to go from assembly to species
    val seq_platform
    val include_whole_genus

    output:
    tuple val(sample_name), path(competitive_mapping_report), emit: cm_report
    tuple val(sample_name), path(competitive_mapping_csv), emit: cm_csv

    script:
    competitive_mapping_report = "species_comparison_report.json"
    competitive_mapping_csv = "species_comparison.csv"
    whole_genera_arg = include_whole_genus ? "--include_whole_genus" : ""
    """
    manifest_builder --sylph_report ${sylph_report} \
        --genome_dirs ${gtdb_genomes_dir} \
        --metadata_files ${assembly_metadata} \
        ${whole_genera_arg} \
        --cpus ${task.cpus} \
        --output_root "out."


    manifest_mapper \
        --seq_platform ${seq_platform} \
        --manifest out.manifest.fasta.gz \
        --reads ${fqs} \
        --cpus ${task.cpus} \
        -o aln.bam

    competitive_mapping \
        --input_bam aln.bam \
        --seq_platform ${seq_platform} \
        --contigs out.contigs.csv \
        --cpus ${task.cpus} \
        --output_root "out."

    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}

    # clean up large intermediate files
    rm aln.bam
    rm out.manifest.fasta.gz
    """
}

process has_enough_reads {
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:bca3398' : params.test_container_cm
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
