process contigMapping {
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:d676ce4' : params.test_container_cm
    }

    cpus 4
    memory {
        params.testing == "" ? ( 8.GB + (8.GB * (task.attempt - 1))) : (8.GB)
    }


    pod label: "name", value: "competitive_mapping_pipeline:contigMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(contigs), path(contig_stats), path(manifest), path(species_list)

    output:
    tuple val(sample_name), path("contig_species_comparison.json"), emit: cm_report
    tuple val(sample_name), path("contig_species_comparison.csv"), emit: cm_csv

    script:
    """
    competitive_mapping --manifest ${manifest} --manifest_contigs ${species_list} \
        --query ${contigs} --platform fasta --query_contig_stats ${contig_stats} \
        --cpus ${task.cpus} --output_root "contig_"
    """
}


process dynamicContigMapping {
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:4627e5b' : params.test_container_cm
    }

    cpus 4
    memory { 8.GB + (12.GB * task.attempt) }

    pod label: "name", value: "competitive_mapping_pipeline:dynamicContigMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(sylph_report), path(contigs), path(contig_stats)
    path gtdb_genomes_dir
    path assembly_metadata
    // Used to go from assembly to species
    val seq_platform
    val include_whole_genus

    output:
    tuple val(sample_name), path(competitive_mapping_report), emit: cm_report
    tuple val(sample_name), path(competitive_mapping_csv), emit: cm_csv
    tuple val(sample_name), path("manifest.fasta.gz"), path("contigs.csv"), emit: manifest

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


    competitive_mapping --platform fasta \
        --query ${contigs} \
        --query_contig_stats ${contig_stats} \
        --manifest out.manifest.fasta.gz \
        --manifest_contigs out.contigs.csv \
        --cpus ${task.cpus} \
        --output_root "out."

    mv out.species_comparison.json ${competitive_mapping_report}
    mv out.species_comparison.csv ${competitive_mapping_csv}

    mv out.manifest.fasta.gz manifest.fasta.gz
    mv out.contigs.csv contigs.csv
    """
}
