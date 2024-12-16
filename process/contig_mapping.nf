process contigMapping {
    container {
        params.test_container_cm == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/competitivemapping_pipeline:d676ce4' : params.test_container_cm
    }

    cpus 4
    memory {
        params.testing == "" ? ( 4.GB + (8.GB * (task.attempt - 1))) : (4.GB)
    }


    pod label: "name", value: "competitive_mapping_pipeline:contigMapping"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(contigs), path(manifest), path(species_list)

    output:
    tuple val(sample_name), path("contig_species_comparison.csv"), emit: species_comparison
    tuple val(sample_name), path("contig_blast_mapping.tsv"), emit: blast_mapping

    script:
    """
    mkdir outputs
    contig_mapping --contigs ${contigs} --manifest ${manifest} --manifest_contigs ${species_list} \
        --cpus ${task.cpus} --output_root "contig_"
    """
}
