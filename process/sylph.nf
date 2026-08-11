process sylph {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    maxRetries 5

    container {
        params.test_container_cm == "" ? params.container_prefix + '/gpas/competitivemapping_pipeline:4.0.1' : params.test_container_cm
    }
    cpus 4
    memory { 14.GB * task.attempt }

    pod label: "name", value: "gatekeeper_pipeline:sylph"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fqs)
    // Use fixed name to avoid name conflicts.
    // Can use multiple databases
    path "database??.syldb"
    path "taxonomies?/*"
    val seq_platform
    val query_ani_threshold
    val profile_ani_threshold
    // normally 90 for query and 95 for profile

    output:
    tuple val(sample_name), path(sylph_report), emit: sylph_report
    tuple val(sample_name), path(sylph_query), emit: sylph_query
    tuple val(sample_name), path(taxonomy_report), emit: taxonomy_report

    script:
    sylph_query = "sylph_query.tsv"
    sylph_report = "sylph_report.tsv"
    taxonomy_report = "taxonomy_report.tsv"
    """
    if [ ${seq_platform} == 'ont' ]
    then
        sylph sketch -t ${task.cpus} -c 100 -r ${fqs}
    elif [ ${seq_platform} == 'illumina' ]
    then
        sylph sketch -t ${task.cpus} -c 100 \
            -1 ${fqs[0]} -2 ${fqs[1]}
    fi
    # either way a file ending .sylsp is created

    sylph query -t ${task.cpus} -m ${query_ani_threshold} \
        -o ${sylph_query} database*.syldb *.sylsp

    sylph profile -u -t ${task.cpus} -m ${profile_ani_threshold} \
        -o ${sylph_report} database*.syldb *.sylsp


    if [[ \$(wc -l <${sylph_report}) -le 1 ]]
    then
        echo "Sylph failed to generate a profile report"
        echo "Generating a dummy taxonomy report"
        if [ ${seq_platform} == 'ont' ]
        then
            echo -e "#SampleID\t${fqs}" > ${taxonomy_report}
        elif [ ${seq_platform} == 'illumina' ]
        then
            echo -e "#SampleID\t${fqs[0]}" > ${taxonomy_report}
        fi
        echo -e "clade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)" >> ${taxonomy_report}
    else
        sylph-tax taxprof ${sylph_report} -t taxonomies*/* -o taxo_
        cat taxo_* > ${taxonomy_report}
    fi
    """
}
