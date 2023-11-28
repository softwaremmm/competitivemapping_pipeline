if [ ${workflow.profile} == 'kubernetes' ]
then
    echo "Running with kubernetes"
    /bin/bash ${projectDir}/lib/s3fs_setup.sh $WORKSPACE
fi

set +e

touch ${competitive_mapping_error}

# Perform competitive mapping
if [ $seq_platform == 'ont' ]
then
    minimap2 -ax map-ont -t $task.cpus ${manifest} ${fqs} > ${sample_name}_alignments.sam 2>>${competitive_mapping_error}
elif [ $seq_platform == 'illumina' ]
then
    minimap2 -ax sr -t $task.cpus ${manifest} ${fqs[0]} ${fqs[1]} > ${sample_name}_alignments.sam 2>>${competitive_mapping_error}
fi

# Sort competitive mapping output
samtools sort -@ $task.cpus -o ${sample_name}_sorted_alignments.bam ${sample_name}_alignments.sam 2>>${competitive_mapping_error}

# Generate a TSV file summarising coverage against each reference genome
samtools coverage ${sample_name}_sorted_alignments.bam > ${cov} 2>>${competitive_mapping_error}

# Index alignments
samtools index ${sample_name}_sorted_alignments.bam index.bai 2>>${competitive_mapping_error}

# Extract reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome
samtools view -X ${sample_name}_sorted_alignments.bam index.bai ${h37rv_rname} -o h37rv.bam 2>>${competitive_mapping_error}

# Sort H37Rv mapped reads
samtools sort -n h37rv.bam -o h37rv_sorted.bam 2>>${competitive_mapping_error}

# Extract unmapped reads
samtools view -X ${sample_name}_sorted_alignments.bam index.bai "*" -o unmapped.bam 2>>${competitive_mapping_error}

# Sort unmapped mapped reads
samtools sort -n unmapped.bam -o unmapped_sorted.bam 2>>${competitive_mapping_error}

# Merge reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome and unmapped reads
samtools merge -o merged.bam h37rv_sorted.bam unmapped_sorted.bam

# Convert BAM output to FASTQ
if [ $seq_platform == 'ont' ]
then
    samtools fastq -@ $task.cpus -0 ${competitive_mapping_file} merged.bam
elif [ $seq_platform == 'illumina' ]
then
    samtools fastq -@ $task.cpus -1 ${competitive_mapping_file_1} -2 ${competitive_mapping_file_2}  -0 /dev/null -s /dev/null merged.bam
fi

# Generate competitive mapping json
process_mapping --coverage ${cov} --species_list ${species_list} --output ${competitive_mapping_report}

if [ ${workflow.profile} == 'kubernetes' ]
then
    /bin/bash ${projectDir}/lib/s3fs_teardown.sh
fi