if [ ${workflow.profile} == 'kubernetes' ]
then
    echo "Running with kubernetes"
    /bin/bash ${projectDir}/lib/s3fs_setup.sh $WORKSPACE
    trap 'PROCESS_EXIT=\$?; /bin/bash $projectDir/lib/s3fs_teardown.sh; exit \$PROCESS_EXIT;' EXIT
fi

set +e

touch ${competitive_mapping_error}

# Perform competitive mapping
# note that sr preset excludes secondary by default, so would need to override to keep secondary alignments
if [ $seq_platform == 'ont' ]
then
    minimap2 -ax map-ont -t $task.cpus --secondary yes -N 1000 ${manifest} ${fqs} > alignments.sam 2>>${competitive_mapping_error}
elif [ $seq_platform == 'illumina' ]
then
    minimap2 -ax sr -t $task.cpus --secondary yes -N 1000 ${manifest} ${fqs[0]} ${fqs[1]} > alignments.sam 2>>${competitive_mapping_error}
fi

# Sort competitive mapping output
samtools sort -@ $task.cpus -o sorted_alignments.bam alignments.sam 2>>${competitive_mapping_error}

# Generate a TSV file summarising coverage against each reference genome
# can use --excl-flags 1540 to include secondary alignments
samtools coverage sorted_alignments.bam > ${cov} 2>>${competitive_mapping_error}

# Index alignments
samtools index sorted_alignments.bam index.bai 2>>${competitive_mapping_error}

# Extract reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome
samtools view -X sorted_alignments.bam index.bai ${h37rv_rname} -o h37rv.bam 2>>${competitive_mapping_error}

# Sort H37Rv mapped reads
samtools sort -n h37rv.bam -o h37rv_sorted.bam 2>>${competitive_mapping_error}

# Extract unmapped reads
samtools view -X sorted_alignments.bam index.bai "*" -o unmapped.bam 2>>${competitive_mapping_error}

# Sort unmapped mapped reads
samtools sort -n unmapped.bam -o unmapped_sorted.bam 2>>${competitive_mapping_error}

# Merge reads aligned to AL123456.3 Mycobacterium tuberculosis H37Rv complete genome and unmapped reads
samtools merge -o merged.bam h37rv_sorted.bam unmapped_sorted.bam

# Convert BAM output to FASTQ
# Default excl-flag is 0x900 which is secondary (0x100) and supplementary (0x800)
# So we need to exclude secondary alignments (0x100) only
if [ $seq_platform == 'ont' ]
then
    samtools fastq --excl-flags 0x100 -@ $task.cpus -0 ${competitive_mapping_file} merged.bam
elif [ $seq_platform == 'illumina' ]
then
    samtools fastq --excl-flags 0x100 -@ $task.cpus -1 ${competitive_mapping_file_1} -2 ${competitive_mapping_file_2}  -0 /dev/null -s /dev/null merged.bam
fi

# Generate alignment summary csv
process_aln_stats --bam sorted_alignments.bam --species_list ${species_list} \
    --output aln_summary.csv \
    

# Generate competitive mapping json
process_mapping --coverage ${cov} --species_list ${species_list} \
    --aln_summary aln_summary.csv --output ${competitive_mapping_report}
