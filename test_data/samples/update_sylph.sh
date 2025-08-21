
db=/home/ubuntu/pipelines/competitivemapping_pipeline/test_data/virus_manifest/part_1/viral_db1_c100.syldb
db2=/home/ubuntu/pipelines/competitivemapping_pipeline/test_data/virus_manifest/part_2/viral_db2_c100.syldb
for sample in viral_mix; do
    sylph query -c 100 $db $db2 ont/$sample/*gz -o ont/$sample/sylph_query.tsv
    sylph profile -u -c 100 $db $db2 ont/$sample/*gz -o ont/$sample/sylph_profile.tsv

    sylph query -c 100 $db $db2 -1 illumina/$sample/*_1.fastq.gz -2 illumina/$sample/*_2.fastq.gz -o illumina/$sample/sylph_query.tsv
    sylph profile -u -c 100 $db $db2 -1 illumina/$sample/*_1.fastq.gz -2 illumina/$sample/*_2.fastq.gz -o illumina/$sample/sylph_profile.tsv
done