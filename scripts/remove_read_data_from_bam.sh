
# For test data sometimes only need the alignment information
# This script removes all the fastq data from a BAM file
# leaving a much smaller file

# get input bam file
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_bam_file> <output_bam_file>"
    exit 1
fi
input_bam_file=$1
output_bam_file=$2
if [ ! -f "$input_bam_file" ]; then
    echo "Error: Input BAM file '$input_bam_file' does not exist."
    exit 1
fi


echo "Removing read data from BAM file: $input_bam_file"

samtools view -h $input_bam_file | \
    awk 'BEGIN {OFS="\t"} /^@/ {print; next} { $10="*"; $11="*"; print }' | \
    samtools view -o $output_bam_file
