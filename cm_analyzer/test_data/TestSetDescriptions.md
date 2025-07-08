# Notes about test data

bam files are large because they contain all the fastq data in them.
But for this test we only really need the alignment infomation.
However, you then lose the read quality score, so we use tb_small to test that aspect.

## tb_small (ONT)
The bam file includes the fastq data so quality score filters are tested as required

## tb_intracellulare (ONT)
Made from a mix of two samples so as to contain reads from both TB and intracellulare.

### with ani grouping
This just changes the contigs file so as to group the intracellulare complex to all be group 77

## Fortuitum (Illumina)
A fortuitum sample from the Tim walker set (ebdf7096-03ed-4c18-893b-abc61b4163fb, SRR32024568).

Key thing is that illumina reads are paired and so making sure that filters are applied to each part of read independently.
