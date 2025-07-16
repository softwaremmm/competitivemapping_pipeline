# Competitive Mapping analyzer
Rust code to do compare all the alignments for a read and select the best alignment, with the result being a single alignment per read.

### Running
Can be run directly with:

```bash
cargo run --release -- \
    --input-bam aln.bam \
    --contigs contigs.csv \
    --threads 10 \
    --parameters params.yml \
    --output-root results
```

Note that `aln.bam` must be sorted by name:
```bash
samtools sort -n input.bam -o name_sorted.bam
```

use `--debug` to print all alignment info to files.

### Testing

Use
```bash
cargo test
```

If you've made changes and want to update test expectations can use:
```bash
UPDATE_EXPECTATIONS=1 cargo test
```

## Design
Aims:
1. Want to assign majority of reads uniquely to a reference, to allow for subsequent analysis
2. Can only disintinguish species if sufficient difference between reference genomes. This can either be by ANI distance or by genomes having enough unique sections of DNA.
3. Medics want NCBI species for familiarity.

For point 2. we group targets which are close enough by ANI and allow reads to map to references in both. Soem future step could disentangle these "complexes".

WARNING: if some genes are similar to both A and B, but missing in our reference for A then will appear unique to B.
Want complete references to avoid this. Or even allow multiple references to account for this.
For example, have found long reads which megablast online correctly identifies it as michigensis, but isn't in aligned to the RefSeq reference!


### Proceed in two rounds
Proceed in two rounds
1. Selecting high quality signals to determine species present and removing dubious ones
2. reassigning reads uniquely to the likely refs, with winner takes all for close ties

Parameter yml files are used for specifying the actual thresholds here.

Alignments are filtered based on:
- Query coverage
- sequence divergence
- Alignment score (which is a composite of the previous two)

Can use the read quality scores to determine the expected sequence divergence.

In general:
- Round 1 looks for 95% nucleotide identity (ONT adjusted based on read error rate) with majority of query covered
- Round 2 only requires 80% nucleotide identity, but ignores alignments 5% under best. It is very relaxed with query coverage.

---
---


# Misc

## Map quality
Note that the mapping quality column can be unhelpful as is set to 0 when reads maps **equally** well to multiple places.
the AS tag (only in the sam) gives the Smith-Waterman alignment score.

Mapping quality score is based on the ratio between the best and second best alignment.
This makes it less helpful for competitive mapping scenarios.

## Expected coverage from possion

If true read mean depth is d can model depth as
```
Y ~ poisson(d)
P(Y = 0) = exp(-d)
so
expected cov = 1 - exp(-d)

or in reverse to estimate depth from coverage
d' = ln(1/ 1-cov)
```

Have copied Sylph's approach for trying to estimate mean depth when only considering unique alignments.

TODO: This needs fixing up. Doesn't seem to work particularly well at the moment. Likely due to not having the independence of K-mers.

## Notes on Minimap2 Alignment
Purely interested in how alignment, chaining and mapping scores relate to each other.

- Chaining is how minimap2 initially finds mapping
- then does basepair alignment and produces AS as well as a max segment score (ms).
- AS is used to determine the primary mapping. For paired reads I think it uses the combination of the two, so one half may actually have a higher score elsewhere.
- For alignments often both parts of pair will map to location, but sometimes only get one half mapping (for both primary and secondary)
- The tag "read mapped in proper pair" or "is_properly_segmented" means that the primary alignment for both parts is to same location.

quirk with supplementary
If a read has two parts so gives rise to a supplementary alignment. The alns which are secondary to the supplementary one are only marked with secondary flag.
As such it is not easy to tell if a secondary alnmnent overlaps with the primary or not.

If AS is very similar then primary may actually be slightly less than max as it is considered a multimap and so primary is random.


## Pileups
I did some testing and when running a pileup to count read depths you find that:
- mismatches/skips/deletions still count a coverage
- samtools flags like supplementary/secondary still counted to read depth

This has been replicated in this code.
