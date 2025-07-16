# Experimental influenza virus manifest and code

* [all_influenza_A.fasta] was pulled from all influenza A sequences in RefSeq https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320
* [all_influenza_A_plus_all_h_all_n.fasta] is the file above with additional high quality segments in order to cover some inherent variation in what people call H and N
* [all_influenza_A_plus_all_h_all_n_filtered.fasta] is the above fasta but filtered for duplicates (both by rname and ANI) to ensure it works.
* [ani_duplicates.txt] is a list of the rnames within the above file which are 100% ANI similar to others according to `skani`. These are filtered out during manifest building to reduce complexity without. Generate with `skani triangle -i -E all_influenza_A_plus_all_h_all_n.fasta | grep -E "fasta\s+100.00" | cut -d $'\t' -f 6 | sort | uniq | cut -d "|" -f 1 > ani_duplicates.txt`. If you add segements to the fasta, it is recommended to run `build_manifest.py` without ANI de-duplication first to spit out any duplicate rnames, which should then be removed. Without this, _both_ of the duplicate segments will be removed as they are ANI identical
* [build_manifest.py] was written to parse the above FASTA file to figure out which segments refer to Hemagglutinin and Neuraminidase, asigning H1/H2 etc based on the HxNy which is in the header. Also de-duplicates based on rname and ANI
* [manifest_20250620.csv] is the manifest generated 20/06/2025
* [manifest_20250624.csv] is the manifest generated 24/06/2025 from this
* [manifest_20250704.csv] is the manifest generated 04/07/2025
