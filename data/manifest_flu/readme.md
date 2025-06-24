# Experimental influenza virus manifest and code

* [all_influenza_A.fasta] was pulled from all influenza A sequences in RefSeq https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320
* [build_manifest.py] was written to parse the above FASTA file to figure out which segments refer to Hemagglutinin and Neuraminidase, asigning H1/H2 etc based on the HxNy which is in the header
* [manifest_20250620.csv] is the manifest generated 20/06/2025
