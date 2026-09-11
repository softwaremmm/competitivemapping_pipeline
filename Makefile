clean:
	find . -type d -name .nextflow | xargs rm -rf
	find . -type d -name work | xargs rm -rf
	find . -type f -regex '.*\.nextflow\.log.*' | xargs rm -f
	find . -type d -name .nf-test | xargs rm -rf

run_myco:
	nextflow run . \
		--seq_platform illumina \
		--input_dir test_data/samples/illumina/chloro \
		--manifest test_data/myco_manifest/manifest.fasta.gz \
		--species_list test_data/myco_manifest/contigs.csv \
		--publish_dir results \
		--refs_for_fastqs "M.tuberculosis,M.pyrenivorans" \
		-profile local_docker \
		-resume

run_flu:
	nextflow run . \
		--seq_platform illumina \
		--input_dir test_data/samples/illumina/flu \
		--manifest test_data/flu_manifest/manifest.fasta.gz \
		--species_list test_data/flu_manifest/contigs.csv \
		--publish_dir flu_results \
		-profile local_docker \
		-resume

run_dynamic:
	nextflow run . \
		--workflow dynamic \
		--seq_platform ont \
		--input_dir test_data/samples/ont/viral_mix \
		--ref_genome_dirs test_data/virus_refs_db/part_1  \
		--sylph_dbs test_data/virus_refs_db/part_1/viral_db1_c100.syldb \
		--taxonomy_files test_data/virus_refs_db/part_1/taxonomy.tsv \
		--publish_dir flu_results \
		-profile local_docker \
		-resume

run_dynamic_list:
	nextflow run . \
		--workflow dynamic \
		--seq_platform ont \
		--input_dir test_data/samples/ont/viral_mix \
		--ref_genome_dirs test_data/virus_refs_db/part_1,test_data/virus_refs_db/part_2   \
		--sylph_dbs test_data/virus_refs_db/part_1/viral_db1_c100.syldb,test_data/virus_refs_db/part_2/viral_db2_c100.syldb  \
		--taxonomy_files test_data/virus_refs_db/part_1/taxonomy.tsv,test_data/virus_refs_db/part_2/taxonomy.tsv \
		--publish_dir flu_results \
		-profile local_docker \
		-resume



container:
	docker build -t test_container_cm .

test:
	pytest tests
	nf-test test tests/nextflow/*.test

test_local:
	find tests/ -type d -name test_outputs | xargs rm -rf
	pytest tests
	docker build -t test_container_cm .
	nf-test test tests/nextflow/*.test --profile local_docker
