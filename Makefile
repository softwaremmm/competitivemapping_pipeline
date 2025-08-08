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

container:
	docker build -t test_container_cm .

test:
	cd tie_break && cargo test
	pytest tests
	nf-test test tests/nextflow/*.test

test_local:
	cd tie_break && cargo test
	pytest tests
	docker build -t test_container_cm .
	nf-test test tests/nextflow/*.test --profile local_docker

clippy:
	cd tie_break && \
	cargo clippy --all --all-features --tests --fix --allow-dirty -- -D warnings
