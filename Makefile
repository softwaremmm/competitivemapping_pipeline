clean:
	find . -type d -name .nextflow | xargs rm -rf
	find . -type d -name work | xargs rm -rf
	find . -type f -regex '.*\.nextflow\.log.*' | xargs rm -f
	find . -type d -name .nf-test | xargs rm -rf

run:
	nextflow run . \
		--seq_platform illumina \
		--input_dir test_data/chloro_10k \
		--manifest data/manifest/manifest_20231001 \
		--species_list test_data/species_list_manifest_20250324.csv \
		--publish_dir results \
		-resume

run_flu:
	nextflow run . \
		--seq_platform illumina \
		--input_dir flu_input/N1_illumina \
		--manifest data/manifest_flu/all_influenza_A_plus_all_h_all_n.fasta \
		--species_list data/manifest_flu/manifest_20250624_groups.csv \
		--publish_dir flu_results \
		-profile local_docker \
		-resume

container:
	docker build -t test_container_cm .

test:
	pytest tests
	nf-test test tests/nextflow/*.test

test_local:
	pytest tests
	docker build -t test_container_cm .
	nf-test test tests/nextflow/*.test --profile local_docker

clippy:
	cd cm_analyzer && \
	cargo clippy --all --all-features --tests --fix --allow-dirty -- -D warnings
