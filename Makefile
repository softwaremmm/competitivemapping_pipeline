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

test:
	pytest tests
	nf-test test tests/nextflow/*.test

test_local:
	pytest tests
	docker build -t test_container_cm .
	nf-test test tests/nextflow/*.test --profile local_docker
