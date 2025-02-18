clean:
	find . -type d -name .nextflow | xargs rm -rf
	find . -type d -name work | xargs rm -rf
	find . -type f -regex '.*\.nextflow\.log.*' | xargs rm -f
	find . -type d -name .nf-test | xargs rm -rf

run:
	nextflow run . \
		--input_dir "/mnt/volume_1/hieu_test" \
		--manifest "/mnt/volume_1/knowledge/manifest/manifest_20231001" \
		--species_list "/mnt/volume_1/knowledge/manifest/species_list_manifest_20240710.csv" \
		--seq_platform "illumina" \
		-resume \
		-profile standard
