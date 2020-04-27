#!/usr/bin/env bash

#
echo " INFO : Creating local Snakemake main"
cp Snakefile Snakefile_main
cp config/Snakefile.yaml config/Snakefile_local.yaml

# Add additional config tags for local build
echo "INFO : Adding tags for local build "
echo "local_sequences: <path to your loacl sequences.fasta>" >> config/Snakefile_local.yaml
echo "local_gisaid_metadata: <path to your local metadata tsv in gisaid format>" >> config/Snakefile_local.yaml
### Add additional tags as needed

# Replace config to local config
sed -i 's/config\/Snakefile.yaml/config\/Snakefile_local.yaml/g' Snakefile_main
sed -i '1i include: "Snakefile_prepare"' Snakefile_main 
# Replace download tags
AWS_SEQ='aws s3 cp s3:\/\/nextstrain\-ncov\-private\/sequences.fasta'
AWS_META='aws s3 cp s3:\/\/nextstrain\-ncov\-private\/metadata.tsv'
LOCAL_SEQ='cp "final\/sequences.fasta"'
LOCAL_META='cp "final\/metadata.tsv"'
sed -i -e "s/${AWS_SEQ}/${LOCAL_SEQ}/g" Snakefile_main
sed -i -e "s/${AWS_META}/${LOCAL_META}/g" Snakefile_main

# Setup done
echo "INFO : Setup Completed Successfully"
echo "#######################################################"
echo " Please update local config tags in config/Snakefile_local.yaml "  
echo "#######################################################"

