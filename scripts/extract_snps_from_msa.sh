#!/usr/bin/env bash

# Script to make vcf from MSA to make filtered vcf and report

set -e

ALN=$1
REF=$2
EXCLUDE_SAMPLES=$3
OUT_VCF=$4
COV=$5

EXE="snp-sites"
SCPATH="/home/idies/workspace/Storage/sramakr4/persistent/nextstrain_scripts/ncov/scripts"
clades_tsv="${HOME}/workspace/covid19/nextstrain/clades.tsv"

if [[ ! -r "$ALN" ]]
then
	echo "$0: input $ALN not found"
	exit 1
fi

if [[ ! -r "$REF" ]]
then
        echo "$0: input $REF not found"
        exit 1
fi

if [[ ! -r "$EXCLUDE_SAMPLES" ]]
then
	echo "$0: input $EXCLUDE_SAMPLES not found"
	exit 1
fi

if [[ ! -r "$clades_tsv" ]]
then
        echo "$0: input $clades_tsv not found in ${HOME}/workspace/covid19/nextstrain/"
        exit 1
fi

if ! [ -x "$(command -v $EXE )" ]; then
  echo "Error: ${EXE} is not installed." >&2
  exit 1
fi

if ! [ -x "$(command -v ${SCPATH}/filter_SNPs.py)" ]; then
  echo "Error: ${SCPATH}/filter_SNPS.py is not present." >&2
  exit 1
fi

if [[ -z "$COV" ]]
then
	echo "Using default cumulative ATCG allele frequency to be 0.5"
	COV=0.5
fi

OUT_DIR=$( dirname $OUT_VCF)
BASE=$(basename $OUT_VCF ".vcf")
PREFIX="snps_msa"
OUT_VCF_FILTERED=${OUT_DIR}/${BASE}_filtered.vcf

mkdir -p ${OUT_DIR}
SNP_VCF=${OUT_DIR}/${PREFIX}.vcf
SNP_VCF_ALL=${OUT_DIR}/${PREFIX}_all.vcf

# Run snp-sites on MSA
echo "INFO : Identifying SNPs from MSA" 
snp-sites -rvp -o ${OUT_DIR}/snps_msa_all ${ALN}
snp-sites -rvpc -o ${OUT_DIR}/snps_msa ${ALN}


if [[ ! -r "$SNP_VCF" ]]
then
        echo "Error running snp-sites ${SNP_VCF}"
        exit 1
fi

# Make formatted VCF 
# Print vcf header

echo "INFO : Formatting VCF from snp-sites"
${SCPATH}/filter_SNPs.py -in ${SNP_VCF} -r ${REF} -l ${EXCLUDE_SAMPLES} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF_FILTERED}
${SCPATH}/filter_SNPs.py -in ${SNP_VCF_ALL} -r ${REF} -l ${EXCLUDE_SAMPLES} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF}

if [[ ! -r "$OUT_VCF" ]]
then
        echo "Error filtering VCF from snp-sites ${OUT_VCF}"
        exit 1
fi

exit 0
