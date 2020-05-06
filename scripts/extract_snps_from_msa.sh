#!/usr/bin/env bash

# Script to make vcf from MSA to make filtered vcf and report

set -e

ALN=$1
REF=$2
COV=$3
OUT_VCF=$4
EXCLUDE_SAMPLES=$5

EXE="snp-sites"
SCPATH="/home/idies/workspace/Storage/sramakr4/persistent/nextstrain_scripts/ncov/scripts"
clades_tsv="${HOME}/workspace/covid19/nextstrain/clades.tsv"
COUNT_Ns=1000

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

#if [[ ! -r "$EXCLUDE_SAMPLES" ]]
#then
#	echo "$0: input $EXCLUDE_SAMPLES not found"
#	exit 1
#fi

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
ALN_BASE=$(basename $ALN ".fasta")
PREFIX="snps_msa"
OUT_VCF_FILTERED=${OUT_DIR}/${BASE}_filtered.vcf

mkdir -p ${OUT_DIR}

# seqtk comp output columns #4 is N
# chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts

# Make filtered fasta
seqtk comp $ALN | awk -v x=$COUNT_Ns '$9 < x { print $1}' > ${OUT_DIR}/filtered_samples.txt
seqtk subseq $ALN ${OUT_DIR}/filtered_samples.txt > ${OUT_DIR}/${ALN_BASE}_filtered.fasta

# Run SnpEff
#SNP_VCF=${OUT_DIR}/${PREFIX}.vcf
SNP_VCF_ALL=${OUT_DIR}/${PREFIX}_all.vcf

# Run snp-sites on MSA
echo "INFO : Identifying SNPs from MSA" 
#snp-sites -rvp -o ${OUT_DIR}/snps_msa_all ${ALN}
#snp-sites -rvpc -o ${OUT_DIR}/snps_msa ${ALN}

snp-sites -rvp -o ${OUT_DIR}/snps_msa_all ${OUT_DIR}/${ALN_BASE}_filtered.fasta
#snp-sites -rvpc -o ${OUT_DIR}/snps_msa ${OUT_DIR}/${ALN_BASE}_filtered.fasta

if [[ ! -r "$SNP_VCF_ALL" ]]
then
        echo "Error running snp-sites ${SNP_VCF_ALL}"
        exit 1
fi

# Make formatted VCF 
# Print vcf header

echo "INFO : Formatting VCF from snp-sites"

if [[ -r "$EXCLUDE_SAMPLES" ]]
then
   echo "INFO : Filtering samples from ${EXCLUDE_SAMPLES}"
   # ${SCPATH}/filter_SNPs.py -in ${SNP_VCF} -r ${REF} -l ${EXCLUDE_SAMPLES} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF_FILTERED}
   ${SCPATH}/filter_SNPs.py -in ${SNP_VCF_ALL} -r ${REF} -l ${EXCLUDE_SAMPLES} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF}
else
   #${SCPATH}/filter_SNPs.py -in ${SNP_VCF} -r ${REF} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF_FILTERED}
   ${SCPATH}/filter_SNPs.py -in ${SNP_VCF_ALL} -r ${REF} -m ${clades_tsv} -c ${COV} -o ${OUT_VCF}
fi 

if [[ ! -r "$OUT_VCF" ]]
then
        echo "Error filtering VCF from snp-sites ${OUT_VCF}"
        exit 1
fi

exit 0
