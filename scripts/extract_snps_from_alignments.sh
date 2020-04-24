#!/usr/bin/env bash

# Script to make vcf from MSA to make filtered vcf and report

set -e

ALN=$1
EXCLUDE_SAMPLES=$2
OUT_VCF=$3
OUT_REPORT=$4

EXE="snp-sites"

if [[ ! -r "$ALN" ]]
then
	echo "$0: input $ALN not found"
	exit 1
fi

if [[ ! -r "$MASK_SAMPLES" ]]
then
	echo "$0: input $ALN not found"
	exit 1
fi

if ! [ -x "$(command -v ${EXE})" ]; then
  echo 'Error: ${EXE} is not installed.' >&2
  exit 1
fi

OUT_DIR=$( dirname $OUT_VCF)
PREFIX="snps_msa"


mkdir -p ${OUT_DIR}
SNP_VCF=${OUT_DIR}/snps_msa.vcf


# Run snp-sites on MSA
echo "INFO : Identifying SNPs from MSA" 
snp-sites -rvp -o ${OUT_DIR}/snps_msa ${ALN}

if [[ ! -r "$SNP_VCF" ]]
then
        echo "Error running snp-sites ${ALN}"
        exit 1
fi

# Make formatted VCF 
# Print vcf header

echo "INFO : Formatting VCF from snp-sites"

awk '/^#/{if (/^#CHROM/ ) { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\tNUM_SAMPLES\tOCCURENCES\t%OCCURENCES" } else { print  }}' ${SNP_VCF} > ${OUT_VCF} 
awk '!/^#/ { snp=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; for (i=10 ; i <= NF ; i++) { sum += $i  } ; if( sum == 0 ){ pcent = 0} else { pcent=sum/(NF-9) * 100  }; print snp"\t"NF"\t"sum"\t"pcent ; sum = 0 }' ${SNP_VCF} >> ${OUT_VCF} 

echo "INFO: Printing SNP report"

# Print SNP report
echo -e "SNP_POS\tNUM_SAMPLES\t%OCCURENCES" > ${OUT_REPORT}
awk '!/^#/{ print $2"\t"$NF}' ${OUT_VCF} >> ${OUT_REPORT}

exit 0
