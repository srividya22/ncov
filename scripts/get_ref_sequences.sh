#!/usr/bin/env bash
set -e
REF_SARSCOV2_IN=$1
REF_META=$2
INCLUDE_FILE=$3
REF_SARSCOV2_OUT=$4
REF_META_OUT=$5

OUT_DIR=$( dirname $REF_SARSCOV2_OUT )
PREFIX=$( basename $REF_SARSCOV2_OUT ".fasta" )

if [[ ! -r "$REF_SARSCOV2_IN" ]]
then
	echo "$0: input $REF_SARSCOV2_IN not found"
	exit 1
fi
if [[ ! -r "$REF_META" ]]
then
        echo "$0: input $REF_META not found"
        exit 1
fi
if [[ ! -r "$INCLUDE_FILE" ]]
then
        echo "$0: input $INCLUDE_FILE not found"
        exit 1
fi 

echo "Getting the reference only seq and metadata from REF file $REF_SARSCOV2_IN"
TP_DIR=${OUT_DIR}/tmp
mkdir -p ${TP_DIR}

seqtk subseq $REF_SARSCOV2_IN $INCLUDE_FILE > $REF_SARSCOV2_OUT

echo "Filtering REF file $REF_SARSCOV2_IN based on metadata"
awk '/^>/ { print } ' $REF_SARSCOV2_OUT | sed 's/>//' > ${TP_DIR}/${PREFIX}_ref_fasta_header.txt
head -1 $REF_META > $REF_META_OUT
awk 'NR==FNR { a[$1] = $0 ; next } ($1 in a ) { OFS="\t" ; print $0 }' ${TP_DIR}/${PREFIX}_ref_fasta_header.txt $REF_META >> $REF_META_OUT
exit 0
