#!/usr/bin/env bash
set -e
GISAID_SARSCOV2_IN=$1
GISAID_META=$2
GISAID_SARSCOV2_OUT=$3
GISAID_META_OUT=$4
MIN_LENGTH=$5

OUT_DIR=$( dirname $GISAID_SARSCOV2_OUT )
PREFIX=$( basename $GISAID_SARSCOV2_OUT ".fasta" )

if [[ ! -r "$GISAID_SARSCOV2_IN" ]]
then
	echo "$0: input $GISAID_SARSCOV2_IN not found"
	exit 1
fi
if [[ ! -r "$GISAID_META" ]]
then
        echo "$0: input $GISAID_META not found"
        exit 1
fi
if [[ -z "$MIN_LENGTH" ]]
then
	echo "Using default minimum length of 25000"
	MIN_LENGTH=25000
fi

echo "Normalizing GISAID file $GISAID_SARSCOV2_IN  (min length $MIN_LENGTH)"
# Check

# Remove leading virus name prefix from sequence names
# Remove embedded spaces in sequence names (Hong Kong sequences)
# Remove trailing |EPI_ISL_id|datestamp from sequence names
# Remove sequences shorter than minimum length
# Eliminate duplicate sequences (keep only the first seen)
# Then match columns to metadata
# Remove metadata without seq
# Remove seq without metadata
TP_DIR=${OUT_DIR}/tmp
mkdir -p ${TP_DIR}
#cat $GISAID_SARSCOV2_IN |

	sed 's/^>[hn][Cc]o[Vv]-19\//>/g' $GISAID_SARSCOV2_IN |	# remove leading prefix
	sed 's/ //g' |					# remove embedded spaces
	sed 's/|.*$//' | 				# remove trailing metadata
	awk "BEGIN{RS=\">\";FS=\"\n\"}length>$MIN_LENGTH{print \">\"\$0}" |	# remove short seqs
	awk 'BEGIN{RS=">";FS="\n"}!x[$1]++{print ">"$0}' | 	# remove duplicates	
        grep -v '^>*$' > ${TP_DIR}/${PREFIX}_norm.fasta
        
echo "Filtering GISAID file $GISAID_SARSCOV2_IN based on metadata"
        awk '/^>/ { print } ' ${TP_DIR}/${PREFIX}_norm.fasta | sed 's/>//' > ${TP_DIR}/${PREFIX}_fasta_header.txt
        head -1 $GISAID_META > $GISAID_META_OUT
        awk 'NR==FNR { a[$1] = $0 ; next } ($1 in a ) { OFS="\t" ; print $0 }' ${TP_DIR}/${PREFIX}_fasta_header.txt $GISAID_META >> $GISAID_META_OUT
        awk 'NR==FNR { a[$1] = $0 ; next } ($1 in a ) { OFS="\t" ; print $1 }' $GISAID_META_OUT ${TP_DIR}/${PREFIX}_fasta_header.txt > ${TP_DIR}/${PREFIX}_filter_seq.txt 
        seqtk subseq ${TP_DIR}/${PREFIX}_norm.fasta ${TP_DIR}/${PREFIX}_filter_seq.txt > $GISAID_SARSCOV2_OUT 
        #rm -rf ${TP_DIR}
exit 0
