#!/usr/bin/env bash
set -e
GISAID_SEQ=$1
GISAID_META=$2
JHU_SEQ=$3
JHU_META=$4
OUT_SEQ=$5
OUT_META=$6

OUT_DIR=$( dirname GISAID_SEQ )
if [[ ! -r "$GISAID_SEQ" ]]
then
	echo "$0: input $GISAID_SEQ not found"
	exit 1
fi
if [[ ! -r "$GISAID_META" ]]
then
        echo "$0: input $GISAID_META not found"
        exit 1
fi
if [[ ! -r "$JHU_SEQ" ]]
then
        echo "$0: input $JHU_SEQ not found"
        exit 1
fi
if [[ ! -r "$JHU_META" ]]
then
        echo "$0: input $JHU_META not found"
        exit 1
fi


echo "Combining Input sequences for nextstrain GISAID"
TP_DIR=${OUT_DIR}/tmp
mkdir -p ${TP_DIR}

cat ${JHU_SEQ} ${GISAID_SEQ} > ${OUT_SEQ}
head -1 ${GISAID_META} > ${OUT_META} && tail -n +2 -q ${JHU_META} ${GISAID_META} >> ${OUT_META}
exit 0
