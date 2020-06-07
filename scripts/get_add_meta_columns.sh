#!/usr/bin/env bash 

tsv=$1
next_cols=$2
OUT=$3

awk -F"\t" -v x=${next_cols} 'NR==1 { for ( i=x+1 ; i <= NF ; i++) { print $i }}' ${tsv} > ${OUT}
