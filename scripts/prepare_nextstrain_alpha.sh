#!/usr/bin/env bash

SEQ_DIR=$1
SUBMISSIONS_MANIFEST=$2
exclude_submitted_files=$3



#TODO
# If file exists in submissions directory ; include submitted sequence
# To do if exists 5-post-filter/*.complete.fasta 
# else include files from 4-draft-consensus/MD*.nanopolish.consensus.fasta
module_4="4-draft-consensus"
module_5="5-post-filter"

s_prefix="MD"
barcode_prefix="_NB"
draft_pat="MD*.nanopolish.consensus.fasta"
complete_pat="MD*.complete.fasta"
#sub_prefix="MD*.complete.updated.fasta"

dated=$( TZ=":US/Eastern" date +%Y%m%d )
allowed_Ns=5000
allowed_length=29903

for i in `cat ${SUBMISSIONS_MANIFEST}`;
do 
  


done

for i in `ls -d ${SEQ_DIR}/[2-9]*/artic-pipeline/` ; 
do 
   if [[ -d ${i}/${module_4} ]] ; then 
      for j in `ls ${i}/${module_4}/draft_pat`; 
      do 
      	  sum_Ns=$( seqtk comp seqtk ${j} | awk '{ print $9}' ); 
          slen=$( seqtk comp seqtk ${j} | awk '{ print $2}' ); 
          if [ ${slen} -eq ${allowed_length} ] && [ ${sum_Ns} < ${allowed_Ns} ]; then
             sed "/^>/ s/${barcode_prefix}[0-9]*/ /g" ${j} >> ${OUTDIR}/${dated}-jhu_sequences_raw.fasta
   fi 
done   


cat ${SEQ_DIR}/artic-pipeline/4-draft-consensus/MD*.nanopolish.consensus.fasta > ${dated}-jhu_sequences_raw.fasta 


grep ">" 20200503-jhu_sequences_raw.fasta  | sed 's/>//g' | awk '{ print $1}' | awk '{ split($1,a,"/") ; split(a[10],b,"_" );   print $1"\t"b[1] }' | sed 's/MDHP-00/MD-HP/2' | sed 's/MDHP-/MD-HP0/2' | awk 'NF>1  { print }' >  rename_headers.txt

# Rename MD-HP006 and MD-HP007 to MD-HP06 and MD-HP07
# Remove any empty sequences
# then run rename_headers.sh
./rename_headers.sh rename_headers.txt 20200503-jhu_sequences_raw.fasta 20200503-jhu_sequences_renamed.fasta

# Get submitted list from beta_folder
grep ">" ${exclude_submitted_files} | sed 's/>//g' > exclude_beta_2020_05_03.txt
grep ">" 20200503-jhu_sequences_renamed.fasta | sed 's/>//g' > 20200503_seq.txt

awk 'NR==FNR { a[$1] += 1 ; next } !( $1 in a ) { print }' exclude_beta_2020_05_03.txt 20200503_seq.txt > 20200503_filt_seq.txt

seqtk subseq 20200503-jhu_sequences_renamed.fasta 20200503_filt_seq.txt > 20200503-jhu_sequences.fasta
