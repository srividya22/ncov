#!/bin/bash

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi


#Setup files
SUB_MANIFEST="/home/idies/workspace/covid19/jhu_sequences/master/submission_manifest.tsv"
SEQ_PATH="/home/idies/workspace/covid19/sequencing_runs/"
FASTA_PATH="artic-pipeline/4-draft-consensus"
SUBMISSIONS_PATH="/home/idies/workspace/covid19/jhu_sequences/submissions"
SUBMITTED_SEQ=${SUBMISSIONS_PATH}"/all_submitted_jhu_sequences.fasta"
MASTER_META="/home/idies/workspace/covid19/jhu_sequences/master/jhu_sequences_metadata_master_paper_isolates.tsv"
OUTPUT_PATH="/home/idies/workspace/covid19/jhu_sequences/alpha"

#SCRIPTS_PATH="/home/idies/workspace/covid19/code/ncov/pipeline_scripts"
SCRIPTS_PATH="/home/idies/workspace/Temporary/sramakr4/scratch/nextstrain_runs/alpha/ncov/scripts"
#SCRIPTS_PATH="/home/idies/workspace/Temporary/sramakr4/scratch/nextstrain_runs/test/ncov/scripts/prepare_nextstrain_alpha.py"
#source /home/idies/workspace/covid19/bashrc
source /home/idies/workspace/Storage/sramakr4/persistent/bashrc
conda activate nextstrain


${SCRIPTS_PATH}/prepare_nextstrain_alpha.py --submission_manifest ${SUB_MANIFEST} --seq-path ${SEQ_PATH} --fasta-path ${FASTA_PATH} --submitted-fasta ${SUBMITTED_SEQ} --master-meta ${MASTER_META} -out ${OUTPUT_PATH}

echo "DONE"
