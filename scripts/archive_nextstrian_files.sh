#!/usr/bin/env bash

####################################
# Script to archive nextstrain files
# Make dated folder to archive files used in nexstrain
# Copy output files to latest folder
####################################


output=$1
dest=$2 # Path to the nextstrain archive folder 
latest_dest=$3 # Path to the latest folder for alpha/beta nextstrain
deploy_dir=$4 # Path to deploy the nextstrain files before restart

mkdir -p $dest

if [[ ! -d "$latest_dest" ]]
then
        echo "$0: input ${latest_dest} not found ; Creating"
        mkdir -p ${latest_dest}
fi

if [[ ! -d "$deploy_dir" ]]
then
        echo "$0: input ${deploy_dir} doesnt exists"
        mkdir -p ${deploy_dir}
fi


# Make dated folder inside destination
f_date=$( TZ=":US/Eastern" date +%Y-%m-%d )

OUT="${dest}/${f_date}"
mkdir -p ${OUT}

echo "INFO: Writing files into archive folder ${OUT}"
mv $output/* ${OUT}

# Update latest folder
echo "Updating files in ${latest_dest} directory"
rm -rf ${latest_dest}/*

if [[ ! -r "${OUT}/final/sequences.fasta" ]]
then
        echo "${OUT}/final/sequences.fasta doesnt exists"
        exit 1 
fi

if [[ ! -r "${OUT}/results/masked.fasta" ]]
then
        echo "${OUT}/results/masked.fasta  doesnt exists"
        exit 1
fi

if [[ ! -r "${OUT}/results/alignments.vcf" ]]
then
        echo "${OUT}/results/alignments.vcf doesnt exists"
        exit 1
fi

if [[ ! -d "${OUT}/auspice" ]]
then
        echo "${OUT}/auspice directory doesnt exists"
        exit 1
fi

cp ${OUT}/final/sequences.fasta ${latest_dest}
cp ${OUT}/results/masked*.fasta ${latest_dest}
cp ${OUT}/results/alignments.vcf ${latest_dest}

# Copy auspice files to dest folder
echo "Copying auspice files to ${deploy_dir} directory"
cp ${OUT}/auspice/* ${deploy_dir}

echo "Done"

