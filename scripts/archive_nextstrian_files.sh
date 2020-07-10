#!/usr/bin/env bash

output=$1
dest=$2

cp -rf $output ${dest}

# Make dated folder inside destination
f_date=$( TZ=":US/Eastern" date +%Y-%m-%d )

OUT="${dest}/${f_date}"

echo "INFO: Writing files into archive folder ${OUT}"
cp -rf ${output}/* ${OUT}
echo "Done"

