#!/usr/bin/env bash

output=$1
dest=$2


# Make dated folder inside destination
f_date=$( TZ=":US/Eastern" date +%Y-%m-%d )

OUT="${dest}/${f_date}"

mkdir -p ${OUT}

echo "INFO: Writing files into archive folder ${OUT}
mv $output ${OUT}

echo "Done"

