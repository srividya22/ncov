#!/usr/bin/env python
"""
Parse GISAID metadata to add additional columns as required

"""

import pandas as pd
import numpy as np 
import json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse GISAID metadata files to add additional columns to the metadata")
    parser.add_argument("-in",dest="xfile",type=str, required=True, help="Path to the GISAID tsv")
    parser.add_argument("-l",dest="meta_list",type=str, required=True, help="Path to the GISAID tsv")
    parser.add_argument("-o",dest="out", type=str, default="gisaid_metadata.tsv", help="Output file path")

    args = parser.parse_args()
    gisaid_file = args.xfile
    meta_list = args.meta_list
    out_file = args.out
    
    data_tsv = pd.read_table(gisaid_file,header=0,na_filter= False)
    with open(args.meta_list) as file:
         new_cols = [line.strip() for line in file]
    for col in new_cols:
        data_tsv[col] = "unknown"
    #data_tsv['run_date'] = 'unknown'
    #data_tsv['run_id'] = 'unknown'
    print("INFO : Added additional metadata to existing metadata file")
    data_tsv.to_csv(out_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
