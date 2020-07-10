#!/usr/bin/env python
"""
Parse GISAID metadata to add additional columns as required

"""

import pandas as pd
import numpy as np 
import json , os,re
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse GISAID metadata files to add additional columns to the metadata")
    parser.add_argument("-in",dest="xfile",type=str, required=True, help="Path to the GISAID tsv")
    parser.add_argument("-l",dest="meta_list",type=str, required=True, help="Path to the GISAID tsv")
    parser.add_argument("-m",dest="all_jhu_meta",type=str, required=True, help="Path to all jhu tsv")
    parser.add_argument("-o",dest="out", type=str, default="gisaid_metadata.tsv", help="Output file path")

    args = parser.parse_args()
    gisaid_file = args.xfile
    meta_list = args.meta_list
    all_jhu_meta = args.all_jhu_meta
    out_file = args.out
    
    data_tsv = pd.read_table(gisaid_file,header=0,index_col=0,na_filter= False)
    all_jhu_tsv = pd.read_table(all_jhu_meta,header=0,index_col=0,na_filter= False)
    print(data_tsv.columns.tolist())
    print(all_jhu_tsv.columns.tolist())

    with open(args.meta_list) as file:
         new_cols = [line.strip() for line in file]
    #for col in new_cols:
    #    data_tsv[col] = "unknown"
    #data_tsv['run_date'] = 'unknown'
    #data_tsv['run_id'] = 'unknown'
    print(data_tsv.shape)
    print(data_tsv.index)
    print(all_jhu_tsv.shape)
    print(all_jhu_tsv.index)
    data_tsv = data_tsv.join(all_jhu_tsv, how='left', lsuffix='_left', rsuffix='_right')
    data_tsv = data_tsv.drop(data_tsv.filter(regex='_right').columns, axis=1)
    data_tsv = data_tsv.rename(columns=lambda x: re.sub('_left','',x))

    #data_tsv = pd.merge(data_tsv,all_jhu_tsv[['strain'] + new_cols],on='strain', how='outer').fillna("unknown",inplace=True)
    #data_tsv = pd.merge(data_tsv,all_jhu_tsv[[ 'strain' ] + new_cols], left_on = 'strain' ,right_on ='strain' , how='left').fillna("unknown",inplace=True)
    #data_tsv = pd.merge(data_tsv,all_jhu_tsv[['strain'] + new_cols],on='strain', how='left').fillna("unknown",inplace=True)
    print(data_tsv.shape)
    #data_tsv = data_tsv.merge(all_jhu_tsv[ ['strain'] + new_cols], on='strain',how="inner")
    data_tsv.fillna(value="unknown", inplace = True)
    #data_tsv = data_tsv.loc[~data_tsv['date'].isin(["unknown"])]

    data_tsv = data_tsv[data_tsv['date'] != 'unknown']
    print(data_tsv.head)

    print("INFO : Added additional metadata to existing metadata file")
    data_tsv.to_csv(out_file, sep='\t', encoding='utf-8',  index=True, line_terminator='\r\n')
