#!/usr/bin/env python
"""
Parse excel to tsv 

"""

## Updates:
# Currently only works with acknowlegement table and normalized fasta from GISAID.
# Acknowledgement table has missing columns age, sex, region_exposure, country_exposure , submission_data
# Currently imputes age to be a random value between 1 - 89
# Currently imputes submission date to be 2 weeks from sample collection data
# Region exposure and country exposure to region and country values  ( only very few samples have updated values )

# Validations performed
# Drops all samples without valid date
# Drops all samples which doesnt have a corresponding normalized fasta sequence in the fasta provided

import pandas as pd
import numpy as np 
import json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize
from Bio import SeqIO


#Order of metadata in nextstrain as of 04/07/2020
COL_ORDER =['strain','virus','gisaid_epi_isl','genbank_accession','date','region','country','division','location','region_exposure','country_exposure','division_exposure','segment','length','host','age','sex','originating_lab','submitting_lab','authors','url','title','date_submitted']

#df = df[[df.columns[i] for i in order]]
AGE_RANGE = range(1,89)
AGE = [ i for i in AGE_RANGE ]
AGE.append('unknown')

SEX = ['Male' , 'Female' , 'unknown']
GISAID_URL = 'https://www.gisaid.org/'
CUR_DATE = pd.datetime.today()
EPI_START =  datetime.strptime('2020-01-03', '%Y-%m-%d')

def get_fasta_lengths(fasta):
    fa_lens = {} 
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if seq_record.id not in fa_lens:
           fa_lens[seq_record.id] =  len(seq_record)
           print(seq_record.id+"\t"+str(len(seq_record))) 
        else:
           print("Duplicate record found for " + seq_record.id )
    return fa_lens

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse augur json files")
    parser.add_argument("-in",dest="xfile",type=str, required=True, help="Path to the Excel file from GISAID")
    parser.add_argument("-g",dest="fa_file",type=str, required=True, help="Path to the fasta file from GISAID")
    parser.add_argument("-o",dest="out", type=str, default="gisaid_metadata.tsv", help="Output file path")
    args = parser.parse_args()
    xls_file = args.xfile
    fa_file = args.fa_file
    out_file = args.out
    data_xlsx = pd.read_excel(xls_file, 'Acknowledgement Table', skiprows=list(range(2)))
    # Drop empty rows in excel sheet
    data_xlsx = data_xlsx.dropna(axis=0, how='any', thresh=None, subset=None, inplace=False)
    
    # Strip whitespaces from columns
    data_xlsx = data_xlsx.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Drop rows with empty virus name 
    data_xlsx['Virus name'] = data_xlsx['Virus name'].str.replace(" ","")
    data_xlsx['Virus name'] = data_xlsx['Virus name'].str.replace("hCoV-19/","")
    
    data_xlsx[['region','country','division','location']] = data_xlsx.Location.str.split("/",n=3,expand=True)
    
    # Create new dataframe
    old_data_xlsx = data_xlsx
    #print(data_xlsx.head(5))
    # Delete unneccessary columns after splitting
    data_xlsx.drop(['Location'], axis=1, inplace=True)
    
    # Rename columns to nextstrain columns
    data_xlsx=data_xlsx.rename(columns = {'Virus name' : 'strain',
                                           'Accession ID':'gisaid_epi_isl' ,
                                           'Originating lab' : 'originating_lab', 
                                           'Submitting lab' : 'submitting_lab', 
                                           'Authors' : 'authors' ,
                                           'Collection date' : 'date' })
    
    data_xlsx['virus'] = 'ncov'
    data_xlsx['genbank_accession'] = "?"
    data_xlsx['region_exposure'] = data_xlsx['region']
    data_xlsx['country_exposure'] = data_xlsx['country']
    data_xlsx['division_exposure'] = data_xlsx['division']
    data_xlsx['segment'] = 'genome'
    print("Total Number of Samples : " + str(data_xlsx.shape[0]))
    
    data_xlsx['date'] = pd.to_datetime(data_xlsx['date'], format='%Y-%m-%d', errors='coerce')
    data_xlsx = data_xlsx[data_xlsx['date'].between(EPI_START,CUR_DATE )]
    #print("Number of filtered Samples : " + str(data_xlsx.shape[0]))
    
    #data_xlsx = data_xlsx.query("@EPI_START < date <= @CUR_DATE")
    
    # Not currently in gisaid
    f_lens = get_fasta_lengths(fa_file)
    data_xlsx['length'] = [ f_lens[i] if i in f_lens else 0 for i in data_xlsx['strain'] ] #ADD Lengths from fasta
    data_xlsx = data_xlsx[data_xlsx.length != 0]
    print("Number of filtered Samples : " + str(data_xlsx.shape[0]))
    data_xlsx['host'] = 'Human'
    
    # Not currently in GISAID xls
    data_xlsx['age'] =  [ random.choice(AGE) for i in range(len(data_xlsx['segment'])) ]
    data_xlsx['sex'] = [ random.choice(SEX) for i in range(len(data_xlsx['segment'])) ]
              
    # GISAID url
    data_xlsx['url'] = GISAID_URL
    data_xlsx['title'] = 'unknown'
    #data_xlsx['date_submitted'] = data_xlsx['date']
    data_xlsx['date_submitted'] = [ i + timedelta(weeks=2) for i in data_xlsx['date'] ]
    data_xlsx.fillna("unknown", inplace=True)
    data_xlsx = data_xlsx.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    
    
    #data_xlsx.columns
    
    #data_xlsx = data_xlsx[[data_xlsx.columns[i] for i in ORDER]]
    
    if (len(COL_ORDER) == len(data_xlsx.columns)):
        data_xlsx = data_xlsx[COL_ORDER]
    
    
    data_xlsx.to_csv(out_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
