#!/usr/bin/env python
"""
Parse GISAID metadata to nextstrain meta tsv 

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
JHU_URL = 'https://www.jhu.edu/'

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
    parser.add_argument("-in",dest="xfile",type=str, required=True, help="Path to the GISAID tsv")
    parser.add_argument("-g",dest="fa_file",type=str, required=True, help="Path to the fasta file from GISAID")
    parser.add_argument("-o",dest="out", type=str, default="gisaid_metadata.tsv", help="Output file path")
    #parser.add_argument("-x",dest="exclude", type=str, default="exclude_samples.txt", help="Output exclude file path")

    args = parser.parse_args()
    gisaid_file = args.xfile
    fa_file = args.fa_file
    out_file = args.out
    #ex_file = args.exclude
    
    #data_tsv = pd.read_excel(gisaid_file, 'Acknowledgement Table', skiprows=list(range(2)))
    data_tsv = pd.read_table(gisaid_file,header=0,na_filter= False)
    
    # Drop empty rows in excel sheet
    data_tsv = data_tsv.dropna(axis=0, how='any', thresh=None, subset=None, inplace=False)
    # Strip whitespaces from columns
    data_tsv = data_tsv.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Drop rows with empty virus name 
    data_tsv['Virus name'] = data_tsv['Virus name'].str.replace(" ","")
    data_tsv['Virus name'] = data_tsv['Virus name'].str.replace("hCoV-19/","")
    data_tsv['Virus name'] = data_tsv['Virus name'].str.replace("hcoV-19/","")
    #print(data_tsv.Location.str.rsplit("/",n=3))
    loc_df = data_tsv.Location.str.split("/",n=3,expand=True)
    data_tsv['region'] = loc_df[0]
    data_tsv['country'] = loc_df[1]
    data_tsv['division'] = loc_df[2]
    if loc_df.shape[1] == 4 :
       data_tsv['location'] = loc_df[3]
    else:
       data_tsv['location'] = loc_df[2] # if location is unknown set division to location
    
    #Uncomment the following for run specific nextstrain
    data_tsv['Run date'] = pd.to_datetime(data_tsv['Run date'], format='%Y-%m-%d', errors='coerce')
    most_recent_date = data_tsv['Run date'].max()
    
    #data_tsv[['region','country','division','location']] = region,country,division,location
    #data_tsv[['region','country','division','location']] = data_tsv.Location.apply(lambda x: pd.Series(str(x).split("/"))) 
    
    #data_tsv[['region','country','division']] = data_tsv.Location.str.split("/",n=2,expand=True)
    #data_tsv['location'] = None
    # Create new dataframe
    old_data_tsv = data_tsv
        #print(data_tsv.head(5))
    # Delete unneccessary columns after splitting
    #data_tsv.drop(['Location'], axis=1, inplace=True)
    DROP_COLS = [ 'Type' , 
                  'Passage details / history' ,
                  'Location' , 
                  'Additional location information' , 
                  'Additional host information', 
                  'Patient status' ,
                  'Specimen source' , 
                  'Outbreak detail' ,
                  'Last vaccinated' ,
                  'Treatment' , 
                  'Sequencing technology',
                  'Assembly method', 
                  'Coverage',
                  'Address',
                  'Address.1' ,
                  'Sample ID given by the Submitting lab' , 
                  'Submitter' , 
                  'Run date' , 
                  'Run folder']
    
    data_tsv.drop(DROP_COLS, axis = 1,inplace=True) 
        
    # Rename columns to nextstrain columns
    data_tsv=data_tsv.rename(columns = {'Virus name' : 'strain',
                                               #'Accession ID':'gisaid_epi_isl' ,
                                               'Originating lab' : 'originating_lab', 
                                               'Submitting lab' : 'submitting_lab', 
                                               'Authors' : 'authors' ,
                                               'Collection date' : 'date' ,
                                               'Patient age' : 'age' ,
                                               'Gender' : 'sex' , 
                                               'Submission date' : 'date_submitted' ,
                                               'Host' : 'host'
                                                })
    
    #data_tsv['gisaid_epi_isl'] = data_tsv['gisaid_epi_isl'].fillna(data_tsv['Sample ID given by the provider'])
    #data_tsv.gisaid_epi_isl.combine_first(data_tsv['Sample ID given by the provider'])
    #data_tsv['gisaid_epi_isl']=data_tsv['gisaid_epi_isl'].mask(pd.isnull, data_tsv['Sample ID given by the provider'])
    data_tsv['gisaid_epi_isl'] = data_tsv['Sample ID given by the provider']
    #print(data_tsv['Sample ID given by the provider'])
    #print(data_tsv['gisaid_epi_isl'].tolist())
    data_tsv.drop('Sample ID given by the provider',axis =1 , inplace = True)
    data_tsv['gisaid_epi_isl'] = ''
    data_tsv['virus'] = 'ncov'
    data_tsv['genbank_accession'] = "?"
    data_tsv['region_exposure'] = data_tsv['region']
    data_tsv['country_exposure'] = data_tsv['country']
    data_tsv['division_exposure'] = data_tsv['division']
    data_tsv['segment'] = 'genome'
    print("Total Number of Samples : " + str(data_tsv.shape[0]))
        
    data_tsv['date'] = pd.to_datetime(data_tsv['date'], format='%Y-%m-%d', errors='coerce')
    data_tsv = data_tsv[data_tsv['date'].between(EPI_START,CUR_DATE )]
    #print("Number of filtered Samples : " + str(data_tsv.shape[0]))
        
    #data_tsv = data_tsv.query("@EPI_START < date <= @CUR_DATE")
        
        # Not currently in gisaid
    #f_lens = get_fasta_lengths(fa_file)
    #data_tsv['length'] = [ f_lens[i] if i in f_lens else 0 for i in data_tsv['strain'] ] #ADD Lengths from fasta
    data_tsv['length'] = 29903
    #print("Number of filtered Samples : " + str(data_tsv.shape[0]))
    #data_tsv['host'] = 'Human'
        
    # Not currently in GISAID xls
    #data_tsv['age'] =  [ random.choice(AGE) for i in range(len(data_tsv['segment'])) ]
    #data_tsv['sex'] = [ random.choice(SEX) for i in range(len(data_tsv['segment'])) ]
                  
    # GISAID url
    data_tsv['url'] = JHU_URL
    data_tsv['title'] = 'unknown'
    #data_tsv['date_submitted'] = data_tsv['date']
    #data_tsv['date_submitted'] = [ i + timedelta(weeks=2) for i in data_tsv['date'] ]
    data_tsv.fillna("unknown", inplace=True)
    data_tsv = data_tsv.applymap(lambda x: x.strip() if isinstance(x, str) else x)
            
    #data_tsv.columns
        
    #data_tsv = data_tsv[[data_tsv.columns[i] for i in ORDER]]
    #print(len(data_tsv.columns.tolist()))
    #print(data_tsv.columns.tolist())
    #print(len(COL_ORDER))
    
    
    if (len(COL_ORDER) == len(data_tsv.columns.tolist())):
        data_tsv = data_tsv[COL_ORDER]
    else: 
        assert len(COL_ORDER) != len(data_tsv.columns.tolist()), 'Problem processing metadata fields from GISAID ; generated only :' + str(len(data_tsv.columns))
    
    print("Processed Samples : " +  str(data_tsv.shape[0]) + " ; Metadata fields : " + str(len(data_tsv.columns.tolist())) )    
    data_tsv.to_csv(out_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
