#!/usr/bin/env python
"""
Script to prepare nextstrain alpha files
"""
###TODO
# Assumptions : Second column sample_name is unique by itself and only the latest run for the sample is retained in submission manifest
#1) Get files in the submission manifest from seq_path 
#2) Exclude submitted files as GISAID will already have it
#3) Add submitted files to local build only
#4) If metadata exists use metadata
#5) Else create fake metadata for the non existing samples


import pandas as pd
import numpy as np 
import json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize
import glob ,re , tempfile ,shutil
import random , time
from datetime import datetime , date , timedelta
from pyfaidx import Fasta
from datetime import date
import warnings
from pytz import timezone
from Bio import SeqIO

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

class PrepareAlphaNextstrain:
    """ 
    Prepare data for Nextstrain
    """
    
    def __init__(self,sub_manifest,SEQ_PATH,FASTA_PATH,MASTER_META,SUBMITTED_SEQ,OUTPUT_DIR):
        """
        Class to prepare alpha data consensus data for nextstrain
        
        :param sub_manifest: Submission Manifest file for a sequencing run
        :param SEQ_PATH: path to sequencing runs
        :param FASTA_PATH: path to fasta file within sequencing runs
        :param MASTER_META: path to master metadata excel sheet
        :param SUBMITTED_SEQ: path to fasta file of already submitted sequences
        :param OUTPUT_DIR: path to alpha input files 
        """
        self.id_mapping_file = sub_manifest
        self.SEQ_PATH = SEQ_PATH
        self.FASTA_PATH = FASTA_PATH
        self.MASTER_META = MASTER_META
        self.SUBMITTED_SEQ = SUBMITTED_SEQ
        self.OUTPUT_DIR = OUTPUT_DIR
        self.CONSENSUS_FASTA = "*.nanopolish.consensus.fasta"
        self.DATE_FMT ="%Y%m%d"
        self.SUB_FMT ="%Y-%m-%d"

        self.REF_LENGTH = 29903
        self.count_Ns = 5000
        
        self.CUR_DATE = datetime.now(timezone('US/Eastern'))
        self.SUB_DATE = self.CUR_DATE.strftime(self.SUB_FMT)
        self.CUR_DATE = self.CUR_DATE.strftime(self.DATE_FMT)

        self.ALPHA_PREFIX= os.path.join(self.OUTPUT_DIR, self.CUR_DATE)
        self.STATE_PREFIX = "MD"
        self.header_pat = re.compile('([A-Z]+-[0-9]*_[A-Z]+[0-9]+)')
       
        ### Setting default metadata for all samples
        self.loc_default = "North America / USA / Maryland"
        self.host = "Human"
        self.specimen_source = "Nasopharyngeal swab"
        self.seq_tech = "Nanopore MinION"
        self.assembly_method = "ARTIC 1.0.0"
        self.dep = "Johns Hopkins Hospital Department of Pathology"
        self.address = "600 N. Wolfe St; Baltimore, MD 21287"
        self.authors = "Peter M. Thielen, Thomas Mehoke, Shirlee Wohl, Srividya Ramakrishnan, Melanie Kirsche, Amanda Ernlund, Oluwaseun Falade-Nwulia, Timothy Gilpatrick, Paul Morris, Norah Sadowski, Nídiá Trovao, Victoria Gniazdowski, Michael Schatz, Stuart C. Ray, Winston Timp, Heba Mostafa"
        self.submitter = "JHU"
        self.next_fields = [ "Virus name",
                             "Type",
                             "Passage details/history",
                             "Collection date",
                             "Location",
                             "Additional location information",
                             "Host",
                             "Additional host information",
                             "Gender",
                             "Patient age",
                             "Patient status",
                             "Specimen source",
                             "Outbreak",
                             "Last vaccinated",
                             "Treatment",
                             "Sequencing technology",
                             "Assembly method",
                             "Coverage",
                             "Originating lab",
                             "Address",
                             "Sample ID given by the sample provider",
                             "Submitting lab",
                             "Address.1",
                             "Sample ID given by the submitting laboratory",
                             "Authors",
                             "Submitter",
                             "Submission date",
                             "Run date",
                             "Run folder" ]
        
    def log(self,message):
        """ Log messages to standard output. """
        print(time.ctime() + ' --- ' + message, flush=True)
   
    def get_fasta_lengths(self,fasta):
        fa_lens = []
        N_counts = []
        for seq_record in SeqIO.parse(fasta, "fasta"):
               A_count = seq_record.seq.count('A')
               C_count = seq_record.seq.count('C')
               G_count = seq_record.seq.count('G')
               T_count = seq_record.seq.count('T')
               N_count = seq_record.seq.count('N')
               fa_lens.append(len(seq_record))
               N_counts.append(N_count)
        return fa_lens, N_counts

    def getKeysByValues(self,dictOfElements, listOfValues):
        listOfKeys = list()
        listOfItems = dictOfElements.items()
        for item  in listOfItems:
            if item[1] in listOfValues:
                listOfKeys.append((item[1], item[0]))
        return  listOfKeys
    
    def get_file_list(self,file_list,file_pattern):
        """
        get file list of consensus fasta from SEQ_PATH and FASTA_PATH
    
        """
        c_list = []
        for i in file_list:
            cfile = glob.glob(i + self.CONSENSUS_FASTA)[0]
            c_list.append(cfile)
        return c_list
        
    def get_fasta_header(self,seq):
        """
         get header from fasta file
        """
        with open(seq,'r') as fd:
           header = []
           for line in fd:
               if line[0]=='>':
                   header.append(line.strip().split(">")[1].split(" ")[0])
        return header
    
    def map_submission_names(self):
        """ 
        Extract submission names from submission_manifest.txt
        """
        sub_file = self.id_mapping_file
        self.sname_dict= {}
        self.sheader_dict = {}
        self.only_sname_dict = {}
        self.only_sheader_dict = {}
        with open(sub_file,'r') as fs:
             for line in fs: 
                 run_id, sample_name, sample_header , renamed_header, coveragex = line.strip().split("\t")  
                 s_key = run_id + ":" + sample_name
                 h_key = run_id + ":" + sample_header
                 sn_key = sample_name
                 sh_key = sample_header
                 if  s_key not in self.sname_dict:
                     self.sname_dict[s_key] = renamed_header
                 if  h_key not in self.sheader_dict:
                     self.sheader_dict[h_key] = renamed_header
                 if  sn_key not in self.only_sname_dict:
                     self.only_sname_dict[sn_key] = renamed_header
                 if  sh_key not in self.only_sheader_dict:
                     self.only_sheader_dict[sh_key] = renamed_header
        return ( self.sname_dict , self.sheader_dict, self.only_sname_dict, self.only_sheader_dict )
    
    def concat_fasta_files(self,fa_list, out):
        """
        Function to concat fasta files from a list
        """
        filt_fnames = [] 
        with open(out, 'w') as outfile:
             for fname in fa_list:
                 ffile = glob.glob(fname)[0]
                 g_lengths, n_counts = self.get_fasta_lengths(ffile)
                 if g_lengths[0] == 29903 and n_counts[0] < 5000:
                    filt_fnames.append(ffile)
                    with open(ffile) as infile:
                         for line in infile:
                             outfile.write(line)
                 else:
                    self.log("INFO : Skipping {0} from alpha Genome Length :  {1} , N_counts : {2}".format(fname,g_lengths[0],n_counts[0]))
        return filt_fnames
 
    def filter_metadata(self,sample_list,out,fill_meta=None):
        """
        Function to filter metadata from a list
        """
        meta_tsv = pd.read_table(self.MASTER_META)
        existing_meta_list = meta_tsv["Virus name"].tolist()
        self.meta_columns = meta_tsv.columns.values.tolist()
        self.no_meta_ids = [ i for i in sample_list if i not in existing_meta_list ]
        meta_tsv = meta_tsv[ meta_tsv["Virus name"].isin(sample_list) ] 
        # Change all unknown collection date to run_date
        for index, row in meta_tsv.iterrows():
            if meta_tsv.loc[index,'Collection date'] in ["", "?", "??","unknown", None] :
               #(meta_tsv.loc[index,'Collection date'])
               meta_tsv.loc[index,'Collection date'] = meta_tsv.loc[index,'Run date']
            if meta_tsv.loc[index,'Submission date'] in ["", "?", "??","unknown", None] :
               meta_tsv.loc[index,'Submission date'] = self.SUB_DATE
            #if meta_tsv.loc[index,'Location'].astype(str).str.contains("??|unknown") :
            #   meta_tsv.loc[index,'Location'].replace("??",)
        meta_tsv[['Collection date', 'Submission date']] = meta_tsv[['Collection date','Submission date']].fillna(value=self.SUB_DATE)


        # Change all unknown State to MT (  unknown )
        meta_tsv.to_csv(out, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
        return out
    
    def rename_submission_fasta(self,alpha_fa_file,rename_dict):
        """
        Get the GISAID submission names for all the files selected for submission
        
        """
        fa_file = Fasta(alpha_fa_file, key_function=lambda key: rename_dict[key])
        renamed_fa_file = alpha_fa_file + ".renamed"
        with open(renamed_fa_file, 'w') as renamed:
             for entry in fa_file:
                 print(">" + entry.name, file=renamed)
                 for line in entry:
                      print(line, file=renamed)
        shutil.copy(renamed_fa_file, alpha_fa_file)
        if os.path.isfile(renamed_fa_file):
            os.remove(renamed_fa_file)
        fa_file = Fasta(alpha_fa_file)
        return alpha_fa_file
    
    def prepare_fasta(self):
        """
        Function to the final Alpha nextstrain fasta files
        - Get sname_dict  , sheader_dict mappings from submission manifest
        - Get existing_seq , headers of the submitted files
        - If key in already submitted , remove the file from alpha_file_list; concat files
        - append submitted files to a new file all_jhu_sequences.fasta
        """
        self.sname_dict , self.sheader_dict , self.only_sname_dict , self.only_header_dict =  self.map_submission_names()
        alpha_tuple_list = [ (i.split(":")[0] , i.split(":")[1]) for i in self.sname_dict.keys() ] 
        self.alpha_file_list = [x[1] for x in alpha_tuple_list]

        self.existing_seq = self.get_fasta_header(self.SUBMITTED_SEQ)
        self.existing_snames = [ name.split(":")[1] for name, rname in self.sname_dict.items() if rname in self.existing_seq ]
        self.alpha_file_list = [e for i, e in enumerate(self.alpha_file_list) if e not in self.existing_snames]
        alpha_file_path_list = [ self.SEQ_PATH + "/" + r + "/" + self.FASTA_PATH + "/" + s for (r,s) in alpha_tuple_list if s in self.alpha_file_list ] 
        alpha_file_path_list = self.get_file_list(alpha_file_path_list, self.CONSENSUS_FASTA)
        

        self.outseq = self.ALPHA_PREFIX + "-jhu_sequences.fasta"
        self.all_outseq = self.ALPHA_PREFIX + "-all_local_jhu_sequences.fasta"
        filt_fnames = self.concat_fasta_files(alpha_file_path_list,self.outseq)
           
        ori_headers = self.get_fasta_header(self.outseq)
        #new_headers = [ i.rsplit("/",3)[1].split(".nanopolish")[0] for i in ori_headers ] 
        new_headers = {}
        for header in ori_headers:
            m = re.search(self.header_pat,header)
            new_head = m.group(0)
            if header not in new_headers:
                new_headers[header] = new_head
        
        #new_headers = { i : re.search(self.header_pat,i).group(0) for i in ori_headers } 
        #new_headers = [ i.split("4-draft-consensus/")[1] for i in new_headers ] 
        self.rename_submission_fasta(self.outseq,new_headers)                       
        # Rename header if not in self.only_header_dict
        ### This part of the code is to deal with names of format MDHP-18 to be renamed to MDHP-00018
        #for k ,i in new_headers.items():
        #    if i not in self.only_header_dict:
        #        parts = i.split("-")
        #        if len(parts[1].split("_")[0]) < 5:
        #            new_value = "-".join([parts[0],parts[1].split("_")[0].zfill(5)+"_"+parts[1].split("_")[1]])
        #            new_headers[k] = new_value
                    
        rename_dict =  { key : self.only_header_dict[key] if key in self.only_header_dict else key for key in new_headers.values() }
        self.rename_submission_fasta(self.outseq,rename_dict)                       
        self.all_jhu_list = [self.SUBMITTED_SEQ ,self.outseq ]
        self.concat_fasta_files(self.all_jhu_list,self.all_outseq)
        return ( self.outseq  , self.all_outseq )
        # return final fasta
        
    def prepare_metadata(self):
        """
        Function to get metadata for the fasta files
        - get all submitted seq metadata
        - Get strain ids , cut by sep "/" index 2, cut "-" [1] to get HP name
        - Append STATE prefix to HP name
        - If present. Filter rows based on ids 
        - If not present, create fake metadata based on run_date, default location (Maryland); unknown for everything else
        """
        sample_names = self.get_fasta_header(self.outseq)
        all_sample_names = self.get_fasta_header(self.all_outseq)
        self.outmeta = self.ALPHA_PREFIX + "-jhu_metadata.tsv"
        self.all_outmeta = self.ALPHA_PREFIX + "-all_local_jhu_metadata.tsv"
        self.filter_metadata(sample_names,self.outmeta)   
        self.filter_metadata(all_sample_names,self.all_outmeta)   
        if len(self.no_meta_ids) > 0:
            self.log("INFO : No metadata for ids : {0}".format(",".join(self.no_meta_ids)))
            self.log("INFO : Adding metadata for the ids")                    
            no_meta_mapped_list = [ (k, i.split(":")[1] , i.split(":")[0]) for k,i in self.getKeysByValues(self.sname_dict,self.no_meta_ids) ]
            missing_df_lists = [] 
            missing_df = pd.DataFrame(columns = self.meta_columns)
            self.add_meta_cols = [ i for i in missing_df.columns.tolist() if i not in self.next_fields ]
            #print("\n".join(self.add_meta_cols))
            add_rows = []
            for i,j,k in no_meta_mapped_list:
                
                run_date = datetime.strptime(k.split("_")[0],self.DATE_FMT)
                run_date = run_date.strftime(self.SUB_FMT)
                row = [ i , 
                       "betacoronavirus", 
                       "Original", 
                       run_date , 
                       self.loc_default , "unknown" , 
                       self.host,
                       "unknown", "unknown" , "unknown" , "unknown" , 
                       self.specimen_source ,
                       "unknown", "unknown" , "unknown" , 
                       self.seq_tech , self.assembly_method ,
                       "unknown" , self.dep  , self.address , j , self.dep  , self.address, j,
                       self.authors , self.submitter , self.SUB_DATE, run_date , k ]
                row = row + ["unknown" for i in self.add_meta_cols]
                add_rows.append(row)
                #missing_df.loc[len(missing_df)] = row
                #missing_df = missing_df.append(row)
                #missing_df = pd.concat([missing_df, pd.DataFrame(row)], axis=0)
                #missing_df.loc[ missing_df.index.max() +  1 ] =  row
                #missing_df[self.meta_columns] = row

                #missing_df = missing_df.append(pd.Series(row),index = self.meta_columns )
            missing_df = pd.DataFrame(add_rows, columns=self.meta_columns )    
            missing_df.to_csv(self.outmeta, mode='a', sep="\t",header=False)
            missing_df.to_csv(self.all_outmeta, mode='a', sep="\t",header=False)
          
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse new sequencing runs from alpha nextstrain")
    parser.add_argument("--submission_manifest",dest="sub_manifest",type=str, required=True, help="Path to the tsv file with run_id,sample_name,sample_header,rename_header,coveragex")
    parser.add_argument("--seq-path",dest="SEQ_PATH",type=str, required=True, help="Path to the sequencing runs folder")
    parser.add_argument("--fasta-path",dest="FASTA_PATH",type=str, required=True, help="Path to the draft consensus fasta folder")
    parser.add_argument("--submitted-fasta",dest="SUBMITTED_SEQ",type=str, required=True, help="Path to a single fasta of already submitted fasta files")
    parser.add_argument("--master-meta",dest="MASTER_META",type=str, required=True, help="File containing the master metadata for GISAID submission")
    parser.add_argument("-out",dest="OUTPUT_DIR",type=str, required=True, help="Path to output directory for alpha nextstrain")
    
    args = parser.parse_args()
    prepAlpha = PrepareAlphaNextstrain(args.sub_manifest,args.SEQ_PATH,args.FASTA_PATH,args.MASTER_META,args.SUBMITTED_SEQ,args.OUTPUT_DIR)
    try:
        prepAlpha.prepare_fasta()
        prepAlpha.prepare_metadata()
    except:
        raise
