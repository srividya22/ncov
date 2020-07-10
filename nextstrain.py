#!/usr/bin/env python

"""
Script to configure nextstrain

"""

#Configuration steps
run_env="alpha / beta / dev"
input_fasta=""
input_meta=""
all_input_fasta = ""
all_input_meta = ""
outdir=""
deploy_dir=""
main_config=""
auspice_config=""
regions = ""
add_metadata=""

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
import yaml

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--submission_manifest",dest="sub_manifest",type=str, required=True, help="Path to the tsv file with run_id,sample_name,sample_header,rename_header,coveragex")
    parser.add_argument("--seq-path",dest="SEQ_PATH",type=str, required=True, help="Path to the sequencing runs folder")
    parser.add_argument("--fasta-path",dest="FASTA_PATH",type=str, required=True, help="Path to the draft consensus fasta folder")
    parser.add_argument("--submitted-fasta",dest="SUBMITTED_SEQ",type=str, required=True, help="Path to a single fasta of already submitted fasta files")
    parser.add_argument("--master-meta",dest="MASTER_META",type=str, required=True, help="File containing the master metadata for GISAID submission")
    parser.add_argument("-out",dest="OUTPUT_DIR",type=str, required=True, help="Path to output directory for alpha nextstrain")
    
    args = parser.parse_args()
    with open("my_file.yaml") as f:
         list_doc = yaml.load(f)
 
    for sense in list_doc:
        if sense["name"] == "sense2":
            sense["value"] = 1234

    with open("my_file.yaml", "w") as f:
         yaml.dump(list_doc, f)
