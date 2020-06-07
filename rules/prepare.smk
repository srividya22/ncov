from os import environ
from socket import getfqdn
from getpass import getuser

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

configfile: "config/config_beta.yaml"

# Number of nextstrain columns
COLS_NEXTSTRAIN=23

if "outdir" not in config:
    OUTDIR = os.getcwd()
else:
    OUTDIR = config["outdir"]

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# simple rule to call snakemake for outsider users
rule all_data:
    input:
        final_sequences = OUTDIR + "final/sequences.fasta",
        final_metadata = OUTDIR + "final/metadata.tsv" ,
        local_sequences = OUTDIR + "final/local_sequences.fasta",
        local_metadata = OUTDIR + "final/local_metadata.tsv"


### For now just copy gisaid seq file; later update with scraping 
rule download_gisaid:
    message: "Downloading metadata and fasta files from GISAID"
    input:
        gi_in_seq = OUTDIR + "final/gisaid_cov2020_sequences.fasta",
        gi_in_meta = OUTDIR + "final/gisaid_metadata_raw.tsv"
    output:
        gi_sequences = OUTDIR + "final/gisaid_raw.fasta",
        gi_metadata = OUTDIR + "final/gisaid_raw_meta.tsv"
    shell:
        """
        ## TODO : Add download script from GISAID ; Temporary only copying
        cp {input.gi_in_seq} {output.gi_sequences}
        cp {input.gi_in_meta} {output.gi_metadata}
        #curl -o {output.gi_metadata:q} https://raw.githubusercontent.com/nextstrain/ncov/master/data/metadata.tsv
        
        """

rule normalize_local:
    message:
        """
        Parse local metadata to nextstrain format and Prepare JHU sequences
        """
    input:
        l_sequences = config["local_sequences"],
        l_metadata = config["local_gisaid_metadata"]
        all_sequences = config["all_local_sequences"],
        all_metadata = config["all_local_gisaid_metadata"]
    params:
        min_length = 29902,
        add_metadata = config["add_metadata"],
        next_cols = 23
        #next_cols = COLS_NEXTSTRAIN
    output:
        lf_sequences = OUTDIR + "final/JHU_sequences.fasta",
        lf_tmp_meta = OUTDIR + "final/JHU_nextstrain.tsv",
        all_lf_tmp_meta = OUTDIR + "final/all_JHU_nextstrain.tsv",
        lf_metadata = OUTDIR + "final/JHU_metadata.tsv" ,
        lf_add_meta = OUTDIR + "final/add_metadata_list.txt"
    shell:
        """
        #echo "Running Cmd : scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -o {output.lf_tmp_meta}"
        echo "Add_metadata : {params.add_metadata}"
        if [ {params.add_metadata} = "yes" ]; then
           echo "Running Cmd : scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -m "yes" -o {output.lf_tmp_meta}"
           scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -m "yes" -o {output.lf_tmp_meta}
           scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.all_metadata} -g {input.all_sequences} -m "yes" -o {output.all_lf_tmp_meta}
           echo "Running Cmd: scripts/get_add_meta_columns.sh {output.lf_tmp_meta} {params.next_cols} {output.lf_add_meta}"
           scripts/get_add_meta_columns.sh {output.lf_tmp_meta} {params.next_cols} {output.lf_add_meta}
        else
           echo "Running Cmd : scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -m "no" -o {output.lf_tmp_meta}"
           scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -m "no" -o {output.lf_tmp_meta}
           scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.all_metadata} -g {input.all_sequences} -m "no" -o {output.all_lf_tmp_meta}
        fi
        echo "Running Cmd : scripts/normalize_and_filter_gisaid_fasta.sh {input.l_sequences} {output.lf_tmp_meta} {output.lf_sequences} {output.lf_metadata} {params.min_length}"
        scripts/normalize_and_filter_gisaid_fasta.sh {input.l_sequences} {output.lf_tmp_meta} {output.lf_sequences} {output.lf_metadata}
        #if [ params.add_metadata == "yes" ]; then
        #   echo "Adding additional metadata run_data and run_id"
        #   scripts/add_additional_metadata.py -in {output.lf_metadata} -o {output.lf_metadata}
        #fi
        """
        
rule normalize_gisaid:
    message:
        """
        Normalize GISAID sequences and metadata
        """
    input:
        n_sequences = rules.download_gisaid.output.gi_sequences,
        n_metadata = rules.download_gisaid.output.gi_metadata,
        meta_list = rules.normalize_local.output.lf_add_meta
    params:
        add_metadata = config['add_metadata'] ## "yes" or "no"
    output:
        on_sequences = OUTDIR + "final/gisaid_sequences.fasta",
        on_metadata = OUTDIR + "final/gisaid_metadata.tsv"
    shell:
        """
        echo "Running Cmd : scripts/normalize_and_filter_gisaid_fasta.sh {input.n_sequences} {input.n_metadata} {output.on_sequences} {output.on_metadata}"
        scripts/normalize_and_filter_gisaid_fasta.sh {input.n_sequences} {input.n_metadata} {output.on_sequences} {output.on_metadata}
        if [ {params.add_metadata} = "yes" ]; then
           echo "Adding additional metadata run_data and run_id"
           scripts/add_additional_metadata.py -in {output.on_metadata} -l {input.meta_list} -o {output.on_metadata}
        fi
        """
        
rule merge_seq:
    message:
        """
        Merge JHU sequences and local sequences
        """
    input:
        gi_seq = rules.normalize_gisaid.output.on_sequences,
        gi_metadata = rules.normalize_gisaid.output.on_metadata,
        local_seq = rules.normalize_local.output.lf_sequences,
        local_metadata = rules.normalize_local.output.lf_metadata
    output:
        final_sequences = OUTDIR + "final/sequences.fasta",
        final_metadata = OUTDIR + "final/metadata.tsv"
    shell:
        """
        echo "Running Cmd : scripts/combine_seq.sh {input.gi_seq} {input.gi_metadata} {input.local_seq} {input.local_metadata} {output.final_sequences} {output.final_metadata}
"
        scripts/combine_seq.sh {input.gi_seq} {input.gi_metadata} {input.local_seq} {input.local_metadata} {output.final_sequences} {output.final_metadata}
        """

rule run_only_jhu_seq:
    message:
        """
        Merge JHU sequences and reference.fasta
        """
    input:
        ref_seq = rules.normalize_gisaid.output.on_sequences,
        ref_metadata = rules.normalize_gisaid.output.on_metadata,
        jhu_seq = rules.normalize_local.output.all_sequences,
        jhu_metadata = rules.normalize_local.output.all_lf_tmp_metadata
        #jhu_seq = config['all_local_sequences'],
        #jhu_metadata = config['all_local_gisaid_metadata']
    params:
        include = config["files"]["include"]
    output:
        only_ref = OUTDIR + "final/only_ref_sequences.fasta" , 
        only_ref_meta = OUTDIR + "final/only_ref_metadata.tsv",
        local_sequences = OUTDIR + "final/local_sequences.fasta",
        local_metadata = OUTDIR + "final/local_metadata.tsv"
    shell:
        """
        ### Ref seq subsample metadata for Wuhan reference
        
        echo "Running Cmd : scripts/get_ref_sequences.sh {input.ref_seq} {input.ref_metadata} {params.include} {output.only_ref} {output.only_ref_meta}"
        scripts/get_ref_sequences.sh {input.ref_seq} {input.ref_metadata} {params.include} {output.only_ref} {output.only_ref_meta} 
        echo " Running Cmd  : scripts/combine_seq.sh {output.only_ref} {output.only_ref_meta} {input.jhu_seq} {input.jhu_metadata} {output.local_sequences} {output.local_metadata}"
        scripts/combine_seq.sh {output.only_ref} {output.only_ref_meta} {input.jhu_seq} {input.jhu_metadata} {output.local_sequences} {output.local_metadata}
        """
