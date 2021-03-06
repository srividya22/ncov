from os import environ
from socket import getfqdn
from getpass import getuser

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

configfile: "config/config_local.yaml"

# simple rule to call snakemake for outsider users
rule all:
    input:
        final_sequences = "final/sequences.fasta",
        final_metadata = "final/metadata.tsv"


### For now just copy gisaid seq file; later update with scraping 
rule download_gisaid:
    message: "Downloading metadata and fasta files from GISAID"
    output:
        gi_sequences = "final/gisaid_raw.fasta",
        gi_metadata = "final/gisaid_raw_meta.tsv"
    shell:
        """
        ## TODO : Add download script from GISAID ; Temporary only copying
        cp "final/gisaid_cov2020_sequences.fasta" {output.gi_sequences}
        cp "final/gisaid_metadata_raw.tsv" {output.gi_metadata}
        #curl -o {output.gi_metadata:q} https://raw.githubusercontent.com/nextstrain/ncov/master/data/metadata.tsv
        
        """
        
rule normalize_gisaid:
    message:
        """
        Normalize GISAID sequences and metadata
        """
    input:
        n_sequences = rules.download_gisaid.output.gi_sequences,
        n_metadata = rules.download_gisaid.output.gi_metadata
    output:
        on_sequences = "final/gisaid_sequences.fasta",
        on_metadata = "final/gisaid_metadata.tsv"
    shell:
        """
        echo "Running Cmd : scripts/normalize_and_filter_gisaid_fasta.sh {input.n_sequences} {input.n_metadata} {output.on_sequences} {output.on_metadata}"
        scripts/normalize_and_filter_gisaid_fasta.sh {input.n_sequences} {input.n_metadata} {output.on_sequences} {output.on_metadata}
        """

rule normalize_local:
    message:
        """
        Parse local metadata to nextstrain format and Prepare JHU sequences
        """
    input:
        l_sequences = config["local_sequences"],
        l_metadata = config["local_gisaid_metadata"]
    params:
        min_length = 29902
    output:
        lf_sequences = "final/JHU_sequences.fasta",
        lf_tmp_meta = "final/JHU_nextstrain.tsv",
        lf_metadata = "final/JHU_metadata.tsv"
    shell:
        """
        echo "Running Cmd : scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -o {output.lf_tmp_meta}"
        scripts/parse_meta_from_GISAID_to_nextstrain.py -in {input.l_metadata} -g {input.l_sequences} -o {output.lf_tmp_meta}
        echo "Running Cmd : scripts/normalize_and_filter_gisaid_fasta.sh {input.l_sequences} {output.lf_tmp_meta} {output.lf_sequences} {output.lf_metadata} {params.min_length}"
        scripts/normalize_and_filter_gisaid_fasta.sh {input.l_sequences} {output.lf_tmp_meta} {output.lf_sequences} {output.lf_metadata}
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
        final_sequences = "final/sequences.fasta",
        final_metadata = "final/metadata.tsv"
    shell:
        """
        echo "Running Cmd : scripts/combine_seq.sh {input.gi_seq} {input.gi_metadata} {input.local_seq} {input.local_metadata} {output.final_sequences} {output.final_metadata}
"
        scripts/combine_seq.sh {input.gi_seq} {input.gi_metadata} {input.local_seq} {input.local_metadata} {output.final_sequences} {output.final_metadata}
        """
