

## Install Dependencies : 

1) Nextstrain Local installation : https://nextstrain.org/docs/getting-started/local-installation </br>
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml </br> 
conda env create -f nextstrain.yml </br>
conda activate nextstrain </br> 
npm install --global auspice </br>
nextstrain check-setup --set-default </br>
</br>

2) Install seqtk : https://github.com/lh3/seqtk </br>
conda install -c bioconda seqtk </br>

## Customizing Local Nextstrain build

git clone https://github.com/srividya22/ncov.git </br>
cd ncov </br>
chmod 755 setup_local.sh </br>
./setup_local.sh </br>

-Update local file paths in config/Snakefile_local.file </br>

## Running Local Nextstrain

## Download GISAID data from GISAID https://www.gisaid.org/
mkdir -p final </br>
cp gisaid_cov2020_sequences.fasta final/gisaid_cov2020_sequences.fasta </br>

## Augur build
snakemake -j 4 --forceall --snakefile Snakefile_prepare </br>
snakemake -j 4 --forceall --snakefile Snakefile_main </br>

## View nextstrain
nextstrain view ./auspice

## To update nextstrain/ncov latest code
git remote add upstream https://github.com/nextstrain/ncov.git</br>
git fetch upstream </br>
git pull upstream master </br>

## Todo List
- Automate GISAID data download in Snakemake_prepare
- Add Snakefile for deploy and archive 
- Additional validation on GISAID metadata

## Setup Remote 

## About  

This is a [Nextstrain](https://nextstrain.org) build for novel coronavirus, alternately known as hCoV-19 or SARS-CoV-2 that is responsible for the COVID-19 pandemic. This build is visible at [nextstrain.org/ncov](https://nextstrain.org/ncov).

### Data

The hCoV-19 / SARS-CoV-2 genomes were generously shared via GISAID. We gratefully acknowledge the Authors, Originating and Submitting laboratories of the genetic sequence and metadata made available through GISAID on which this research is based. For a full list of attributions please see the [metadata file](data/metadata.tsv).

### Bioinformatics notes

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of [8 &times; 10^-4 subs per site per year](http://virological.org/t/phylodynamic-analysis-176-genomes-6-mar-2020/356). There were SNPs present in the nCoV samples in the first and last few bases of the alignment that were masked as likely sequencing artifacts.

If you'd like to customize and run the analysis yourself, please see the [developer docs](./DEV_DOCS.md).

## How to run using your own data

Please see
- [How to format the metadata](./docs/metadata.md)
- [Running a SARS-CoV-2 analysis](./docs/running.md)


## Contributing

We welcome contributions from the community! Please note that we strictly adhere to the [Contributor Covenant Code of Conduct](https://github.com/nextstrain/.github/blob/master/CODE_OF_CONDUCT.md).

### Contributing to translations of our situation reports  
Please see the [translations repo](https://github.com/nextstrain/translations) to get started!

### Contributing to software or documentation   
Please see our [Contributor Guide](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md) to get started!

### Contributing data  
**Please note that we automatically pick up any hCoV-19 data that is submitted to GISAID.**  

If you're a lab and you'd like to get started sequencing, please see:  
* [Protocols from the ARTIC network](https://www.protocols.io/groups/artic/publications)  
* [Funding opportunities for sequencing efforts](https://twitter.com/firefoxx66/status/1242147905768751106)  
* Or, if these don't meet your needs, [get in touch](mailto:hello@nextstrain.org)  
