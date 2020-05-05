#!/usr/bin/env python
"""

Filter vcf from snp-sites and report

"""

import pandas as pd
import numpy as np 
import json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize
import pysam 
from Bio import SeqIO
import logging

def get_ref_allele(ref_file,snps):
    ref_alleles = [] 
    for seq_record in SeqIO.parse(ref_file, "fasta"):
       for x in snps:
              # get the base at position x
              ref_alleles.append(seq_record.seq[x-1])
    return ref_alleles


def get_cov(total,count):
    if count > 0 : 
       return round((float(count) / float(total)) * 100  , 2)
    else:
       return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse augur json files")
    parser.add_argument("-in",dest="vcf",type=str, required=True, help="Path to the vcf from snp-sites")
    parser.add_argument("-r",dest="ref",type=str, required=True, help="Path to the reference sequence")
    parser.add_argument("-l",dest="exclude_samples",type=str, required=True, help="Path to list of samples to exclude")
    parser.add_argument("-m",dest="clades",type=str,required=True, help="Path to tab limited file of major clade specific SNPs")
    parser.add_argument("-c", metavar="15",dest="cov", type=float, default=15, help='Coverage threshold to determine SNP importance')
    parser.add_argument("-o",dest="out", type=str, default="output filtered vcf", help="Output vcf file")
    args = parser.parse_args()
    vcf_file = args.vcf
    exclude_samples = args.exclude_samples
    out_file = args.out
    ref = args.ref
    clades = args.clades
    cov = float(args.cov)
    
    # Count coverage based on the bases
    nt = ['A' ,'T' ,'C' ,'G']
    #CLADE_SPECIFIC_SITES = [ 514, 18060, 29095] 

    with open(vcf_file) as f1:
        with open(out_file, 'w') as out:
           head = [next(f1) for x in range(3)]
           out.write("".join(head))
    f1.close()
    out.close()

    # Read clades and major defining SNPS 
    clades = pd.read_table(clades)

    clade_nuc = clades.loc[clades['gene'] == 'nuc']['site'].to_list()
    print("Total Number of SNPs in major clades :" + str(len(clade_nuc)))

    snps_vcf = pd.read_table(vcf_file, skiprows=3)
    exclude_list = pd.read_table(exclude_samples,header=None)
    print("Total SNPs :" +str(snps_vcf.shape[0]) + " ; Samples : " + "\t" + str(snps_vcf.shape[1]-9) )

    ex_list =  [ x for x in exclude_list[0] if x in snps_vcf.columns ] 

    # get SNP postions
    snps = snps_vcf.POS.to_list()
    msa_alleles = snps_vcf.REF.to_list()
    msa_alt_alleles = snps_vcf.ALT.to_list()
    ref_alleles = get_ref_allele(ref,snps)
    ref_dict = dict(zip(snps ,ref_alleles))
    msa_dict = dict(zip(snps ,msa_alleles))
    msa_alt_dict = dict(zip(snps ,msa_alt_alleles))
    if  set(ref_dict.keys()) == set(msa_dict.keys()):
        swap_ref = {} 
        row_ids = []
        for key in set(msa_dict.keys()):
            if ref_dict[key] != msa_dict[key] :
                
               print ("INFO: REF ALLELE : " +  str(ref_dict[key]) + str(key) + " does not match MAJOR ALLELE : " + str(msa_dict[key] + str(key)))
               ref_allele_pos = msa_alt_dict[key].split(",").index(ref_dict[key]) + 1
               print("INFO : REF ALLELE : " + str(key) + ":" + str(ref_allele_pos) + " in " + msa_alt_dict[key] ) 
               swap_ref[key] = ref_allele_pos
               # Replace major allele to ALT allele
               old = msa_alt_dict[key]
               msa_alt_dict[key] = msa_alt_dict[key].replace(ref_dict[key],msa_dict[key])
               print( old + " to " + msa_alt_dict[key])
  
    #snps_vcf.drop(exclude_list[0],axis=1, inplace=True)
    snps_vcf.drop(ex_list,axis=1, inplace=True)
    snps_vcf['REF'] = ref_alleles
    snps_vcf['ALT'] = [ msa_alt_dict[i] for i in sorted(msa_alt_dict.keys())]

    # Filter SNP positions with just IUPAC codes

    snps_vcf = snps_vcf[snps_vcf.ALT.str.contains('A|C|G|T|a|c|g|t') == True]
    num_samples = (snps_vcf.shape[1]-9)

    counts= {}
    aln_counts = []
    m_allele = []
    alt_list = []
    alele_freq = []
    drop_snps = []
    atcg_cov_thr = [] 
    conf_flag = []

    for i in range(snps_vcf.shape[0]):
       snp_pos = snps_vcf['POS'].values[i]
       snp_ref = snps_vcf['REF'].values[i]
       snp_alt = snps_vcf['ALT'].values[i]
       snp_alt_alleles = snp_alt.split(',')
       counts = snps_vcf.iloc[i,9:snps_vcf.shape[1]].value_counts().to_dict()
       sorted_counts={k: counts[k] for k in sorted(counts)}
       
       # initiate conf score and flag to false
       c_flag = 'NO'
    
       # Check if the ref ALLELE is different switch counts
       if snp_pos in swap_ref:
          counts[swap_ref[snp_pos]] = counts[0]
        
       # Check length of alt allele counts is greater than 1
       if len(sorted_counts) >  1:
          maf_index = np.argmax([counts[i] for i in sorted_counts if i != 0])
          old_counts = ",".join([ str(counts[i]) for i in sorted_counts ])   
          counts = ",".join([ str(counts[i]) for i in sorted_counts if i != 0])
            
          # Update alt column if corresponding count is completely dropped to 0 after removing certain columns 
          if  len(snp_alt_alleles) != len(counts.split(',')):
               tmp_snp_alt = snp_alt
               snp_alt =  ",".join([ snp_alt_alleles[int(i-1)] for i in sorted_counts if i != 0 ])
               print("Updated Alleles based on counts  SNP : " +  str(snp_pos) +  ": " + snp_alt + " ; " + tmp_snp_alt )
 
               # Keep only alt alleles based on index from sorted dict
               #snps_vcf.at[i,"ALT"] = ",".join([ snp_alt_alleles[int(i-1)] for i in sorted_counts if i != 0 ])
               #print("DEBUG : "  + str(snps_vcf.shape[0]))
          
          # Get counts for ATGC bases
          snps = snp_alt.split(',')
          cnts = [ int(i) for i in counts.split(',') ] 
          nt_dict = dict(zip(snps, cnts))
          atcg_count = [ float(nt_dict.get(k, 0)) for k in nt ]
          atcg_cov = get_cov(num_samples, sum(atcg_count))
          # Update flag
          if (atcg_cov >= float(cov)) or (snp_pos in clade_nuc) : c_flag = 'YES'
          # Calculate frequency of alt alleles      
          freq = ",".join([ str(round(int(i) / int(num_samples),3)) if i != 0 else 0  for i in counts.split(',')])
          alt_list.append(snp_alt)
          aln_counts.append(counts)
          alele_freq.append(freq)
          atcg_cov_thr.append(atcg_cov)
          conf_flag.append(c_flag)
       else:
          drop_snps.append(int(snp_pos))

    print("INFO: Dropping SNPS : " + ",".join([ str(i) for i in drop_snps]))
    

    snps_vcf = snps_vcf[~snps_vcf['POS'].isin(drop_snps)]
    snps_vcf['ALT'] = alt_list
    print("Number of filtered rows : " +  str(snps_vcf.shape[0])) 
    out_vcf = snps_vcf.iloc[:,0:8]
    out_vcf['TOTAL_SAMPLES']= num_samples   
    out_vcf['OCCURENCES'] = aln_counts
    out_vcf['ALELE_FREQ'] = alele_freq
    out_vcf['COV'] = atcg_cov_thr
    out_vcf['CONF_FLAG'] = conf_flag

    out_vcf = out_vcf[out_vcf.OCCURENCES != '']
    print("Filtered SNPs : " +str(out_vcf.shape[0]) + " ; Samples : " + "\t" + str(snps_vcf.shape[1]-9) )
    out_vcf.to_csv(out_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n', mode = 'a')
