#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 21:38:03 2020

@author: francesco
"""

# 30 min per vcf di chr1 -> 6gb
import gzip
import sys
import json
import time
import os

# argv1 is ALLchrx.gzip --> Create the json file containing chr, pos, ref, alt, list of samples (HG001,HG002...)
# argv2 is chr number (eg 1)
vcf = sys.argv[1]
out = sys.argv[2]
start_time = time.time()
with open(out+".vcf.indels_only", 'w') as out:
    with gzip.open(vcf, 'rt', encoding='utf-8') as targets: #gzip.   'rt', encoding='utf-8'
        #Skip vcf header
        for line in targets:
            #line = line.decode('ascii')
            if ('#CHROM') in line:
                column_vcf = line.strip().split('\t')   #Save this header for retrieving sample id
                break
        out.write('\t'.join(column_vcf)+'\n')
        for i, line in enumerate(targets):                #Save CHROM [0], POS[1], RSID[2], REF [3], ALT [4], AF[7], List of Samples [9:]
            line = line.strip().split('\t') #.decode('ascii')
            
            if len(line[3]) != len(line[4]) and '<' not in line[3] and '<' not in line[4]:
                if ',' in line[4]:
                    continue
                    #splitted = line[4].split(',')
                    #filtered_splitted = [item for item in splitted if len(item) != len(line[3])]
                    #line[4] = ','.join(filtered_splitted)
                    
                out.write('\t'.join(line)+'\n')
       
print('Filtered VCF in ', time.time() - start_time)
