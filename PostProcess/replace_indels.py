#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 13:39:25 2020

@author: francesco
"""

import pandas as pd
import numpy as np
import time 
import re
import os
import sys

print("Substituting indels in fasta files")

start_time = time.time()
indels_pos = sys.argv[1] 
indel_samples = pd.read_csv(sys.argv[2], sep='\t', dtype=np.bool)
ref_alt = pd.read_csv(sys.argv[3], sep='\t')
folder_fastas = sys.argv[4]
chrom = sys.argv[5]

with open(folder_fastas+"/log"+chrom+".txt", 'w') as log:
    log.write("CHR\tSAMPLES\n")
#df = pd.read_csv(indels_file)
r = 0
with open(indels_pos, 'r') as in_file:
    for line in in_file:
        samples = indel_samples.columns[np.where(indel_samples.iloc[r,:] == True)]
        splitted = line.rstrip().split("\t")
        with open(folder_fastas+"/"+splitted[0]+"_"+splitted[1]+"-"+splitted[2]+"_"+str(r+1)+".fa", 'r') as fasta:
            first = fasta.readline().rstrip()
            second = fasta.readline().rstrip()
            
        if len(samples) > 0:
            first_out = first+"\t"+','.join(samples)
        else:
            first_out = first+"\t"+'NO_SAMPLES'
        second_out = second[0:25] + re.sub(ref_alt.iloc[r,2], ref_alt.iloc[r,3], second[25:], 1, flags=re.IGNORECASE) 
        
        with open(folder_fastas+"/"+splitted[0]+"_"+splitted[1]+"-"+splitted[2]+"_"+str(r+1)+".fa", 'w') as fasta:
            #fasta.write(first_out+"\n")
            fasta.write(first+"\n")
            fasta.write(second_out+"\n")
            
        with open(folder_fastas+"/log"+chrom+".txt", 'a') as log:
            log.write(first_out+"\n")
        
        r += 1

print("Done in "+str(time.time()-start_time))    
