#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:39:09 2020

@author: francesco
"""
import os
import sys
import pandas as pd

#log = pd.read_csv(sys.argv[2], sep='\t')
#log["rsID"] = ['n'] * log.shape[0]
#log["AF"] = [0] * log.shape[0]
#log["indel"] = ['n'] * log.shape[0]

with open(sys.argv[1], 'r') as targets:
    with open(sys.argv[2], 'r') as log:
        with open(sys.argv[2]+".tmp", 'w') as logout:
            column_vcf = targets.readline()
            header = log.readline().strip()
            logout.write(header+"\trsID\tAF\tindel\n")
            first_line = True
            for i, line in enumerate(targets):                #Save CHROM [0], POS[1], RSID[2], REF [3], ALT [4], AF[7], List of Samples [9:]
                #print(i)
                line = line.strip().split('\t') 
                if first_line:
                    first_line = False
                    splitted = line[7].split(";")
                    for pos, ele in enumerate(splitted):
                        if ele[0:2] == "AF":
                            pos_AF = pos
                            #print(pos_AF)
                            break
                if "chr" not in line[0]:
                    line[0] = "chr"+line[0]
                chr_pos_string = line[0] + '_' + line[1] + '_' + line[3] + '_' +line[4]
                
                rsID = line[2]
                af = line[7].split(";")[pos_AF][3:] #retrieve AF=number --> number
                
                line_log = log.readline().strip()
                logout.write(line_log+'\t'+rsID+'\t'+af+'\t'+chr_pos_string+'\n')
                #log["AF"][i] = af
                #log["rsID"][i] = rsID
                #log["indel"][i] = chr_pos_string

os.system('mv '+sys.argv[2]+".tmp "+sys.argv[2])
#log.to_csv(sys.argv[2], sep='\t', index=False)
        
