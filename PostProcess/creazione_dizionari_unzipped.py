#!/usr/bin/env python

'''
Generate json dictionaries for sample extraction
'''
# 30 min per vcf di chr1 -> 6gb
import gzip
import sys
import json
import time
import os

# argv1 is ALLchrx.gzip --> Create the json file containing chr, pos, ref, alt, list of samples (HG001,HG002...)
# argv2 is chr number (eg 1)
chr_dict = dict()
start_time = time.time()
with open(sys.argv[1], 'r') as targets: #gzip.open(sys.argv[1], 'rb')
    #Skip vcf header
    for line in targets:
        #line = line.decode('ascii')
        if ('#CHROM') in line:
            column_vcf = line.strip().split('\t')   #Save this header for retrieving sample id
            break
    first_line = True
    pos_AF = 0
    for line in targets:                #Save CHROM [0], POS[1], RSID[2], REF [3], ALT [4], AF[7], List of Samples [9:]
        line = line.strip().split('\t') #.decode('ascii')
        if first_line:
            first_line = False
            splitted = line[7].split(";")
            for pos, ele in enumerate(splitted):
                if ele[0:2] == "AF":
                    pos_AF = pos
                    break
        list_samples = []
        list_chars = []
        if len(line[3]) == 1 and len(line[4]) == 1:
            for pos, i in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
                if ('1' in i):
                    list_samples.append(column_vcf[ pos + 9])
            if "chr" not in line[0]:
                line[0] = "chr"+line[0]
            chr_pos_string = line[0] + ',' + line[1]
            #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
            #If no sample was found, the dict is chrX,100 -> ;A,T
            rsID = line[2]
            list_chars.append(line[3])
            list_chars.append(line[4])
            af = line[7].split(";")[pos_AF][3:] #retrieve AF=number --> number 
            if len(list_samples) > 0:
                chr_dict[chr_pos_string] = ','.join(sorted(list_samples)) + ';' + ','.join(list_chars) + ";" + rsID + ";" + af
            else:
                chr_dict[chr_pos_string] = ';' + ','.join(list_chars) + ";" + rsID + ";" + af #None
        elif len(line[3]) == 1:
            variants = line[4].split(",")
            snps = []
            for var in variants:
                if len(var) == 1: 
                    snps.append(var)
            if len(snps) > 0:		
                for pos, i in enumerate(line[9:]):
                    if ('1' in i):
                        list_samples.append(column_vcf[ pos + 9])
                if "chr" not in line[0]:
                    line[0] = "chr"+line[0]
                chr_pos_string = line[0] + ',' + line[1]
                rsID = line[2]
                af = line[7].split(";")[pos_AF][3:]
            for snp in snps:
                list_chars = [line[3]]
                list_chars.append(snp)
                if len(list_samples) > 0:
                    chr_dict[chr_pos_string] = ','.join(sorted(list_samples)) + ';' + ','.join(list_chars) + ";" + rsID + ";" + af
                else:
                    chr_dict[chr_pos_string] = ';' + ','.join(list_chars) + ";" + rsID + ";" + af #None
        

with open('my_dict_chr' + sys.argv[2] + '.json', 'w') as f:
    json.dump(chr_dict, f) 
print('Created ' + 'my_dict_chr' + sys.argv[2] + '.json' + ' in', time.time() - start_time)

