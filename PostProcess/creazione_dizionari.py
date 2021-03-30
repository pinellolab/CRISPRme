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
with gzip.open(sys.argv[1], 'rb') as targets:
    #Skip vcf header
    for line in targets:
        line = line.decode('ascii')
        if ('#CHROM') in line:
            column_vcf = line.strip().split('\t')   #Save this header for retrieving sample id
            break

    for line in targets:                #Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        line = line.decode('ascii').strip().split('\t')
        list_samples = []
        list_chars = []
        for pos, i in enumerate(line[9:]):          #if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
            if ('1' in i):
                list_samples.append(column_vcf[ pos + 9])
        if "chr" not in line[0]:
            line[0] = "chr"+line[0]
        chr_pos_string = line[0] + ',' + line[1]
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
        #If no sample was found, the dict is chrX,100 -> ;A,T
        list_chars.append(line[3])
        list_chars.append(line[4])
        try:
            chr_dict[chr_pos_string] = ','.join(sorted(list_samples)) + ';' + ','.join(list_chars)
        except:
            chr_dict[chr_pos_string] = ';' + ','.join(list_chars) #None
        
with open(os.path.dirname(sys.argv[2]) + "/my_dict_" + str(line[0]) + '.json', 'w') as f:
    json.dump(chr_dict, f) 
print('Created ' + sys.argv[2] + '.json' + ' in', time.time() - start_time)

