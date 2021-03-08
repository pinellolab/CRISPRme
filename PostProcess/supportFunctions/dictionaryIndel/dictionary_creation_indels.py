'''
Creation of dictionaries for single sample indels. Takes in input a chromosome and a sample id, and selects the variants 
only from that sample. The position written in the dictionary is adjusted accordingly to the addition or removal of nucleotides.
The name of the fasta file and the name in the #CHROM column MUST be the same.
Use http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ as vcf reference
WARNING! About 20 minutes per vcf file, for a total of 80GB of disk space
'''
# argv 1 is All.chrX.gzip, the vcf file
# argv 2 is selected sample
# argv 3 is directory to save the dictionary

import gzip
import sys
import json
import time

vcf_file = sys.argv[1]
chr_name = vcf_file.split('.')
for i in chr_name:
    if 'chr' in i:
        chr_name = i
        break
selected_sample = sys.argv[2]
save_directory = sys.argv[3]

dictionary_dict = dict()
start_time = time.time()
add_to_name = ''     #string to add to the chr number, eg in VCF hg38 -> is already chr1, so add_to_name is '';
                        #VCF hg19 -> is 1, so add_to_name is 'chr'
offset = 0      #Sum this value to the position of the variant ad adjust this value for every indel found
with gzip.open(sys.argv[1], 'rb') as targets:
    for line in targets:
        line = line.decode('ascii')
        if ('#CHROM') in line:
            column_vcf = line.strip().split('\t')
            sample_header_pos = column_vcf.index(selected_sample)
            break

    first_line = next(targets).decode('ascii').strip().split('\t')
    if 'chr' not in first_line[0]:
        add_to_name = 'chr'
    list_samples = []
    list_chars = []
    
    hap = first_line[sample_header_pos].split('|')
    for h in hap:   #NOTE posso avere anche 2|1 #TODO da finire
        if '0' == h:
            continue
        variant_of_sample = int(h) - 1      #If in VAR i have TT,TTA,TTT and in sample 0|2, the variant is TTA
        list_samples.append(selected_sample)
        offset = offset + len(line[3]) - len(line[4].split(',')[variant_of_sample])
        if len(line[3]) != 1 or len(line[4].split(',')[variant_of_sample]) != 1:    #Save only the SNP
            break
        chr_pos_key = add_to_name + line[0] + ',' + str(offset + int(line[1]))
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
        #If no sample was found, the dict is chrX,100 -> ;A,T
        list_chars.append(line[3])  #REF char
        list_chars.append(line[4].split(',')[variant_of_sample])  #VAR char
        try:
            dictionary_dict[chr_pos_key] = ','.join(list_samples) + ';' + ','.join(list_chars)
        except:
            dictionary_dict[chr_pos_key] = ';' + ','.join(list_chars) #None
    
    for line in targets:                #Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        line = line.decode('ascii').strip().split('\t')
        list_samples = []
        list_chars = []
        if ('1' in first_line[sample_header_pos]):
            list_samples.append(selected_sample)
        else:           #This variant is not found in the selected sample
            continue
        chr_pos_key = add_to_name + line[0] + ',' + line[1]
        #Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
        #If no sample was found, the dict is chrX,100 -> ;A,T
        list_chars.append(line[3])
        list_chars.append(line[4])
        try:
            dictionary_dict[chr_pos_key] = ','.join(list_samples) + ';' + ','.join(list_chars)
        except:
            dictionary_dict[chr_pos_key] = ';' + ','.join(list_chars) #None
        #result.write(line[0] + '\t' + line[1] + '\t' + line[3] + '\t' + line[4] + '\t' + ','.join(list_samples) + '\n')
        
with open(save_directory + '_' + selected_sample + '/my_dict_' + chr_name + '.json', 'w') as f:
    json.dump(dictionary_dict, f) 
print('Created ' + 'my_dict_' + chr_name + '.json for sample ' + selected_sample + ' in', time.time() - start_time)