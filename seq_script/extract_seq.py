#Use bedtool to extract sequence
import subprocess
import os
from os.path import isfile, isdir,join      #for getting lst of chr to know file extension and if enriched
from os import listdir

#Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    current_working_directory = os.getcwd() + '/'
    chrom = input_range.split(':')[0]
    start_position = input_range.split(':')[1].split('-')[0].replace(',','').replace('.','').replace(' ','')
    end_position = input_range.split(':')[1].split('-')[1].replace(',','').replace('.','').replace(' ','')

    list_chr = [f for f in listdir(current_working_directory + 'Genomes/' + genome_selected) if isfile(join(current_working_directory + 'Genomes/' + genome_selected, f)) and not f.endswith('.fai')]
    add_ext = '.fa'
    if '.fasta' in list_chr[0]:
        add_ext = '.fasta'
    with open(current_working_directory + name + '.bed','w') as b:
        b.write(chrom + '\t' + start_position + '\t' + end_position)

    output_extract = subprocess.check_output(['bedtools getfasta -fi ' + current_working_directory + 'Genomes/' + genome_selected + '/' + chrom + add_ext + ' -bed ' + current_working_directory + name + '.bed'], shell=True).decode("utf-8") 
    try:
        os.remove(current_working_directory + 'Genomes/' + genome_selected + '/' + chrom + '.fa.fai')
    except:
        pass
    try:
        os.remove(current_working_directory + name + '.bed')
    except:
        pass
    ret_string = output_extract.split('\n')[1].strip() 
    return ret_string