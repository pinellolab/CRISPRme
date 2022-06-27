#!/usr/bin/env python

import sys
import os
import subprocess
import itertools
from Bio.Seq import Seq
import re
# for getting lst of chr to know file extension and if enriched
from os.path import isfile, isdir, join
from os import listdir

script_path = os.path.dirname(os.path.abspath(__file__))
origin_path = os.path.dirname(os.path.abspath(__file__))
# path where this file is located
# origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crisprme/PostProcess/"
# path corrected to use with conda
corrected_origin_path = script_path[:-3]+conda_path
corrected_web_path = origin_path[:-3]+"opt/crisprme/"
# corrected_web_path = os.getcwd()

script_path = corrected_origin_path
current_working_directory = os.getcwd() + '/'
# script_path = corrected_web_path+"/PostProcess/"

input_args = sys.argv

if '--debug' in input_args:
    print('DEBUG MODE')
    script_path = current_working_directory+'PostProcess/'
    corrected_web_path = current_working_directory

VALID_CHARS = {'a', 'A', 't', 'T', 'c', 'C', 'g', 'G',
               "R",
               "Y",
               "S",
               "W",
               "K",
               "M",
               "B",
               "D",
               "H",
               "V",
               "r",
               "y",
               "s",
               "w",
               "k",
               "m",
               "b",
               "d",
               "h",
               "v"
               }


# Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    name = '_'.join(name.split())
    current_working_directory = os.getcwd() + '/'
    chrom = input_range.split(':')[0]
    start_position = input_range.split(':')[1].split(
        '-')[0].replace(',', '').replace('.', '').replace(' ', '')
    end_position = input_range.split(':')[1].split(
        '-')[1].replace(',', '').replace('.', '').replace(' ', '')

    list_chr = [f for f in listdir(current_working_directory + 'Genomes/' + genome_selected) if isfile(
        join(current_working_directory + 'Genomes/' + genome_selected, f)) and not f.endswith('.fai')]
    add_ext = '.fa'
    if '.fasta' in list_chr[0]:
        add_ext = '.fasta'
    with open(current_working_directory + name + '.bed', 'w') as b:
        b.write(chrom + '\t' + start_position + '\t' + end_position)

    output_extract = subprocess.check_output(['bedtools getfasta -fi ' + current_working_directory + 'Genomes/' + genome_selected +
                                              '/' + chrom + add_ext + ' -bed ' + current_working_directory + name + '.bed'], shell=True).decode("utf-8")
    try:
        os.remove(current_working_directory + 'Genomes/' +
                  genome_selected + '/' + chrom + '.fa.fai')
    except:
        pass
    try:
        os.remove(current_working_directory + name + '.bed')
    except:
        pass
    ret_string = output_extract.split('\n')[1].strip()
    return ret_string


def getGuides(extracted_seq, pam, len_guide, pam_begin):
    len_pam = len(pam)
    # dict
    len_guide = int(len_guide)
    pam_dict = {
        'A':  "ARWMDHV",
        'C':  "CYSMBHV",
        'G':  "GRSKBDV",
        'T':  "TYWKBDH",
        'R':  "ARWMDHVSKBG",
        'Y':  "CYSMBHVWKDT",
        'S':  "CYSMBHVKDRG",
        'W':  "ARWMDHVYKBT",
        'K':  "GRSKBDVYWHT",
        'M':  "ARWMDHVYSBC",
        'B':  "CYSMBHVRKDGWT",
        'D':  "ARWMDHVSKBGYT",
        'H':  "ARWMDHVYSBCKT",
        'V':  "ARWMDHVYSBCKG",
        'N':  "ACGTRYSWKMBDHV",
    }
    list_prod = []
    for char in pam:
        list_prod.append(pam_dict[char])

    iupac_pam = []  # NNNNNNN NGG
    for element in itertools.product(*list_prod):
        iupac_pam.append(''.join(element))

    rev_pam = str(Seq(pam).reverse_complement())
    list_prod = []
    for char in rev_pam:
        list_prod.append(pam_dict[char])

    # CCN NNNNNNN  -> results found with this pam must be reverse complemented
    iupac_pam_reverse = []
    for element in itertools.product(*list_prod):
        iupac_pam_reverse.append(''.join(element))

    extracted_seq = extracted_seq.upper()
    len_sequence = len(extracted_seq)
    guides = []
    for pam in iupac_pam:
        pos = ([m.start()
                for m in re.finditer('(?=' + pam + ')', extracted_seq)])
        if pos:
            for i in pos:
                if pam_begin:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    guides.append(
                        extracted_seq[i + len_pam: i + len_pam + len_guide])
                else:
                    if i < len_guide:
                        continue
                    # guides.append(extracted_seq[i-len_guide:i+len_pam])           # i is position where first char of pam is found, eg the N char in NNNNNN NGG
                    # print('1 for:' , extracted_seq[i-len_guide:i])
                    guides.append(extracted_seq[i-len_guide:i])
    for pam in iupac_pam_reverse:  # Negative strand
        pos = ([m.start()
                for m in re.finditer('(?=' + pam + ')', extracted_seq)])
        if pos:
            for i in pos:
                if pam_begin:
                    if i < len_guide:
                        continue
                    guides.append(
                        str(Seq(extracted_seq[i-len_guide:i]).reverse_complement()))
                else:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    # guides.append(str(Seq(extracted_seq[i:i+len_pam+len_guide]).reverse_complement()))         # i is position where first char of pam is found, eg the first C char in CCN NNNNNN
                    # print('2 for:', str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
                    guides.append(str(
                        Seq(extracted_seq[i + len_pam: i + len_guide + len_pam]).reverse_complement()))
    return guides
    # return guides for when adding to app.py


def directoryCheck():
    # function to check the main directory status, if some directory is missing, create it
    directoryList = ['Genomes', 'Results', 'Dictionaries',
                     'VCFs', 'Annotations', 'PAMs', 'samplesIDs']
    for directory in directoryList:
        if not os.path.exists(current_working_directory+directory):
            os.makedirs(current_working_directory+directory)


def complete_search():
    variant = True
    if "--help" in input_args:
        print("This is the automated search process that goes from raw input up to the post-analysis of results.")
        print("These are the flags that must be used in order to run this function:")
        print("\t--genome, used to specify the reference genome folder")
        print(
            "\t--vcf, used to specify the file containing a list of VCF folders (one per line) [OPTIONAL!]")
        print(
            "\t--guide, used to specify the file that contains guides used for the search [IF NOT --sequence]")
        print(
            "\t--sequence, used to specify the file containing DNA sequences or bed coordinates to extract guides [IF NOT --guide]")
        print("\t--pam, used to specify the file that contains the pam")
        print("\t--be-window, used to specify the window to search for susceptibilty to certain base editor (e.g., --be-window 4,8)")
        print("\t--be-base, used to specify the base(s) to check for the choosen editor (e.g., --be-base A,C)")
        print("\t--annotation, used to specify the file that contains annotations of the reference genome")
        print("\t--personal_annotation, used to specify the file that contains personal annotations of the reference genome")
        print(
            "\t--samplesID, used to specify the file with a list of files (one per line) containing the information about samples present in VCF files [OPTIONAL!]")
        print(
            "\t--gene_annotation, used to specify a gencode or similar annotation to find nearest gene for each target found [OPTIONAL]")
        print(
            "\t--mm, used to specify the number of mismatches permitted in the search phase")
        print(
            "\t--bDNA, used to specify the number of DNA bulges permitted in the search phase [OPTIONAL!]")
        print(
            "\t--bRNA, used to specify the number of RNA bulges permitted in the search phase [OPTIONAL!]")
        print("\t--output, used to specify the output name for the results (these results will be saved into Results/<name>)")
        print("\t--thread, used to set the number of thread used in the process (default is 8)")
        exit(0)

    # check if all directories are found, if not, create them
    directoryCheck()

    #check for base and window in base editor
    if '--be-window' in input_args and '--be-base' not in input_args:
        print('Please input the base(s) editor to check in specified window')
        exit(1)    
    if '--be-base' in input_args and '--be-window' not in input_args:
        print('Please input the base window to check for the specified base')
        exit(1)    
    
    #check guide and sequence existence
    if '--guide' not in input_args and '--sequence' not in input_args:
        print('Please input a guide file or a sequence file')
        exit(1)
    if '--guide' in input_args and '--sequence' in input_args:
        print('Please select only ONE input type, either --guide or --sequence')
        exit(1)
        
    # base editor input check
    base_start=1
    base_end=0
    base_set="none"
    if '--be-window' in input_args:
        try:
            base_window = input_args[input_args.index("--be-window")+1]
            try:
                base_start=int(base_window.strip().split(',')[0])
                base_end=int(base_window.strip().split(',')[1])
            except:
                print("Please input a valid set of numbers for flag --be-window")
                exit(1)
        except IndexError:
            print("Please input some parameter for flag --be-window")
            exit(1)
    if '--be-base' in input_args:
        try:
            base_set = input_args[input_args.index("--be-base")+1]
            for base in base_set.strip().split(','):
                if base not in VALID_CHARS:
                    print('Please input a set of valid nucleotides (A,C,G,T)')
                    exit(1)
        except IndexError:
            print("Please input some parameter for flag --be-base")
            exit(1)
            
    #guide input check
    if "--guide" in input_args:
        try:
            guidefile = os.path.abspath(
                input_args[input_args.index("--guide")+1])
        except IndexError:
            print("Please input some parameter for flag --guide")
            exit(1)
        if not os.path.isfile(guidefile):
            print("The file specified for --guide does not exist")
            exit(1)
    
    # sequence input check
    sequence_use = False
    if '--sequence' in input_args:
        try:
            sequence_file = os.path.abspath(
                input_args[input_args.index("--sequence")+1])
            sequence_use = True
        except IndexError:
            print("Please input some parameter for flag --sequence")
            exit(1)
        if not os.path.isfile(sequence_file):
            print("The file specified for --sequence does not exist")
            exit(1)

    #check input genome
    if "--genome" not in input_args:
        print("--genome must be contained in the input")
        exit(1)
    else:
        try:
            genomedir = os.path.abspath(
                input_args[input_args.index("--genome")+1])
        except IndexError:
            print("Please input some parameter for flag --genome")
            exit(1)
        if not os.path.isdir(genomedir):
            print("The folder specified for --genome does not exist")
            exit(1)

    #check input thread
    if "--thread" not in input_args:
        thread = 8  # set to avoid errors in following procedures
    else:
        try:
            thread = input_args[input_args.index("--thread")+1]
        except IndexError:
            print("Please input some parameter for flag --thread")
            exit(1)
        try:
            thread = int(thread)
        except:
            print("Please input a number for flag --thread")
            exit(1)
        if thread <= 0:
            print("thread is set to default (8) ")
            thread = 8
    
    #check input vcf
    if "--vcf" not in input_args:
        variant = False
        vcfdir='_'
    else:
        try:
            vcfdir = os.path.realpath(input_args[input_args.index("--vcf")+1])
        except IndexError:
            print("Please input some parameter for flag --vcf")
            exit(1)
        if not os.path.isfile(vcfdir):
            print("The file specified for --vcf does not exist")
            exit(1)

     #check input gene-annotation
    if "--gene_annotation" not in input_args:
        gene_annotation = script_path+'vuoto.txt'
    else:
        try:
            gene_annotation = os.path.abspath(
                input_args[input_args.index("--gene_annotation")+1])
        except IndexError:
            print("Please input some parameter for flag --gene_annotation")
            exit(1)
        if not os.path.isfile(gene_annotation):
            print("The file specified for --gene_annotation does not exist")
            exit(1)

     #check input pam
    if "--pam" not in input_args:
        print("--pam must be contained in the input")
        exit(1)
    else:
        try:
            pamfile = os.path.abspath(input_args[input_args.index("--pam")+1])
        except IndexError:
            print("Please input some parameter for flag --pam")
            exit(1)
        if not os.path.isfile(pamfile):
            print("The file specified for --pam does not exist")
            exit(1)

     #check input functional annotation
    if "--annotation" not in input_args:
        print("--annotation not used")
        annotationfile = script_path+'vuoto.txt'
        # exit(1)
    else:
        try:
            annotationfile = os.path.abspath(
                input_args[input_args.index("--annotation")+1])
        except IndexError:
            print("Please input some parameter for flag --annotation")
            exit(1)
        if not os.path.isfile(annotationfile):
            print("The file specified for --annotation does not exist")
            exit(1)
        if '--personal_annotation' in input_args:
            try:
                personal_annotation_file = os.path.abspath(
                    input_args[input_args.index("--personal_annotation")+1])
            except:
                pass
            if not os.path.isfile(personal_annotation_file):
                print("The file specified for --personal_annotation does not exist")
                exit(1)
            os.system(
                f'awk \'$4 = $4\"_personal\"\' {personal_annotation_file} | sed "s/ /\t/g" | sed "s/,/_personal,/g" > {personal_annotation_file}.tmp')
            os.system(
                f'cat {personal_annotation_file}.tmp {annotationfile} > {annotationfile}+personal.bed')
            os.system(f'rm -f {personal_annotation_file}.tmp')
            annotationfile = annotationfile+'+personal.bed'

     #check input personal annotation
    if '--personal_annotation' in input_args and '--annotation' not in input_args:
        try:
            personal_annotation_file = os.path.abspath(
                input_args[input_args.index("--personal_annotation")+1])
        except:
            pass
        if not os.path.isfile(personal_annotation_file):
            print("The file specified for --personal_annotation does not exist")
            exit(1)
        os.system(
            f'awk \'$4 = $4\"_personal\"\' {personal_annotation_file} | sed "s/ /\t/g" | sed "s/,/_personal,/g" > {personal_annotation_file}.tmp')
        os.system(
            f'cat {personal_annotation_file}.tmp {annotationfile} > {annotationfile}+personal.bed')
        os.system(f'rm -f {personal_annotation_file}.tmp')
        annotationfile = annotationfile+'+personal.bed'

    #check input for variant search (existance of all necessary file)
    samplefile=script_path+'vuoto.txt' #use void file for samples if variant not used
    if variant and "--samplesID" not in input_args:
        print("--samplesID must be contained in the input to perform variant search")
        exit(1)
    elif not variant and "--samplesID" in input_args:
        print("--samplesID was in the input but no VCF directory was specified")
        exit(1)
    elif "--samplesID" in input_args:
        try:
            samplefile = os.path.abspath(
                input_args[input_args.index("--samplesID")+1])
        except IndexError:
            print("Please input some parameter for flag --samplesID")
            exit(1)
        if not os.path.isfile(samplefile):
            print("The file specified for --samplesID does not exist")
            exit(1)
            
    #check input bMax
    # if "--bMax" not in input_args:
    #     print("--bMax must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         bMax = input_args[input_args.index("--bMax")+1]
    #     except IndexError:
    #         print("Please input some parameter for flag --bMax")
    #         exit(1)
    #     try:
    #         bMax = int(bMax)
    #     except:
    #         print("Please input a number for flag bMax")
    #         exit(1)
    #     # if bMax < 0 or bMax > 2:
    #     #     print("The range for bMax is from 0 to 2")
    #     #     exit(1)

    #check input mm
    if "--mm" not in input_args:
        print("--mm must be contained in the input")
        exit(1)
    else:
        try:
            mm = input_args[input_args.index("--mm")+1]
        except IndexError:
            print("Please input some parameter for flag --mm")
            exit(1)
        try:
            mm = int(mm)
        except:
            print("Please input a number for flag mm")
            exit(1)

    #check input bDNA
    if "--bDNA" not in input_args:
        # print("--bDNA must be contained in the input")
        # exit(1)
        bDNA = 0
    else:
        try:
            bDNA = input_args[input_args.index("--bDNA")+1]
        except IndexError:
            print("Please input some parameter for flag --bDNA")
            exit(1)
        try:
            bDNA = int(bDNA)
        except:
            print("Please input an integer number for flag --bDNA")
            exit(1)
        # if bDNA > bMax:
        #     print("The number of bDNA must be equal or less than bMax")
        #     exit(1)
        # elif bDNA < 0 or bDNA > bMax:
        #     print("The range for bDNA is from 0 to", bMax)
        #     exit(1)

    #check input bRNA
    if "--bRNA" not in input_args:
        # print("--bRNA must be contained in the input")
        # exit(1)
        bRNA = 0
    else:
        try:
            bRNA = input_args[input_args.index("--bRNA")+1]
        except IndexError:
            print("Please input some parameter for flag --bRNA")
            exit(1)
        try:
            bRNA = int(bRNA)
        except:
            print("Please input an integer number for flag --bRNA")
            exit(1)
        # if bRNA > bMax:
        #     print("The number of bRNA must be equal or less than bMax")
        #     exit(1)
        # elif bRNA < 0 or bRNA > 2:
        #     print("The range for bRNA is from 0 to", bMax)
        #     exit(1)
    
    #set bMAX to generate index as max value (bDNA,bRNA)
    bMax=max(bDNA,bRNA)
    
    #check input merge window
    if "--merge" not in input_args:
        merge_t = 3  # default merge is 3 nt
    else:
        try:
            merge_t = input_args[input_args.index("--merge")+1]
        except IndexError:
            print("Please input some parameter for flag --merge")
            exit(1)
        try:
            merge_t = int(merge_t)
        except:
            print("Please input a number for flag merge")
            exit(1)
        if merge_t < 0:
            print("Please specify a positive number for --merge")
            exit(1)

    #check input output directory
    if "--output" not in input_args:
        print("--output must be contained in the input")
        exit(1)
    else:
        try:
            outputfolder = current_working_directory+'Results/' + \
                input_args[input_args.index("--output")+1]
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
            # outputfolder = os.path.abspath(
            #     input_args[input_args.index("--output")+1])
        except IndexError:
            print("Please input some parameter for flag --output")
            exit(1)
        if not os.path.isdir(outputfolder):
            print("The folder specified for --output does not exist")
            exit(1)

    #extract pam seq from file
    pam_len = 0
    total_pam_len = 0
    with open(pamfile, 'r') as pam_file:
        pam_char = pam_file.readline()
        total_pam_len = len(pam_char.split(' ')[0])
        index_pam_value = pam_char.split(' ')[-1]
        if int(pam_char.split(' ')[-1]) < 0:
            end_idx = int(pam_char.split(' ')[-1]) * (-1)
            pam_char = pam_char.split(' ')[0][0: end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(' ')[-1])
            pam_char = pam_char.split(' ')[0][end_idx * (-1):]
            pam_len = end_idx
            pam_begin = False

    genome_ref = os.path.basename(genomedir)
    annotation_name = os.path.basename(annotationfile)
    nuclease = os.path.basename(pamfile).split('.')[0].split('-')[2]
    if bMax != 0:
        search_index = True
    else:
        search_index = False
    if variant:
        genome_idx_list = []
        with open(vcfdir, 'r') as vcfs:
            for line in vcfs:
                if line.strip():
                    if line[-2] == "/":
                        line = line[:-2]
                    base_vcf = os.path.basename(line)
                    genome_idx_list.append(
                        pam_char + '_' + str(bMax) + '_' + genome_ref + '+' + base_vcf.strip())
        genome_idx = ','.join(genome_idx_list)
        ref_comparison = True
    else:
        genome_idx = pam_char + '_' + str(bMax) + '_' + genome_ref
        ref_comparison = False
    # os.chdir(script_path)
    with open(outputfolder + '/Params.txt', 'w') as p:
        p.write('Genome_selected\t' + genome_ref.replace(' ', '_') + '\n')
        p.write('Genome_ref\t' + genome_ref + '\n')
        if search_index:
            p.write('Genome_idx\t' + genome_idx + '\n')
        else:
            p.write('Genome_idx\t' + 'None\n')
        p.write('Pam\t' + pam_char + '\n')
        p.write('Max_bulges\t' + str(bMax) + '\n')
        p.write('Mismatches\t' + str(mm) + '\n')
        p.write('DNA\t' + str(bDNA) + '\n')
        p.write('RNA\t' + str(bRNA) + '\n')
        p.write('Annotation\t' + str(annotation_name) + '\n')
        p.write('Nuclease\t' + str(nuclease) + '\n')
        # p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write('Ref_comp\t' + str(ref_comparison) + '\n')
        p.close()
    len_guide_sequence = total_pam_len - pam_len
    if sequence_use:
        guides = list()
        text_sequence = str()
        for line in open(sequence_file, 'r'):
            text_sequence += line
        for name_and_seq in text_sequence.split('>'):
            if '' == name_and_seq:
                continue
            name = name_and_seq[:name_and_seq.find('\n')]
            seq = name_and_seq[name_and_seq.find('\n'):]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if 'chr' in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split('\n'):
                    if '' == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = pieces_of_row[0]+":" + \
                        pieces_of_row[1]+"-"+pieces_of_row[2]
                    extracted_seq = extractSequence(
                        name, seq_to_extract, genome_ref.replace(' ', '_'))
                    guides.extend(getGuides(
                        extracted_seq, pam_char, len_guide_sequence, pam_begin))
            else:
                seq = seq.split()
                seq = ''.join(seq)
                extracted_seq = seq.strip()
                guides.extend(getGuides(
                    extracted_seq, pam_char, len_guide_sequence, pam_begin))
        temp_guides = list()
        for guide in guides:
            addN = 'N'*pam_len
            if pam_begin:
                temp_guides.append(addN+guide)
            else:
                temp_guides.append(guide+addN)
        if len(temp_guides) > 1000000000:
            temp_guides = temp_guides[:1000000000]
        guides = temp_guides
        extracted_guides_file = open(outputfolder+'/guides.txt', 'w')
        for guide in guides:
            extracted_guides_file.write(guide+'\n')
        extracted_guides_file.close()
    # print(guides)
    # exit(0)
    void_mail='_'
    if sequence_use == False:
        os.system(f'cp {guidefile} {outputfolder}/guides.txt')
    print(
        f"Launching job {outputfolder}. The stdout is redirected in log_verbose.txt and stderr is redirected in log_error.txt")
    # start search with set parameters
    with open(f"{outputfolder}/log_verbose.txt", 'w') as log_verbose:
        with open(f"{outputfolder}/log_error.txt", 'w') as log_error:
            subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), str(vcfdir), str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(
                samplefile), str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # else:
    #     with open(f"{outputfolder}/log_verbose.txt", 'w') as log_verbose:
    #         with open(f"{outputfolder}/log_error.txt", 'w') as log_error:
    #             subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), '_', str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(script_path+'vuoto.txt'),
    #                             str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # change name of guide and param files to hidden
    os.system(f"mv {outputfolder}/guides.txt {outputfolder}/.guides.txt")
    os.system(f"mv {outputfolder}/Params.txt {outputfolder}/.Params.txt")


def target_integration():
    if "--help" in input_args:
        print("This is the automated integration process that process the final result file to generate a usable target panel.")
        print("These are the flags that must be used in order to run this function:")
        print("\t--targets, used to specify the final result file to use in the panel creation process")
        print("\t--empirical_data, used to specify the file that contains empirical data provided by the user to assess in-silico targets"),
        print("\t--output, used to specify the output folder for the results")
        exit(0)

    if "--targets" not in input_args:
        print("--targets must be contained in the input")
        exit(1)
    else:
        try:
            target_file = os.path.abspath(
                input_args[input_args.index("--targets")+1])
        except IndexError:
            print("Please input some parameter for flag --targets")
            exit(1)
        if not os.path.isfile(target_file):
            print("The file specified for --target_file does not exist")
            exit(1)

    # if "--vcf_dir" not in input_args:
    #     print("--vcf_dir non in input, multi-variant haplotype will not be calculated")
    #     vcf_dir = script_path+'vuota/'
    #     # exit(1)
    # else:
    #     try:
    #         vcf_dir = os.path.abspath(
    #             input_args[input_args.index("--vcf_dir")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --vcf_dir")
    #         exit(1)
    #     if not os.path.isdir(vcf_dir):
    #         print("The folder specified for --vcf_dir does not exist")
    #         exit(1)

    # if "--genome_version" not in input_args:
    #     print("--genome_version must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         genome_version = input_args[input_args.index(
    #             "--genome_version")+1]
    #     except IndexError:
    #         print("Please input some parameter for flag --genome")
    #         exit(1)

    # if "--guide" not in input_args:
    #     guidefile = script_path+'vuoto.txt'
    #     # print("--guide must be contained in the input")
    #     # exit(1)
    # else:
    #     try:
    #         guidefile = os.path.abspath(
    #             input_args[input_args.index("--guide")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --guide")
    #         exit(1)
    #     if not os.path.isfile(guidefile):
    #         print("The file specified for --guide does not exist")
    #         exit(1)

    if "--empirical_data" not in input_args:
        print("--empirical_data not in input, proceeding without empirical data")
        empiricalfile = script_path+'vuoto.txt'
        # exit(1)
    else:
        try:
            empiricalfile = os.path.abspath(
                input_args[input_args.index("--empirical_data")+1])
        except IndexError:
            print("Please input some parameter for flag --empirical_data")
            exit(1)
        if not os.path.isfile(empiricalfile):
            print("The file specified for --empirical_data does not exist")
            exit(1)

    # if "--gencode" not in input_args:
    #     print("--gencode must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         gencode_file = os.path.abspath(
    #             input_args[input_args.index("--gencode")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --gencode")
    #         exit(1)
    #     if not os.path.isfile(gencode_file):
    #         print("The file specified for --gencode does not exist")
    #         exit(1)

    if "--output" not in input_args:
        print("--output must be contained in the input")
        exit(1)
    else:
        try:
            outputfolder = os.path.abspath(
                input_args[input_args.index("--output")+1])
        except IndexError:
            print("Please input some parameter for flag --output")
            exit(1)
        if not os.path.isdir(outputfolder):
            print("The folder specified for --output does not exist")
            exit(1)

    os.system(
        f'{script_path}./empirical_integrator.py {target_file} {empiricalfile} {outputfolder}')


def gnomAD_converter():
    if "--help" in input_args:
        print("This is the VCF gnomAD converter provided to convert all gnomADv3.1 VCFs into CRISPRme supported VCFs")
        print("These are the flags that must be used in order to run this function:")
        print("\t--gnomAD_VCFdir, used to specify the directory containing gnomADv3.1 original VCFs")
        print("\t--samplesID, used to specify the pre-generated samplesID file necessary to introduce samples into gnomAD variant")
        print(
            "\t--thread, used to specify the number of core used to process VCFs in parallel (DEFAULT is ALL available minus 2) [OPTIONAL]")
        exit(0)

    if "--gnomAD_VCFdir" not in input_args:
        print("--gnomAD_VCFdir not in input, MANDATORY TO CONVERT DATA")
        # vcf_dir = script_path+'vuota/'
        exit(1)
    else:
        try:
            vcf_dir = os.path.abspath(
                input_args[input_args.index("--gnomAD_VCFdir")+1])
        except IndexError:
            print("Please input some parameter for flag --gnomAD_VCFdir")
            exit(1)
        if not os.path.isdir(vcf_dir):
            print("The folder specified for --gnomAD_VCFdir does not exist")
            exit(1)

    if "--thread" not in input_args:
        # print("--thread must be contained in the input")
        # exit(1)
        thread = len(os.sched_getaffinity(0))-2
    else:
        try:
            thread = input_args[input_args.index("--thread")+1]
        except IndexError:
            print("Please input some parameter for flag --thread")
            exit(1)
        try:
            thread = int(thread)
        except:
            print("Please input a number for flag thread")
            exit(1)
        if thread <= 0 or thread > len(os.sched_getaffinity(0))-2:
            print("thread is set to default (ALL available minus 2)")
            thread = len(os.sched_getaffinity(0))-2
            # exit(1)

    if "--samplesID" not in input_args:
        print("--samplesID not in input, MANDATORY TO CONVERT DATA")
        exit(1)
    elif "--samplesID" in input_args:
        try:
            samplefile = os.path.abspath(
                input_args[input_args.index("--samplesID")+1])
        except IndexError:
            print("Please input some parameter for flag --samplesID")
            exit(1)
        if not os.path.isfile(samplefile):
            print("The file specified for --samplesID does not exist")
            exit(1)

    os.system(script_path+"./convert_gnomAD.py " +
              vcf_dir+" "+samplefile + " " + str(thread))


def personal_card():
    if "--help" in input_args:
        print("This is the personal card generator that creates a files with all the private targets for the input sample")
        print("These are the flags that must be used in order to run this function:")
        print("\t--result_dir, directory containing the result from which extract the targets to generate the card")
        print(
            "\t--guide_seq, sequence of the guide to use in order to exctract the targets")
        print("\t--sample_id, ID of the sample to use in order to generate the card")
        exit(0)

    if "--result_dir" not in input_args:
        print("--result_dir not in input, please input a result directory")
        exit(1)
    else:
        try:
            result_dir = os.path.abspath(
                input_args[input_args.index("--result_dir")+1])
        except IndexError:
            print("Please input some parameter for flag --result_dir")
            exit(1)
        if not os.path.isdir(result_dir):
            print("The folder specified for --result_dir does not exist")
            exit(1)

    if "--guide_seq" not in input_args:
        print("--guide_seq must be contained in the input, e.g. CTAACAGTTGCTTTTATCACNNN")
        exit(1)
    else:
        try:
            guide = input_args[input_args.index(
                "--guide_seq")+1]
        except IndexError:
            print("Please input some parameter for flag --guide_seq")
            exit(1)
    if "--sample_id" not in input_args:
        print("--sample_id must be contained in the input, e.g. HG00001")
        exit(1)
    else:
        try:
            sample_id = input_args[input_args.index(
                "--sample_id")+1]
        except IndexError:
            print("Please input some parameter for flag --sample_id")
            exit(1)

    os.system(script_path+"./generate_sample_card.py "+result_dir+" "+guide +
              " "+sample_id+" "+script_path)


def web_interface():
    if "--help" in input_args:
        print("This function must be launched without input, it starts a local server to use the web interface.")
        print("Open your web-browser and write 127.0.0.1:8080 in the search bar if you are executing locally, if you are executing on an external server write <yourserverip>:8080 in search bar")
        exit(0)
    subprocess.run(corrected_web_path+'/./index.py')


# HELP FUNCTION
def callHelp():
    print("help:\n",
          "\nALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND CHROMOSOME SEPARATED, ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED AND CHROMOSOME SEPARATED",
          "\n",
          "\ncrisprme.py complete-search FUNCTION SEARCHING THE WHOLE GENOME (REFERENCE AND VARIANT IF REQUESTED) AND PERFORM CFD ANALYSIS AND TARGET SELECTION",
          "\ncrisprme.py targets-integration FUNCTION THAT INTEGRATES IN-SILICO TARGETS WITH EMPIRICAL DATA GENERATING A USABLE PANEL",
          "\ncrisprme.py gnomAD-converter FUNCTION THAT CONVERTS ALL gnomADv3.1 vcf.bgz FILES INTO COMPATIBLE VCFs",
          "\ncrisprme.py generate-personal-card FUNCTION TO GENERATE PERSONAL CARD FOR A SPECIFIC SAMPLE EXTRACTING ALL THE PRIVATE TARGETS",
          "\ncrisprme.py web-interface FUNCTION TO ACTIVATE WEB INTERFACE OF CRISPRme TO USE WITH A BROWSER LOCALLY"
          "\n\nADD help TO ANY FUNCTION TO VISUALIZE A BRIEF HELP PAGE (example: crisprme.py complete-search --help)\n")


if len(sys.argv) < 2:
    directoryCheck()
    callHelp()
elif sys.argv[1] == 'complete-search':
    complete_search()
elif sys.argv[1] == 'gnomAD-converter':
    gnomAD_converter()
# elif sys.argv[1] == 'search-only':
#     search_only()
# elif sys.argv[1] == 'post-analysis-only':
#     post_analysis_only()
elif sys.argv[1] == 'targets-integration':
    target_integration()
elif sys.argv[1] == 'web-interface':
    web_interface()
elif sys.argv[1] == 'generate-personal-card':
    personal_card()
else:
    print("ERROR! \"" + sys.argv[1] + "\" is not an allowed!")
