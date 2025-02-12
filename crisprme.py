#!/usr/bin/env python

from Bio.Seq import Seq

import subprocess
import itertools
import sys
import os
import re


version = "2.1.7"  #  CRISPRme version; TODO: update when required
__version__ = version

script_path = os.path.dirname(os.path.abspath(__file__))
origin_path = os.path.dirname(os.path.abspath(__file__))
# path where this file is located
# origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crisprme/PostProcess/"
# path corrected to use with conda
corrected_origin_path = script_path[:-3] + conda_path
corrected_web_path = f"{origin_path[:-3]}opt/crisprme/"
# corrected_web_path = os.getcwd()

script_path = corrected_origin_path
current_working_directory = f"{os.getcwd()}/"
# script_path = corrected_web_path+"/PostProcess/"

input_args = sys.argv

if "--debug" in input_args:
    print("DEBUG MODE")
    script_path = current_working_directory + "PostProcess/"
    corrected_web_path = current_working_directory

VALID_CHARS = {
    "a",
    "A",
    "t",
    "T",
    "c",
    "C",
    "g",
    "G",
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
    "v",
}


# Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    name = "_".join(name.split())
    current_working_directory = os.getcwd() + "/"
    chrom = input_range.split(":")[0]
    start_position = (
        input_range.split(":")[1]
        .split("-")[0]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )
    end_position = (
        input_range.split(":")[1]
        .split("-")[1]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )

    list_chr = [
        f
        for f in os.listdir(current_working_directory + "Genomes/" + genome_selected)
        if os.path.isfile(
            os.path.join(current_working_directory + "Genomes/" + genome_selected, f)
        )
        and not f.endswith(".fai")
    ]
    add_ext = ".fa"
    if ".fasta" in list_chr[0]:
        add_ext = ".fasta"
    with open(current_working_directory + name + ".bed", "w") as b:
        b.write(chrom + "\t" + start_position + "\t" + end_position)

    output_extract = subprocess.check_output(
        [
            "bedtools getfasta -fi "
            + current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + add_ext
            + " -bed "
            + current_working_directory
            + name
            + ".bed"
        ],
        shell=True,
    ).decode("utf-8")
    try:
        os.remove(
            current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + ".fa.fai"
        )
    except:
        pass
    try:
        os.remove(current_working_directory + name + ".bed")
    except:
        pass
    ret_string = output_extract.split("\n")[1].strip()
    return ret_string


def getGuides(extracted_seq, pam, len_guide, pam_begin):
    len_pam = len(pam)
    # dict
    len_guide = int(len_guide)
    pam_dict = {
        "A": "ARWMDHV",
        "C": "CYSMBHV",
        "G": "GRSKBDV",
        "T": "TYWKBDH",
        "R": "ARWMDHVSKBG",
        "Y": "CYSMBHVWKDT",
        "S": "CYSMBHVKDRG",
        "W": "ARWMDHVYKBT",
        "K": "GRSKBDVYWHT",
        "M": "ARWMDHVYSBC",
        "B": "CYSMBHVRKDGWT",
        "D": "ARWMDHVSKBGYT",
        "H": "ARWMDHVYSBCKT",
        "V": "ARWMDHVYSBCKG",
        "N": "ACGTRYSWKMBDHV",
    }
    list_prod = []
    for char in pam:
        list_prod.append(pam_dict[char])

    iupac_pam = []  # NNNNNNN NGG
    for element in itertools.product(*list_prod):
        iupac_pam.append("".join(element))

    rev_pam = str(Seq(pam).reverse_complement())
    list_prod = []
    for char in rev_pam:
        list_prod.append(pam_dict[char])

    # CCN NNNNNNN  -> results found with this pam must be reverse complemented
    iupac_pam_reverse = []
    for element in itertools.product(*list_prod):
        iupac_pam_reverse.append("".join(element))

    extracted_seq = extracted_seq.upper()
    len_sequence = len(extracted_seq)
    guides = []
    for pam in iupac_pam:
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    guides.append(extracted_seq[i + len_pam : i + len_pam + len_guide])
                else:
                    if i < len_guide:
                        continue
                    # guides.append(extracted_seq[i-len_guide:i+len_pam])           # i is position where first char of pam is found, eg the N char in NNNNNN NGG
                    # print('1 for:' , extracted_seq[i-len_guide:i])
                    guides.append(extracted_seq[i - len_guide : i])
    for pam in iupac_pam_reverse:  # Negative strand
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i < len_guide:
                        continue
                    guides.append(
                        str(Seq(extracted_seq[i - len_guide : i]).reverse_complement())
                    )
                else:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    # guides.append(str(Seq(extracted_seq[i:i+len_pam+len_guide]).reverse_complement()))         # i is position where first char of pam is found, eg the first C char in CCN NNNNNN
                    # print('2 for:', str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
                    guides.append(
                        str(
                            Seq(
                                extracted_seq[i + len_pam : i + len_guide + len_pam]
                            ).reverse_complement()
                        )
                    )
    return guides
    # return guides for when adding to app.py


def directoryCheck():
    # function to check the main directory status, if some directory is missing, create it
    directoryList = [
        "Genomes",
        "Results",
        "Dictionaries",
        "VCFs",
        "Annotations",
        "PAMs",
        "samplesIDs",
    ]
    for directory in directoryList:
        if not os.path.exists(current_working_directory + directory):
            os.makedirs(current_working_directory + directory)


def complete_search():
    variant = True
    if "--help" in input_args:
        print(
            "This is the automated search process that goes from raw input up to the post-analysis of results."
        )
        print("These are the flags that must be used in order to run this function:")
        print("\t--genome, used to specify the reference genome folder")
        print(
            "\t--vcf, used to specify the file containing a list of VCF folders (one per line) [OPTIONAL!]"
        )
        print(
            "\t--guide, used to specify the file that contains guides used for the search [IF NOT --sequence]"
        )
        print(
            "\t--sequence, used to specify the file containing DNA sequences or bed coordinates to extract guides [IF NOT --guide]"
        )
        print("\t--pam, used to specify the file that contains the pam")
        print(
            "\t--be-window, used to specify the window to search for susceptibilty to certain base editor (e.g., --be-window 4,8)"
        )
        print(
            "\t--be-base, used to specify the base(s) to check for the choosen editor (e.g., --be-base A,C)"
        )
        print(
            "\t--annotation, used to specify the file that contains annotations of the reference genome"
        )
        print(
            "\t--personal_annotation, used to specify the file that contains personal annotations of the reference genome"
        )
        print(
            "\t--samplesID, used to specify the file with a list of files (one per line) containing the information about samples present in VCF files [OPTIONAL!]"
        )
        print(
            "\t--gene_annotation, used to specify a gencode or similar annotation to find nearest gene for each target found [OPTIONAL]"
        )
        print(
            "\t--mm, used to specify the number of mismatches permitted in the search phase"
        )
        print(
            "\t--bDNA, used to specify the number of DNA bulges permitted in the search phase [OPTIONAL!]"
        )
        print(
            "\t--bRNA, used to specify the number of RNA bulges permitted in the search phase [OPTIONAL!]"
        )
        print(
            "\t--merge, used to specify the window (# of nucleotides) within which to merge candidate off-targets, using the off-target with the highest score as the pivot [default 3]"
        )
        print(
            "\t--sorting-criteria-scoring, specify target sorting criteria using a comma-separated list: 'mm' for mismatches, 'bulges' for bulges, or 'mm+bulges' for both (scoring has highest priority) [default 'mm+bulges']"
        )
        print(
            "\t--sorting-criteria, specify target sorting criteria using a comma-separated list: 'mm' for mismatches, 'bulges' for bulges, or 'mm+bulges' for both [default 'mm+bulges,mm']"
        )
        print(
            "\t--output, used to specify the output name for the results (these results will be saved into Results/<name>)"
        )
        print(
            "\t--thread, used to set the number of thread used in the process [default 8]"
        )
        exit(0)

    # check if all directories are found, if not, create them
    directoryCheck()

    # check for base and window in base editor
    if "--be-window" in input_args and "--be-base" not in input_args:
        print("Please input the base(s) editor to check in specified window")
        exit(1)
    if "--be-base" in input_args and "--be-window" not in input_args:
        print("Please input the base window to check for the specified base")
        exit(1)

    # check guide and sequence existence
    if "--guide" not in input_args and "--sequence" not in input_args:
        print("Please input a guide file or a sequence file")
        exit(1)
    if "--guide" in input_args and "--sequence" in input_args:
        print("Please select only ONE input type, either --guide or --sequence")
        exit(1)

    # base editor input check
    base_start = 1
    base_end = 0
    base_set = "none"
    if "--be-window" in input_args:
        try:
            base_window = input_args[input_args.index("--be-window") + 1]
            try:
                base_start = int(base_window.strip().split(",")[0])
                base_end = int(base_window.strip().split(",")[1])
            except:
                print("Please input a valid set of numbers for flag --be-window")
                exit(1)
        except IndexError:
            print("Please input some parameter for flag --be-window")
            exit(1)
    if "--be-base" in input_args:
        try:
            base_set = input_args[input_args.index("--be-base") + 1]
            for base in base_set.strip().split(","):
                if base not in VALID_CHARS:
                    print("Please input a set of valid nucleotides (A,C,G,T)")
                    exit(1)
        except IndexError:
            print("Please input some parameter for flag --be-base")
            exit(1)

    # guide input check
    if "--guide" in input_args:
        try:
            guidefile = os.path.abspath(input_args[input_args.index("--guide") + 1])
        except IndexError:
            print("Please input some parameter for flag --guide")
            exit(1)
        if not os.path.isfile(guidefile):
            print("The file specified for --guide does not exist")
            exit(1)

    # sequence input check
    sequence_use = False
    if "--sequence" in input_args:
        try:
            sequence_file = os.path.abspath(
                input_args[input_args.index("--sequence") + 1]
            )
            sequence_use = True
        except IndexError:
            print("Please input some parameter for flag --sequence")
            exit(1)
        if not os.path.isfile(sequence_file):
            print("The file specified for --sequence does not exist")
            exit(1)

    # check input genome
    if "--genome" not in input_args:
        print("--genome must be contained in the input")
        exit(1)
    else:
        try:
            genomedir = os.path.abspath(input_args[input_args.index("--genome") + 1])
        except IndexError:
            print("Please input some parameter for flag --genome")
            exit(1)
        if not os.path.isdir(genomedir):
            print("The folder specified for --genome does not exist")
            exit(1)

    # check input thread
    if "--thread" not in input_args:
        thread = 8  # set to avoid errors in following procedures
    else:
        try:
            thread = input_args[input_args.index("--thread") + 1]
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

    # check input vcf
    if "--vcf" not in input_args:
        variant = False
        vcfdir = "_"
    else:
        try:
            vcfdir = os.path.realpath(input_args[input_args.index("--vcf") + 1])
        except IndexError:
            print("Please input some parameter for flag --vcf")
            exit(1)
        if not os.path.isfile(vcfdir):
            print("The file specified for --vcf does not exist")
            exit(1)

    # check input gene-annotation
    if "--gene_annotation" not in input_args:
        gene_annotation = script_path + "vuoto.txt"
    else:
        try:
            gene_annotation = os.path.abspath(
                input_args[input_args.index("--gene_annotation") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --gene_annotation")
            exit(1)
        if not os.path.isfile(gene_annotation):
            print("The file specified for --gene_annotation does not exist")
            exit(1)

    # check input pam
    if "--pam" not in input_args:
        print("--pam must be contained in the input")
        exit(1)
    else:
        try:
            pamfile = os.path.abspath(input_args[input_args.index("--pam") + 1])
        except IndexError:
            print("Please input some parameter for flag --pam")
            exit(1)
        if not os.path.isfile(pamfile):
            print("The file specified for --pam does not exist")
            exit(1)

    # check input functional annotation
    if "--annotation" not in input_args:
        print("--annotation not used")
        annotationfile = script_path + "vuoto.txt"
        # exit(1)
    else:
        try:
            annotationfile = os.path.abspath(
                input_args[input_args.index("--annotation") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --annotation")
            exit(1)
        if not os.path.isfile(annotationfile):
            print("The file specified for --annotation does not exist")
            exit(1)
        if "--personal_annotation" in input_args:
            try:
                personal_annotation_file = os.path.abspath(
                    input_args[input_args.index("--personal_annotation") + 1]
                )
            except:
                pass
            if not os.path.isfile(personal_annotation_file):
                print("The file specified for --personal_annotation does not exist")
                exit(1)
            os.system(
                f'awk \'$4 = $4"_personal"\' {personal_annotation_file} | sed "s/ /\t/g" | sed "s/,/_personal,/g" > {personal_annotation_file}.tmp'
            )
            os.system(
                f"cat {personal_annotation_file}.tmp {annotationfile} > {annotationfile}+personal.bed"
            )
            os.system(f"rm -f {personal_annotation_file}.tmp")
            annotationfile = annotationfile + "+personal.bed"

    # check input personal annotation
    if "--personal_annotation" in input_args and "--annotation" not in input_args:
        try:
            personal_annotation_file = os.path.abspath(
                input_args[input_args.index("--personal_annotation") + 1]
            )
        except:
            pass
        if not os.path.isfile(personal_annotation_file):
            print("The file specified for --personal_annotation does not exist")
            exit(1)
        os.system(
            f'awk \'$4 = $4"_personal"\' {personal_annotation_file} | sed "s/ /\t/g" | sed "s/,/_personal,/g" > {personal_annotation_file}.tmp'
        )
        os.system(
            f"cat {personal_annotation_file}.tmp {annotationfile} > {annotationfile}+personal.bed"
        )
        os.system(f"rm -f {personal_annotation_file}.tmp")
        annotationfile = annotationfile + "+personal.bed"

    # check input for variant search (existance of all necessary file)
    samplefile = (
        script_path + "vuoto.txt"
    )  # use void file for samples if variant not used
    if variant and "--samplesID" not in input_args:
        print("--samplesID must be contained in the input to perform variant search")
        exit(1)
    elif not variant and "--samplesID" in input_args:
        print("--samplesID was in the input but no VCF directory was specified")
        exit(1)
    elif "--samplesID" in input_args:
        try:
            samplefile = os.path.abspath(
                input_args[input_args.index("--samplesID") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --samplesID")
            exit(1)
        if not os.path.isfile(samplefile):
            print("The file specified for --samplesID does not exist")
            exit(1)

    # check input bMax
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

    # check input mm
    if "--mm" not in input_args:
        print("--mm must be contained in the input")
        exit(1)
    else:
        try:
            mm = input_args[input_args.index("--mm") + 1]
        except IndexError:
            print("Please input some parameter for flag --mm")
            exit(1)
        try:
            mm = int(mm)
        except:
            print("Please input a number for flag mm")
            exit(1)

    # check input bDNA
    if "--bDNA" not in input_args:
        # print("--bDNA must be contained in the input")
        # exit(1)
        bDNA = 0
    else:
        try:
            bDNA = input_args[input_args.index("--bDNA") + 1]
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

    # check input bRNA
    if "--bRNA" not in input_args:
        # print("--bRNA must be contained in the input")
        # exit(1)
        bRNA = 0
    else:
        try:
            bRNA = input_args[input_args.index("--bRNA") + 1]
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

    # set bMAX to generate index as max value (bDNA,bRNA)
    bMax = max(bDNA, bRNA)

    # check input merge window
    if "--merge" not in input_args:
        merge_t = 3  # default merge is 3 nt
    else:
        try:
            merge_t = input_args[input_args.index("--merge") + 1]
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

    # check sorting criteria while merging targets on score
    if "--sorting-criteria-scoring" not in input_args:
        sorting_criteria_scoring = "mm+bulges"
    else:
        try:
            sorting_criteria_scoring = input_args[
                input_args.index("--sorting-criteria-scoring") + 1
            ]
        except IndexError as e:
            sys.stderr.write(
                "Please input some parameter for flag --sorting-criteria-scoring\n"
            )
            exit(1)
        if len(sorting_criteria_scoring.split(",")) > len(
            set(sorting_criteria_scoring.split(","))
        ):
            sys.stderr.write("Repeated sorting criteria\n")
            exit(1)
        if len(sorting_criteria_scoring.split(",")) > 3:
            sys.stderr.write("Forbidden or repeated sorting criteria\n")
            exit(1)
        if any(
            c not in ["mm+bulges", "mm", "bulges"]
            for c in sorting_criteria_scoring.split(",")
        ):
            sys.stderr.write("Forbidden sorting criteria selected\n")
            exit(1)

    # check sorting criteria while  merging targets (fewest mm+bulges)
    if "--sorting-criteria" not in input_args:
        sorting_criteria = "mm+bulges,mm"
    else:
        try:
            sorting_criteria = input_args[input_args.index("--sorting-criteria") + 1]
        except IndexError as e:
            sys.stderr.write(
                "Please input some parameter for flag --sorting-criteria\n"
            )
            exit(1)
        if len(sorting_criteria.split(",")) > len(set(sorting_criteria.split(","))):
            sys.stderr.write("Repeated sorting criteria\n")
            exit(1)
        if len(sorting_criteria.split(",")) > 3:
            sys.stderr.write("Forbidden or repeated sorting criteria\n")
            exit(1)
        if any(
            c not in ["mm+bulges", "mm", "bulges"] for c in sorting_criteria.split(",")
        ):
            sys.stderr.write("Forbidden sorting criteria selected\n")
            exit(1)

    # check input output directory
    if "--output" not in input_args:
        print("--output must be contained in the input")
        exit(1)
    else:
        try:
            outputfolder = (
                current_working_directory
                + "Results/"
                + input_args[input_args.index("--output") + 1]
            )
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

    # extract pam seq from file
    pam_len = 0
    total_pam_len = 0
    with open(pamfile, "r") as pam_file:
        pam_char = pam_file.readline()
        total_pam_len = len(pam_char.split(" ")[0])
        index_pam_value = pam_char.split(" ")[-1]
        if int(pam_char.split(" ")[-1]) < 0:
            end_idx = int(pam_char.split(" ")[-1]) * (-1)
            pam_char = pam_char.split(" ")[0][0:end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(" ")[-1])
            pam_char = pam_char.split(" ")[0][end_idx * (-1) :]
            pam_len = end_idx
            pam_begin = False

    genome_ref = os.path.basename(genomedir)
    annotation_name = os.path.basename(annotationfile)
    nuclease = os.path.basename(pamfile).split(".")[0].split("-")[2]
    if bMax != 0:
        search_index = True
    else:
        search_index = False
    if variant:
        genome_idx_list = []
        with open(vcfdir, "r") as vcfs:
            for line in vcfs:
                if line.strip():
                    if line[-2] == "/":
                        line = line[:-2]
                    base_vcf = os.path.basename(line)
                    genome_idx_list.append(
                        pam_char
                        + "_"
                        + str(bMax)
                        + "_"
                        + genome_ref
                        + "+"
                        + base_vcf.strip()
                    )
        genome_idx = ",".join(genome_idx_list)
        ref_comparison = True
    else:
        genome_idx = pam_char + "_" + str(bMax) + "_" + genome_ref
        ref_comparison = False
    # os.chdir(script_path)
    # write crisprme version to file
    with open(outputfolder + "/.command_line.txt", "w") as p:
        p.write("input_command\t" + " ".join(sys.argv[:]))
        p.write("\n")
        p.close()
    with open(outputfolder + "/.version.txt", "w") as p:
        p.write("crisprme_version\t" + __version__)
        p.write("\n")
        p.close()
    # write parameters to file
    with open(outputfolder + "/Params.txt", "w") as p:
        p.write("Genome_selected\t" + genome_ref.replace(" ", "_") + "\n")
        p.write("Genome_ref\t" + genome_ref + "\n")
        if search_index:
            p.write("Genome_idx\t" + genome_idx + "\n")
        else:
            p.write("Genome_idx\t" + "None\n")
        p.write("Pam\t" + pam_char + "\n")
        p.write("Max_bulges\t" + str(bMax) + "\n")
        p.write("Mismatches\t" + str(mm) + "\n")
        p.write("DNA\t" + str(bDNA) + "\n")
        p.write("RNA\t" + str(bRNA) + "\n")
        p.write("Annotation\t" + str(annotation_name) + "\n")
        p.write("Nuclease\t" + str(nuclease) + "\n")
        # p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write("Ref_comp\t" + str(ref_comparison) + "\n")
        p.close()
    len_guide_sequence = total_pam_len - pam_len
    if sequence_use:
        guides = list()
        text_sequence = str()
        for line in open(sequence_file, "r"):
            text_sequence += line
        for name_and_seq in text_sequence.split(">"):
            if "" == name_and_seq:
                continue
            name = name_and_seq[: name_and_seq.find("\n")]
            seq = name_and_seq[name_and_seq.find("\n") :]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if "chr" in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split("\n"):
                    if "" == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = (
                        pieces_of_row[0]
                        + ":"
                        + pieces_of_row[1]
                        + "-"
                        + pieces_of_row[2]
                    )
                    extracted_seq = extractSequence(
                        name, seq_to_extract, genome_ref.replace(" ", "_")
                    )
                    guides.extend(
                        getGuides(
                            extracted_seq, pam_char, len_guide_sequence, pam_begin
                        )
                    )
            else:
                seq = seq.split()
                seq = "".join(seq)
                extracted_seq = seq.strip()
                guides.extend(
                    getGuides(extracted_seq, pam_char, len_guide_sequence, pam_begin)
                )
        temp_guides = list()
        for guide in guides:
            addN = "N" * pam_len
            if pam_begin:
                temp_guides.append(addN + guide)
            else:
                temp_guides.append(guide + addN)
        if len(temp_guides) > 1000000000:
            temp_guides = temp_guides[:1000000000]
        guides = temp_guides
        extracted_guides_file = open(outputfolder + "/guides.txt", "w")
        for guide in guides:
            extracted_guides_file.write(guide + "\n")
        extracted_guides_file.close()
    # print(guides)
    # exit(0)
    void_mail = "_"
    if sequence_use == False:
        os.system(f"cp {guidefile} {outputfolder}/guides.txt")
    print(
        f"Launching job {outputfolder}. The stdout is redirected in log_verbose.txt and stderr is redirected in log_error.txt"
    )
    # start search with set parameters
    with open(f"{outputfolder}/log_verbose.txt", "w") as log_verbose:
        with open(f"{outputfolder}/log_error.txt", "w") as log_error:
            crisprme_run = (
                f"{os.path.join(script_path, 'submit_job_automated_new_multiple_vcfs.sh')} "
                f"{genomedir} {vcfdir} {os.path.join(outputfolder, 'guides.txt')} "
                f"{pamfile} {annotationfile} {samplefile} {bMax} {mm} {bDNA} {bRNA} "
                f"{merge_t} {outputfolder} {script_path} {thread} {current_working_directory} "
                f"{gene_annotation} {void_mail} {base_start} {base_end} {base_set} "
                f"{sorting_criteria_scoring} {sorting_criteria}"
            )
            code = subprocess.call(
                crisprme_run, shell=True, stderr=log_error, stdout=log_verbose
            )
            if code != 0:
                raise OSError(
                    f"\nCRISPRme run failed! See {os.path.join(outputfolder, 'log_error.txt')} for details\n"
                )
            # subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), str(vcfdir), str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(
            #     samplefile), str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
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
        print(
            "This is the automated integration process that process the final result file to generate a usable target panel."
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--targets, used to specify the final result file to use in the panel creation process"
        )
        print(
            "\t--empirical_data, used to specify the file that contains empirical data provided by the user to assess in-silico targets"
        )
        print("\t--output, used to specify the output folder for the results")
        exit(0)

    if "--targets" not in input_args:
        print("--targets must be contained in the input")
        exit(1)
    else:
        try:
            target_file = os.path.abspath(input_args[input_args.index("--targets") + 1])
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
        empiricalfile = script_path + "vuoto.txt"
        # exit(1)
    else:
        try:
            empiricalfile = os.path.abspath(
                input_args[input_args.index("--empirical_data") + 1]
            )
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
            outputfolder = os.path.abspath(input_args[input_args.index("--output") + 1])
        except IndexError:
            print("Please input some parameter for flag --output")
            exit(1)
        if not os.path.isdir(outputfolder):
            print("The folder specified for --output does not exist")
            exit(1)

    os.system(
        f"{script_path}./empirical_integrator.py {target_file} {empiricalfile} {outputfolder}"
    )


def print_help_gnomad_converter():
    """
    Prints the help information for the gnomAD converter functionality, providing
    details on the conversion process from gnomAD VCFs to VCFs compatible with
    CRISPRme. It outlines the options available for specifying directories,
    sample IDs, variant filtering, multiallelic site handling, and thread usage
    during the conversion process.

    Raises:
        SystemExit: If the help information is displayed to guide users on using
        the gnomAD converter functionality.
    """

    # functionality description
    sys.stderr.write(
        "The gnomAD converter functionality simplifies the conversion process "
        "of gnomAD VCFs (versions 3.1 and 4.0) into VCFs supported by CRISPRme. "
        "It ensures a seamless transition while maintaining compatibility with "
        "CRISPRme's requirements, focusing on the structure and content of "
        "precomputed sample IDs file \n\n"
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--gnomAD_VCFdir, specifies the directory containing gnomAD VCFs. "
        "Files must have the BGZ extension\n"
        "\t--samplesID, specifies the precomputed sample IDs file necessary "
        "for incorporating population-specific information into the output "
        "VCFs\n"
        "\t--joint, optional flag to specify the input GnomAD VCF contain joint "
        "allele frequencies\n"
        "\t--keep, optional flag to retain all variants, regardless of their "
        "filter flag. By default, variants with a filter flag different from "
        "PASS are discarded\n"
        "\t--multiallelic, optional flag to merge variants mapped to the "
        "same position, creating multiallelic sites in the output VCFs. By "
        "default, each site remains biallelic\n"
        "\t--thread, used to set the number of thread used in the conversion "
        "process [default 8]\n"
    )
    sys.exit(1)


def gnomAD_converter():
    """
    Runs the gnomAD converter functionality based on specified arguments, converting
    gnomAD VCF files into formats compatible with CRISPRme.

    Raises:
        ValueError: If mandatory arguments are missing or have incorrect values.
        FileExistsError: If the specified gnomAD VCF directory cannot be located.
        FileNotFoundError: If the specified sample IDs file cannot be found.
        subprocess.SubprocessError: If an error occurs during the gnomAD VCF
            conversion process.
    """

    args = input_args[2:]  # reover gnomAD converter args
    if "--help" in args or not args:  # print help
        print_help_gnomad_converter()
        sys.exit(1)
    if "--gnomAD_VCFdir" not in args:
        raise ValueError(
            "--gnomAD_VCFdir is a mandatory argument required for the conversion "
            "process. Please specify the directory containing gnomAD VCFs using "
            "this option\n"
        )
    if "--samplesID" not in args:
        raise ValueError(
            "--samplesID is a mandatory argument required for the conversion "
            "process. Please specify the sample IDs file this option\n"
        )
    # read gnomAD directory arg
    try:
        gnomad_dir = args[args.index("--gnomAD_VCFdir") + 1]
        if gnomad_dir.startswith("--"):
            raise ValueError("Please input some parameter for flag --gnomAD_VCFdir\n")
        gnomad_dir = os.path.abspath(gnomad_dir)  # first sanity check passed
        if not os.path.isdir(gnomad_dir):
            raise FileExistsError(f"Unable to locate {gnomad_dir}")
    except IndexError as e:
        raise ValueError(
            "Please input some parameter for flag --gnomAD_VCFdir\n"
        ) from e
    # read samples ids arg
    try:
        samples_ids = args[args.index("--samplesID") + 1]
        if samples_ids.startswith("--"):
            raise ValueError("Please input some parameter for flag --samplesID")
        samples_ids = os.path.abspath(samples_ids)  # first sanity check passed
        if not os.path.isfile(samples_ids):
            raise FileNotFoundError(f"Unable to locate {samples_ids}")
    except IndexError as e:
        raise ValueError("Please input some parameter for flag --samplesID") from e
    # read joint gnomad vcf files
    joint = "--joint" in args
    # read keep arg
    keep = "--keep" in args  # keep all variants regardless of filter label
    # read multiallelic arg
    multiallelic = "--multiallelic" in args  # merge variants in multiallelic sites
    # read threads arg
    threads = 8
    if "--threads" in args:
        try:
            threads = int(args[args.index("--threads") + 1])
            if threads <= 0:
                raise ValueError(f"Forbidden number of threads ({threads})")
        except IndexError as e:
            raise ValueError("Missing or forbidden threads value") from e
    # run gnom AD converter
    gnomad_converter_script = os.path.join(script_path, "convert_gnomAD_vcfs.py")
    cmd = (
        f"python {gnomad_converter_script} {gnomad_dir} {samples_ids} {joint} "
        f"{keep} {multiallelic} {threads}"
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise subprocess.SubprocessError(
            f"An error occurred while converting gnomAD VCFs in {gnomad_dir}"
        )


def personal_card():
    if "--help" in input_args:
        print(
            "This is the personal card generator that creates a files with all the private targets for the input sample"
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--result_dir, directory containing the result from which extract the targets to generate the card"
        )
        print(
            "\t--guide_seq, sequence of the guide to use in order to exctract the targets"
        )
        print("\t--sample_id, ID of the sample to use in order to generate the card")
        exit(0)

    if "--result_dir" not in input_args:
        print("--result_dir not in input, please input a result directory")
        exit(1)
    else:
        try:
            result_dir = os.path.abspath(
                input_args[input_args.index("--result_dir") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --result_dir")
            exit(1)
        if not os.path.isdir(result_dir):
            print("The folder specified for --result_dir does not exist")
            exit(1)

    if "--guide_seq" not in input_args:
        print(
            "--guide_seq must be contained in the input, e.g. CTAACAGTTGCTTTTATCACNNN"
        )
        exit(1)
    else:
        try:
            guide = input_args[input_args.index("--guide_seq") + 1]
        except IndexError:
            print("Please input some parameter for flag --guide_seq")
            exit(1)
    if "--sample_id" not in input_args:
        print("--sample_id must be contained in the input, e.g. HG00001")
        exit(1)
    else:
        try:
            sample_id = input_args[input_args.index("--sample_id") + 1]
        except IndexError:
            print("Please input some parameter for flag --sample_id")
            exit(1)

    os.system(
        script_path
        + "./generate_sample_card.py "
        + result_dir
        + " "
        + guide
        + " "
        + sample_id
        + " "
        + script_path
    )


def web_interface():
    if "--help" in input_args:
        print(
            "This function must be launched without input, it starts a local server to use the web interface."
        )
        print(
            "Open your web-browser and write 127.0.0.1:8080 in the search bar if you are executing locally, if you are executing on an external server write <yourserverip>:8080 in search bar"
        )
        exit(0)
    subprocess.run(corrected_web_path + "/./index.py")


def crisprme_version():
    if len(input_args) != 2:
        sys.stderr.write("Wrong number of arguments for crisprme.py version\n")
        sys.exit(1)
    sys.stdout.write(f"v{__version__}\n")


def print_help_complete_test():
    """
    Prints the help information for executing comprehensive testing of the
    complete-search functionality provided by CRISPRme.

    Raises:
        SystemExit: If the help information is displayed to guide users on
            executing comprehensive testing.
    """

    # write intro message to stdout
    sys.stderr.write(
        "Execute comprehensive testing for complete-search functionality "
        "provided by CRISPRme\n"
    )
    # list functionality options
    sys.stderr.write(
        "Options:\n"
        "\t--chrom, test the complete-search functionality on the specified "
        "chromosome (e.g., chr22). By default, the test is conducted on all "
        "chromosomes\n"
        "\t--vcf_dataset, VCFs dataset to be used during CRISPRme testing. "
        "Available options include 1000 Genomes (1000G) and Human Genome "
        "Diversity Project (HGDP). To use the combined dataset type '1000G+HGDP' "
        "The default dataset is 1000 Genomes.\n"
        "\t--thread, number of threads.\n"
        "\t--debug, debug mode.\n"
    )
    sys.exit(1)


def complete_test_crisprme():
    """
    Executes comprehensive testing for the complete-search functionality provided
    by CRISPRme based on specified arguments.

    Raises:
        OSError: If the CRISPRme test fails, indicating an issue with the testing
            process.
    """

    if "--help" in input_args or len(input_args) < 3:
        print_help_complete_test()
        sys.exit(1)
    chrom = "all"
    if "--chrom" in input_args:  # individual chrom to test
        try:
            chrom = input_args[input_args.index("--chrom") + 1]
            if chrom.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --chrom\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --chrom\n")
            sys.exit(1)
    vcf_dataset = "1000G"
    if "--vcf_dataset" in input_args:  # specified variant dataset
        try:
            vcf_dataset = input_args[input_args.index("--vcf_dataset") + 1]
            if vcf_dataset.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
            sys.exit(1)
    threads = 4
    if "--thread" in input_args:  # number of threads to use during test
        try:
            threads = input_args[input_args.index("--thread") + 1]
            if threads.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --thread\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some value for flag --thread\n")
            sys.exit(1)
    debug = "--debug" in input_args  # run local or via conda/Docker
    # begin crisprme test
    script_test = os.path.join(script_path, "complete_test.py")
    code = subprocess.call(
        f"python {script_test} {chrom} {vcf_dataset} {threads} {debug}", shell=True
    )
    if code != 0:
        raise OSError(
            "\nCRISPRme test failed! See Results/crisprme-test-out/log_error.txt for details\n"
        )


# HELP FUNCTION
def callHelp():
    # print general help, listing all available functions with a brief
    # description of their purpose
    sys.stderr.write(
        "Help:\n\n"
        "- ALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND SEPARATED BY CHROMOSOME\n"
        "- ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED (WITH BGZIP) AND SEPARATED BY CHROMOSOME\n\n"
        "Functionalities:\n\n"
        "crisprme.py complete-search\n"
        "\tPerforms genome-wide off-targets search (reference and variant, if "
        "specified), including CFD and CRISTA analysis, and target selection\n\n"
        "crisprme.py complete-test\n"
        "\tTest the complete CRISPRme pipeline on single chromosomes or complete "
        "genomes\n\n"
        "crisprme.py targets-integration\n"
        "\tIntegrates in-silico targets with empirical data to generate a usable "
        "panel\n\n"
        "crisprme.py gnomAD-converter\n"
        "\tConverts gnomAD VCF files into CRISPRme compatible VCFs (supports "
        "gnomAD >= v3.1)\n\n"
        "crisprme.py generate-personal-card\n"
        "\tGenerates a personal card for specific samples by extracting all "
        "private targets\n\n"
        "crisprme.py web-interface\n"
        "\tActivates CRISPRme's web interface for local browser use\n\n"
        "crisprme.py --version\n"
        "\tPrints CRISPRme version to stdout and exit\n\n"
        "For additional information on each CRISPRme functionality type <function> "
        "--help (e.g. 'crisprme.py complete-search --help')\n"
    )
    sys.exit(1)  # stop execution


if len(sys.argv) < 2:
    directoryCheck()
    callHelp()
elif sys.argv[1] == "complete-search":
    complete_search()
elif sys.argv[1] == "complete-test":
    complete_test_crisprme()
elif sys.argv[1] == "targets-integration":
    target_integration()
elif sys.argv[1] == "gnomAD-converter":
    gnomAD_converter()
elif sys.argv[1] == "generate-personal-card":
    personal_card()
elif sys.argv[1] == "web-interface":
    web_interface()
elif sys.argv[1] == "--version":
    crisprme_version()
else:
    sys.stderr.write('ERROR! "' + sys.argv[1] + '" is not an allowed!\n\n')
    callHelp()  # print help when wrong input command