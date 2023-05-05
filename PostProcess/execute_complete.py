#!/usr/bin/env python
#
import sys
import datetime
import os
import subprocess
import shutil

# set -e # capture any failure

# file for automated search of guide+pam in reference and variant genomes
ref_folder = sys.argv[1]
vcf_list = open(sys.argv[2]).readlines()
guide_file = sys.argv[3]
pam_file = sys.argv[4]
annotation_file = sys.argv[5]
sampleID_list = sys.argv[6]
bMax = sys.argv[7]
mm = sys.argv[8]
bDNA = sys.argv[9]
bRNA = sys.argv[10]
merge_t = sys.argv[11]
output_folder = sys.argv[12]
processes_dir = sys.argv[13]
ncpus = sys.argv[14]
current_working_directory = sys.argv[15]
gene_proximity = sys.argv[16]
email = sys.argv[17]
# used to solve base editor check in resultintegration phase
base_check_start = sys.argv[18]
base_check_end = sys.argv[19]
base_check_set = sys.argv[20]

##MANDATORY PATHS
genomes_folder = os.path.join(current_working_directory, "Genomes")
genomes_libraries_folder = os.path.join(current_working_directory, "genome_library")
dictionaries_folder = os.path.join(current_working_directory, "Dictionaries")
results_folder = os.path.join(current_working_directory, "Results")
vcfs_folder = os.path.join(current_working_directory, "VCFs")

##GLOBAL VARIABLES
chr_list = list()
fake_chr_list = list()
vcf_list_checked = list()
vcf_process = ""
log_file = open(os.path.join(output_folder, "log.txt"), "w")
pam_complete = ""
pam_seq = ""
pam_position = 0
ref_name = os.path.basename(ref_folder).replace("/", "")


##USER FUNCTIONS
def write_to_log(message):
    ##write to log file with termination
    log_file.write(message + "\n")


def write_to_verbose(message):
    ##write to log verbose file with print autotermination
    print(message)


def write_to_error(message):
    ##write to log error file with print autotermination
    print(message, file=sys.stderr)


def pre_process():
    ## WRITE TO LOG_VERBOSE
    write_to_verbose("Starting pre-processing")
    write_to_verbose(f"input mail is: {email}")
    write_to_verbose(f"input ncpus is: {ncpus}")

    ##GENERATE LOG FILE AND START TIME
    write_to_log(f"Job\tStart\t" + str(datetime.datetime.now()))

    ##CREATE DUMMY FILE WITH ONE LINE
    dummy_file = open(os.path.join(output_folder, ".dummy.txt"), "w")
    dummy_file.write("dummy_file\n")
    dummy_file.close()
    ##CREATE EMPTY FILE
    empty_file = open(os.path.join(output_folder, ".empty.txt"), "w")
    empty_file.close()
    ##CREATE EMPTY DIR
    empty_dir = os.path.join(output_folder, ".empty")
    os.makedirs(empty_dir, exist_ok=True)

    ##extract list of chromosomes from reference genome
    tmp_list = os.listdir(ref_folder)  # type: ignore
    for f in tmp_list:
        if ".fa" in f and ".fai" not in f:
            chr_list.append(f.replace(".fa", ""))

    ##check if vcf_list is empty
    vcf_process = True
    for elem in vcf_list:
        if elem == "NULL":
            vcf_process = False
            break
        else:
            vcf_list_checked.append(elem.strip().replace("/", ""))
    write_to_verbose(f"vcf_process is: {vcf_process}")
    write_to_verbose(f"vcf_list is: {vcf_list_checked}")

    ##generate pam_seq and pam_position
    pam_complete = open(pam_file).readlines()[0].strip()  # type: ignore
    pam_seq = pam_complete.split(" ")[0]
    pam_position = int(pam_complete.split(" ")[1])
    if pam_position > 0:
        pam_seq = pam_seq[::-pam_position]
    else:
        pam_seq = pam_seq[:pam_position]

    return 0


def generate_index(genome_folder, process_indels=False):
    ##generate index for input genome
    genome_name = os.path.basename(genome_folder)
    write_to_log(f"Index-genome {genome_name}\tStart\t" + str(datetime.datetime.now()))
    ##check if genome is already indexed for all bulges from bMax to bMax+20
    for bulge in range(int(bMax), int(bMax) + 20):
        if os.path.isdir(
            os.path.join(
                genomes_libraries_folder, pam_seq, "_", str(bulge), "_", genome_name
            )
        ):
            write_to_verbose("genome already indexed")
            return 0

    index_run = f"'crispritz.py' 'index-genome' {genome_name} {genome_folder} {pam_file} '-bMax' {bMax} '-th' {ncpus}"
    code = subprocess.run(index_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("index-genome failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    ##if process_indels is True, index genome for indels
    if process_indels:
        index_run = [
            "crispritz.py",
            "index-genome",
            genome_name + "_INDELS",
            os.path.join(genome_folder, "_INDELS"),
            pam_file,
            "-bMax",
            bMax,
            "-th",
            ncpus,
        ]

    write_to_log(f"Index-genome {genome_name}\tEnd\t" + str(datetime.datetime.now()))
    return 0


def generate_dict(vcf_data):
    vcf_name = vcf_data
    write_to_log(
        f"Add-variants for VCF {vcf_name}\tStart\t" + str(datetime.datetime.now())
    )

    if os.path.isdir(os.path.join(genomes_folder, ref_name + "+" + vcf_name)):
        write_to_verbose("variants already added")
        write_to_log(
            f"Add-variants for VCF {vcf_name}\tEnd\t" + str(datetime.datetime.now())
        )
        return 0

    write_to_verbose(f"name of genome is: {ref_name}")

    write_to_verbose(f"Starting dictionary generation for vcf: {vcf_name}")
    variant_run = f"'crispritz.py' 'add-variants' {os.path.join(vcfs_folder,vcf_data)} {ref_folder} 'true'"
    # variant_run = ["crispritz.py","add-variants",os.path.join(vcfs_folder,vcf_data),ref_folder,"true"]
    code = subprocess.run(variant_run, shell=True, capture_output=True)
    if code.returncode != 0:
        write_to_error("add-variants failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)
    ##generate fake chr list for INDELS analysis
    tmp_list = os.listdir(os.path.join(vcfs_folder, vcf_data))
    for f in tmp_list:
        if ".vcf.gz" in f:
            split = f.split(".")
            for elem in split:
                if "chr" in elem:
                    fake_chr_list.append("fake" + elem)

    ##rename indexed variant genome folder
    shutil.move(
        os.path.join(current_working_directory, "variants_genome"),
        os.path.join(genomes_folder, "variants_genome"),
    )
    os.rename(
        os.path.join(
            genomes_folder, "variants_genome", "SNPs_genome", f"{ref_name}_enriched"
        ),
        os.path.join(genomes_folder, ref_name + "+" + vcf_name),
    )
    os.makedirs(
        os.path.join(genomes_folder, f"{ref_name}+{vcf_name}_INDELS"), exist_ok=True
    )
    os.makedirs(
        os.path.join(dictionaries_folder, f"dictionaries_{vcf_name}"), exist_ok=True
    )
    os.makedirs(
        os.path.join(dictionaries_folder, f"log_indels_{vcf_name}"), exist_ok=True
    )

    list_files = os.listdir(
        os.path.join(genomes_folder, "variants_genome", "SNPs_genome")
    )
    for file in list_files:
        if ".json" in file:
            os.rename(
                os.path.join(genomes_folder, "variants_genome", "SNPs_genome", file),
                os.path.join(dictionaries_folder, f"dictionaries_{vcf_name}", file),
            )
        if "log" in file:
            os.rename(
                os.path.join(genomes_folder, "variants_genome", "SNPs_genome", file),
                os.path.join(dictionaries_folder, f"log_indels_{vcf_name}", file),
            )

    for fakechr in fake_chr_list:
        os.rename(
            os.path.join(
                genomes_folder,
                "variants_genome",
                f"fake_{vcf_name}_{fakechr.replace('fake','')}",
                f"{fakechr}.fa",
            ),
            os.path.join(
                genomes_folder, f"{ref_name}+{vcf_name}_INDELS", f"{fakechr}.fa"
            ),
        )
    ##remove temporary folder for variant genome
    os.removedirs(os.path.join(genomes_folder, "variants_genome"))
    write_to_log(
        f"Add-variants for VCF {vcf_name}\tEnd\t" + str(datetime.datetime.now())
    )


##START PROCESS FROM SCRATCH AND CHECK IF ANY STEP CAN BE SKIPPED
code = pre_process()
if code != 0:
    write_to_error("pre_process failed")
    sys.exit(1)
generate_index(ref_folder, False)  ##generate index for reference genome
for vcf_data in vcf_list_checked:
    if len(vcf_data):
        pass
    else:
        continue
    generate_dict(vcf_data)  ##generate dictionary for vcf
    generate_index(
        os.path.join(genomes_folder, ref_name, "+", vcf_data), True
    )  ##generate index for vcf genome
