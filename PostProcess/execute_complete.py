#!/usr/bin/env python
#
import sys
import datetime
import os
import subprocess
import shutil
import pandas as pd
import new_simple_analysis as nsa
import adjust_cols as ac
import analisi_indels_NNN as ain
import remove_bad_indel_targets as rindel
import remove_contiguous_samples_cfd as merge
import threading


# set -e # capture any failure

# file for automated search of guide+pam in reference and variant genomes
ref_folder = sys.argv[1]
vcf_list = sys.argv[2]
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
log_verbose = os.path.join(output_folder, "log_verbose.txt")
log_error = os.path.join(output_folder, "log_error.txt")
pam_seq = ""
pam_name = ""
guide_name = ""
annotation_name = ""
ref_name = os.path.basename(ref_folder).replace("/", "")
output_folder_name = os.path.basename(output_folder).replace("/", "")
bestCFD_file = os.path.join(output_folder, output_folder_name + ".bestCFD.txt")
bestCRISTA_file = os.path.join(output_folder, output_folder_name + ".bestCRISTA.txt")
bestMMBUL_file = os.path.join(output_folder, output_folder_name + ".bestMMBLG.txt")
altCFD_file = os.path.join(output_folder, output_folder_name + ".altCFD.txt")
altCRISTA_file = os.path.join(output_folder, output_folder_name + ".altCRISTA.txt")
altMMBUL_file = os.path.join(output_folder, output_folder_name + ".altMMBLG.txt")
chr_df_dict = dict()
bestCFD_df = pd.DataFrame()
bestCRISTA_df = pd.DataFrame()
bestMMBUL_df = pd.DataFrame()
altCFD_df = pd.DataFrame()
altCRISTA_df = pd.DataFrame()
altMMBUL_df = pd.DataFrame()
header = [
    "#Bulge_type",
    "crRNA",
    "DNA",
    "Chromosome",
    "Position",
    "Cluster_Position",
    "Direction",
    "Mismatches",
    "Bulge_Size",
    "Total",
    "PAM_gen",
    "Var_uniq",
    "Samples",
    "Annotation_Type",
    "Real_Guide",
    "rsID",
    "AF",
    "SNP",
    "Reference",
    "CFD_ref",
    "CFD",
    "#Seq_in_cluster",
]


##USER FUNCTIONS
def write_to_log(message):
    ##write to log file with termination
    log_file.write(message + "\n")


def write_to_verbose(message):
    ##write to log verbose file with print autotermination
    print(message, file=sys.stdout)


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
    global chr_list
    for f in tmp_list:
        if ".fa" in f and ".fai" not in f:
            chr_list.append(f.replace(".fa", ""))

    ## process vcf list to extract only non empty lines
    global vcf_process
    vcf_process = False
    for elem in open(vcf_list).readlines():
        if len(elem.strip().replace("/", "")) > 0:
            vcf_list_checked.append(elem.strip().replace("/", ""))
            vcf_process = True

    write_to_verbose(f"vcf_process is: {vcf_process}")
    write_to_verbose(f"vcf_list is: {vcf_list_checked}")

    ##generate pam_seq and pam_position
    pam_complete = open(pam_file).readlines()[0].strip()  # type: ignore
    global pam_seq
    pam_seq = pam_complete.split(" ")[0]
    pam_position = int(pam_complete.split(" ")[1])
    if pam_position > 0:
        pam_seq = pam_seq[-pam_position:]
    else:
        pam_seq = pam_seq[:pam_position]

    global pam_name
    global guide_name
    pam_name = os.path.basename(pam_file).replace("/", "")
    guide_name = os.path.basename(guide_file).replace("/", "")

    return 0


def generate_index(genome_folder, process_indels=False):
    ##generate index for input genome
    genome_name = os.path.basename(genome_folder)
    write_to_log(f"Index-genome {genome_name}\tStart\t" + str(datetime.datetime.now()))
    ##check if genome is already indexed for all bulges from bMax to bMax+20
    for bulge in range(int(bMax), int(bMax) + 20):
        write_to_verbose(
            f"checking presence of indexed genome for {pam_seq}_{str(bulge)}_{genome_name}"
        )
        if os.path.isdir(
            os.path.join(
                genomes_libraries_folder, f"{pam_seq}_{str(bulge)}_{genome_name}"
            )
        ):
            write_to_verbose("genome already indexed")
            write_to_log(
                f"Index-genome {genome_name}\tEnd\t" + str(datetime.datetime.now())
            )
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
        index_run = f"'crispritz.py' 'index-genome' {genome_name}_INDELS {genome_folder}_INDELS {pam_file} '-bMax' {bMax} '-th' {ncpus}"
        code = subprocess.run(index_run, shell=True, capture_output=True)
        write_to_verbose(code.stdout.decode("utf-8"))
        if code.returncode != 0:
            write_to_error("index-genome for indels failed")
            write_to_error(code.stderr.decode("utf-8"))
            sys.exit(1)

    write_to_log(f"Index-genome {genome_name}\tEnd\t" + str(datetime.datetime.now()))
    return 0


def generate_dict(vcf_data):
    vcf_name = vcf_data
    global fake_chr_list

    write_to_log(
        f"Add-variants for VCF {vcf_name}\tStart\t" + str(datetime.datetime.now())
    )

    if os.path.isdir(os.path.join(dictionaries_folder, f"dictionaries_{vcf_name}")):
        write_to_verbose("variants already added")
        write_to_log(
            f"Add-variants for VCF {vcf_name}\tEnd\t" + str(datetime.datetime.now())
        )
        fake_chr_list = os.listdir(
            os.path.join(dictionaries_folder, f"dictionaries_{vcf_name}")
        )
        fake_chr_list = [x.replace("my_dict_", "") for x in fake_chr_list]
        fake_chr_list = [x.replace(".json", "") for x in fake_chr_list]

        return 0

    write_to_verbose(f"name of genome is: {ref_name}")

    os.chdir(output_folder)
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

    os.chdir(current_working_directory)
    ##rename indexed variant genome folder
    shutil.move(
        os.path.join(output_folder, "variants_genome"),
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
    shutil.rmtree(os.path.join(genomes_folder, "variants_genome"), ignore_errors=True)
    write_to_log(
        f"Add-variants for VCF {vcf_name}\tEnd\t" + str(datetime.datetime.now())
    )


def search(ref_name, vcf_data, pam_seq, bMax, ncpus, mm, pam_name, do_ref=False):
    ##move to output folder to save results
    os.chdir(output_folder)
    ##extract name for index ref
    idx_ref = ""
    for bulge in range(int(bMax), int(bMax) + 20):
        if os.path.isdir(
            os.path.join(genomes_libraries_folder, f"{pam_seq}_{str(bulge)}_{ref_name}")
        ):
            idx_ref = os.path.join(
                genomes_libraries_folder, f"{pam_seq}_{str(bulge)}_{ref_name}"
            )
            break
    if do_ref:  ##if ref not already done, search it
        write_to_verbose(
            f"Starting search for genome: {ref_name}, vcf process is: {vcf_process}"
        )
        write_to_log(f"Search Reference\tStart\t" + str(datetime.datetime.now()))
        write_to_verbose(f"idx_ref is: {idx_ref}")
        ref_search_run = f"'crispritz.py' 'search' {idx_ref} {pam_file} {guide_file} {ref_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} '-mm' {mm} '-bDNA' {bDNA} '-bRNA' {bRNA} '-th' {ncpus} '-t'"
        code = subprocess.run(ref_search_run, shell=True, capture_output=True)
        write_to_verbose(code.stdout.decode("utf-8"))
        if code.returncode != 0:
            write_to_error("reference search failed")
            write_to_error(code.stderr.decode("utf-8"))
            sys.exit(1)
        write_to_log(f"Search Reference\tEnd\t" + str(datetime.datetime.now()))

    if vcf_process and len(vcf_data):  ##vcf is process and not empty, search variant
        write_to_verbose(
            f"Starting search for genome: {ref_name}+{vcf_data}, vcf process is: {vcf_process}"
        )
        write_to_log(f"Search Variant\tStart\t" + str(datetime.datetime.now()))
        idx_var = idx_ref.replace(ref_name, ref_name + "+" + vcf_data)
        write_to_verbose(f"idx_var is: {idx_var}")
        var_search_run = f"'crispritz.py' 'search' {idx_var} {pam_file} {guide_file} {ref_name}_{vcf_data}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} '-mm' {mm} '-bDNA' {bDNA} '-bRNA' {bRNA} '-th' {ncpus} '-t'"
        code = subprocess.run(var_search_run, shell=True, capture_output=True)
        write_to_verbose(code.stdout.decode("utf-8"))
        if code.returncode != 0:
            write_to_error("variant search failed")
            write_to_error(code.stderr.decode("utf-8"))
            sys.exit(1)
        write_to_log(f"Search Variant\tEnd\t" + str(datetime.datetime.now()))

        write_to_log(f"Search Indels\tStart\t" + str(datetime.datetime.now()))
        idx_indels = idx_var.replace(
            ref_name + "+" + vcf_data, ref_name + "+" + vcf_data + "_INDELS"
        )
        indel_search_run = f"'crispritz.py' 'search' {idx_indels} {pam_file} {guide_file} {ref_name}_{vcf_data}_INDELS_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA} '-mm' {mm} '-bDNA' {bDNA} '-bRNA' {bRNA} '-th' {ncpus} '-t'"
        code = subprocess.run(indel_search_run, shell=True, capture_output=True)
        write_to_verbose(code.stdout.decode("utf-8"))
        if code.returncode != 0:
            write_to_error("indels search failed")
            write_to_error(code.stderr.decode("utf-8"))
            sys.exit(1)

        write_to_log(f"Search Indels\tEnds\t" + str(datetime.datetime.now()))

    os.chdir(current_working_directory)
    return 0


def variant_analisys(target_df: pd.DataFrame, chr: str, vcf_data: str):
    global chr_df_dict
    target_df_chr = target_df.loc[target_df["Chromosome"] == chr]
    target_df_chr["PAM_gen"] = "n"
    target_df_chr["Var_uniq"] = "n"
    target_df_chr["Samples"] = "n"
    target_df_chr["Annotation_type"] = "n"
    target_df_chr["Real_Guide"] = target_df_chr["crRNA"].str.replace("-", "")
    target_df_chr["rsID"] = "n"
    target_df_chr["AF"] = "n"
    target_df_chr["SNP_position"] = "n"

    ## convert df to list to be processed
    target_df_chr = target_df_chr.values.tolist()
    data_to_process = nsa.init(
        fasta_file=os.path.join(ref_folder, chr + ".fa"),
        pam_file=pam_file,
        dictionary_file=os.path.join(
            dictionaries_folder,
            "dictionaries_" + vcf_data,
            "my_dict_" + chr + ".json",
        ),
        allowed_mms=int(mm),
    )
    ##return list of lists with targets scored by CFD,MMBUL,CRISTA
    lists_of_targets_list = nsa.start_processing(target_df_chr, data_to_process)

    ##convert list of lists to df
    df_CFD = pd.DataFrame(lists_of_targets_list[0], columns=header)
    df_MMBUL = pd.DataFrame(lists_of_targets_list[1], columns=header)
    df_CRISTA = pd.DataFrame(lists_of_targets_list[2], columns=header)  # type: ignore
    chr_df_dict[chr + "_CFD"] = df_CFD
    chr_df_dict[chr + "_MMBUL"] = df_MMBUL
    chr_df_dict[chr + "_CRISTA"] = df_CRISTA


def post_process(
    target_file: str,
    vcf_data: str,
    ref_only: bool = False,
) -> None:
    """_summary_

    Args:
        target_file (str): _description_
        vcf_data (str): _description_
        bestCFD_df (pd.DataFrame): _description_
        bestCRISTA_df (pd.DataFrame): _description_
        bestMMBUL_df (pd.DataFrame): _description_
        ref_only (bool, optional): _description_. Defaults to False.
    """

    write_to_verbose(f"Starting post process")
    write_to_verbose(f"target_file is: {target_file}")
    write_to_log(f"Post Process\tStart\t" + str(datetime.datetime.now()))

    ##scope variables
    target_df = pd.read_csv(os.path.join(output_folder, target_file), sep="\t")
    main_thread = threading.main_thread()
    ##global variables
    global bestCFD_df
    global bestCRISTA_df
    global bestMMBUL_df
    global chr_df_dict

    for chr in chr_list:
        t = threading.Thread(
            target=variant_analisys,
            args=(target_df, chr, vcf_data),
        )
        t.start()
        # target_df_chr = target_df.loc[target_df["Chromosome"] == chr]
        # target_df_chr["PAM_gen"] = "n"
        # target_df_chr["Var_uniq"] = "n"
        # target_df_chr["Samples"] = "n"
        # target_df_chr["Annotation_type"] = "n"
        # target_df_chr["Real_Guide"] = target_df_chr["crRNA"].str.replace("-", "")
        # target_df_chr["rsID"] = "n"
        # target_df_chr["AF"] = "n"
        # target_df_chr["SNP_position"] = "n"

        # ## convert df to list to be processed
        # target_df_chr = target_df_chr.values.tolist()
        # data_to_process = nsa.init(
        #     fasta_file=os.path.join(ref_folder, chr + ".fa"),
        #     pam_file=pam_file,
        #     dictionary_file=os.path.join(
        #         dictionaries_folder,
        #         "dictionaries_" + vcf_data,
        #         "my_dict_" + chr + ".json",
        #     ),
        #     allowed_mms=int(mm),
        # )
        # ##return list of lists with targets scored by CFD,MMBUL,CRISTA
        # lists_of_targets_list = nsa.start_processing(target_df_chr, data_to_process)

        # ##convert list of lists to df
        # df_CFD = pd.DataFrame(lists_of_targets_list[0], columns=header)
        # df_MMBUL = pd.DataFrame(lists_of_targets_list[1], columns=header)
        # df_CRISTA = pd.DataFrame(lists_of_targets_list[2], columns=header)  # type: ignore
        # chr_df_dict[chr + "_CFD"] = df_CFD
        # chr_df_dict[chr + "_MMBUL"] = df_MMBUL
        # chr_df_dict[chr + "_CRISTA"] = df_CRISTA

    for t in threading.enumerate():
        if t is main_thread:
            continue
        # logging.debug('joining %s', t.getName())
        t.join()

    to_concat = [chr_df_dict[key + "_CFD"] for key in chr_list]
    to_concat.append(bestCFD_df)
    bestCFD_df = pd.concat(to_concat, axis=0)

    to_concat = [chr_df_dict[key + "_CFD"] for key in chr_list]
    to_concat.append(bestCRISTA_df)
    bestCRISTA_df = pd.concat(to_concat, axis=0)

    to_concat = [chr_df_dict[key + "_CFD"] for key in chr_list]
    to_concat.append(bestMMBUL_df)
    bestMMBUL_df = pd.concat(to_concat, axis=0)

    write_to_log(f"Post Process\tEnd\t" + str(datetime.datetime.now()))
    write_to_verbose(f"Post Process END")


def post_process_indels(
    target_file: str,
    vcf_data: str,
    ref_only: bool = False,
) -> None:
    write_to_verbose(f"Starting post process indels")
    write_to_verbose(f"target_file is: {target_file}")
    write_to_log(f"Post Process Indels\tStart\t" + str(datetime.datetime.now()))

    global bestCFD_df
    global bestCRISTA_df
    global bestMMBUL_df

    ##tmp varaibles
    target_df = pd.read_csv(os.path.join(output_folder, target_file), sep="\t")
    chr_df_dict = dict()
    bestCFD_df_indel = pd.DataFrame()
    bestCRISTA_df_indel = pd.DataFrame()
    bestMMBUL_df_indel = pd.DataFrame()

    for chr in fake_chr_list:
        target_df_chr = target_df.loc[target_df["Chromosome"] == chr]
        target_df_chr["PAM_gen"] = "n"
        target_df_chr["Var_uniq"] = "n"
        target_df_chr["Samples"] = "n"
        target_df_chr["Annotation_type"] = "n"
        target_df_chr["Real_Guide"] = target_df_chr["crRNA"].str.replace("-", "")
        target_df_chr["rsID"] = "n"
        target_df_chr["AF"] = "n"
        target_df_chr["SNP_position"] = "n"

        ## convert df to list to be processed
        target_df_chr = target_df_chr.values.tolist()
        data_to_process = ain.init(
            fasta_path=os.path.join(ref_folder, chr + ".fa"),
            indel_dict_path=os.path.join(
                dictionaries_folder,
                "log_indels_" + vcf_data,
                "log" + chr.replace("fake", "") + ".txt",
            ),
            max_mm=int(mm),
            dna_bulges=int(bDNA),
            rna_bulges=int(bRNA),
            pam_file=pam_file,
        )
        ##return list of lists with targets scored by CFD,MMBUL,CRISTA
        lists_of_targets_list = ain.start_processing(target_df_chr, data_to_process)
        ##convert list of lists to df
        df_CFD = pd.DataFrame(lists_of_targets_list[0], columns=header)
        df_MMBUL = pd.DataFrame(lists_of_targets_list[1], columns=header)
        df_CRISTA = pd.DataFrame(lists_of_targets_list[2], columns=header)
        ##add to dict
        chr_df_dict[chr + "_CFD"] = df_CFD
        chr_df_dict[chr + "_MMBUL"] = df_MMBUL
        chr_df_dict[chr + "_CRISTA"] = df_CRISTA

    ##concat all chr dfs to one indels df
    bestCFD_df_indel = pd.concat(
        [chr_df_dict[key + "_CFD"] for key in fake_chr_list], axis=0
    )
    bestCRISTA_df_indel = pd.concat(
        [chr_df_dict[key + "_CRISTA"] for key in fake_chr_list], axis=0
    )
    bestMMBUL_df_indel = pd.concat(
        [chr_df_dict[key + "_MMBUL"] for key in fake_chr_list], axis=0
    )

    # append indels to best df
    bestCFD_df_indel = pd.DataFrame(
        rindel.remove_bad_indels(bestCFD_df_indel.values.tolist()),
        columns=bestCFD_df_indel.columns.tolist(),
    )
    bestCFD_df = pd.concat([bestCFD_df, bestCFD_df_indel], axis=0)

    bestCRISTA_df_indel = pd.DataFrame(
        rindel.remove_bad_indels(bestCRISTA_df_indel.values.tolist()),
        columns=bestCRISTA_df_indel.columns.tolist(),
    )
    bestCRISTA_df = pd.concat([bestCRISTA_df, bestCRISTA_df_indel], axis=0)

    bestMMBUL_df_indel = pd.DataFrame(
        rindel.remove_bad_indels(bestMMBUL_df_indel.values.tolist()),
        columns=bestMMBUL_df_indel.columns.tolist(),
    )
    bestMMBUL_df = pd.concat([bestMMBUL_df, bestMMBUL_df_indel], axis=0)

    write_to_verbose(f"Post process INDELs END")
    write_to_log(f"Post Process Indels\tEnd\t" + str(datetime.datetime.now()))


def fix_columns():
    write_to_log(f"Fix Columns\tStart\t" + str(datetime.datetime.now()))
    write_to_verbose(f"Starting fix columns in best files")

    global bestCFD_df
    global bestCRISTA_df
    global bestMMBUL_df

    # adjust cols to final df
    bestCFD_df = ac.order_cols(bestCFD_df)
    bestCRISTA_df = ac.order_cols(bestCRISTA_df)
    bestMMBUL_df = ac.order_cols(bestMMBUL_df)

    write_to_log(f"Fix Columns\tEnd\t" + str(datetime.datetime.now()))
    write_to_verbose(f"Fix columns END")


def merge_results():
    write_to_verbose(f"Starting merge results")
    write_to_log(f"Merge Results\tStart\t" + str(datetime.datetime.now()))

    global bestCFD_df
    global bestCRISTA_df
    global bestMMBUL_df
    global altCFD_df
    global altCRISTA_df
    global altMMBUL_df

    ##tmp variables
    best_list = list()
    discard_list = list()
    header_dict = dict()
    for count, key in enumerate(bestCFD_df.columns.tolist()):
        header_dict[key] = count

    ##CFD PROCESSING
    bestCFD_df.sort_values(
        by=["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    best_list, discard_list = merge.merge_results(
        bestCFD_df.values.tolist(),
        tau=int(merge_t),
        sort_order="score",
        header=header_dict,
    )
    bestCFD_df = pd.DataFrame(best_list, columns=bestCFD_df.columns.tolist())
    altCFD_df = pd.DataFrame(discard_list, columns=bestCFD_df.columns.tolist())
    bestCFD_df.to_csv(bestCFD_file, sep="\t", index=False, mode="w")
    altCFD_df.to_csv(altCFD_file, sep="\t", index=False, mode="w")

    ##CRISTA PROCESSING
    bestCRISTA_df.sort_values(
        by=["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    best_list, discard_list = merge.merge_results(
        bestCRISTA_df.values.tolist(),
        tau=int(merge_t),
        sort_order="score",
        header=header_dict,
    )
    bestCRISTA_df = pd.DataFrame(best_list, columns=bestCRISTA_df.columns.tolist())
    altCRISTA_df = pd.DataFrame(discard_list, columns=bestCRISTA_df.columns.tolist())
    bestCRISTA_df.to_csv(bestCRISTA_file, sep="\t", index=False, mode="w")
    altCRISTA_df.to_csv(altCRISTA_file, sep="\t", index=False, mode="w")

    ##MMBUL PROCESSING
    bestMMBUL_df.sort_values(
        by=["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    best_list, discard_list = merge.merge_results(
        bestMMBUL_df.values.tolist(),
        tau=int(merge_t),
        sort_order="score",
        header=header_dict,
    )
    bestMMBUL_df = pd.DataFrame(best_list, columns=bestMMBUL_df.columns.tolist())
    altMMBUL_df = pd.DataFrame(discard_list, columns=bestMMBUL_df.columns.tolist())
    bestMMBUL_df.to_csv(bestMMBUL_file, sep="\t", index=False, mode="w")
    altMMBUL_df.to_csv(altMMBUL_file, sep="\t", index=False, mode="w")

    write_to_verbose(f"Merge results END")
    write_to_log(f"Merge Results\tEnd\t" + str(datetime.datetime.now()))


##START PROCESS FROM SCRATCH AND CHECK IF ANY STEP CAN BE SKIPPED
code = pre_process()
if code != 0:
    write_to_error("pre_process failed")
    sys.exit(1)

generate_index(ref_folder, process_indels=False)  ##generate index for reference genome
search(
    ref_name, "", pam_seq, bMax, ncpus, mm, pam_name, do_ref=True
)  ##search on reference genome

target_file = f"{ref_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
post_process(
    target_file,
    vcf_data="",
    ref_only=True,
)  ##post process for reference genome

##start process for vcf data if any
for vcf_data in vcf_list_checked:
    if len(vcf_data):
        pass
    else:
        continue
    generate_dict(vcf_data)  ##generate dictionary for vcf
    generate_index(
        os.path.join(genomes_folder, f"{ref_name}+{vcf_data}"), True
    )  ##generate index for vcf genome
    search(
        ref_name, vcf_data, pam_seq, bMax, ncpus, mm, pam_name, False
    )  ##search on vcf genome
    post_process(
        target_file.replace(ref_name, ref_name + "_" + vcf_data),
        vcf_data,
        ref_only=False,
    )
    post_process_indels(
        target_file.replace(ref_name, ref_name + "_" + vcf_data + "_INDELS"),
        vcf_data,
        ref_only=False,
    )

fix_columns()
merge_results()
