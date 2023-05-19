#!/usr/bin/env python
#
import sys
import datetime
import os
import subprocess
import shutil
import pandas as pd
import new_simple_analysis as nsa

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
pam_seq = ""
pam_name = ""
guide_name = ""
annotation_name = ""
ref_name = os.path.basename(ref_folder).replace("/", "")
output_folder_name = os.path.basename(output_folder).replace("/", "")
bestCFD_file = os.path.join(output_folder, output_folder_name + ".bestCFD.txt")
bestCRISTA_file = os.path.join(output_folder, output_folder_name + ".bestCRISTA.txt")
bestMMBUL_file = os.path.join(output_folder, output_folder_name + ".bestmmblg.txt")


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

    ##prepare output files
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

    bestCFD = open(
        os.path.join(output_folder, output_folder_name + ".bestCFD.txt"), "w"
    )
    bestCRISTA = open(
        os.path.join(output_folder, output_folder_name + ".bestCRISTA.txt"), "w"
    )
    bestMMBUL = open(
        os.path.join(output_folder, output_folder_name + ".bestmmblg.txt"), "w"
    )
    bestCFD.write("\t".join(header) + "\n")
    bestCRISTA.write("\t".join(header) + "\n")
    bestMMBUL.write("\t".join(header) + "\n")
    bestCFD.close()
    bestCRISTA.close()
    bestMMBUL.close()

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
    write_to_log(
        f"Add-variants for VCF {vcf_name}\tStart\t" + str(datetime.datetime.now())
    )

    if os.path.isdir(os.path.join(dictionaries_folder, f"dictionaries_{vcf_name}")):
        write_to_verbose("variants already added")
        write_to_log(
            f"Add-variants for VCF {vcf_name}\tEnd\t" + str(datetime.datetime.now())
        )
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
    global fake_chr_list
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


def post_process(target_file: str, vcf_data: str, ref_only: bool = False) -> None:
    """function to post process the results of crispritz

    Args:
        target_file (str): file generated by crispritz containing all targets
        vcf_data (str): vcf name input by user to process variants
        ref_only (bool, optional): check if analysis is ref only (no variant). Defaults to False.

    Returns:
        None: [description]
    """

    write_to_verbose(f"Starting post process")
    write_to_verbose(f"target_file is: {target_file}")

    write_to_log(f"Post Process\tStart\t" + str(datetime.datetime.now()))
    lists_of_targets_list = [[], [], []]

    target_df = pd.read_csv(os.path.join(output_folder, target_file), sep="\t")
    for chr in chr_list:
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
            allowed_mms=mm,
        )
        ##return list of lists with targets scored by CFD,MMBUL,CRISTA
        lists_of_targets_list = nsa.start_processing(target_df_chr, data_to_process)

        # write to bestFILES to check
        file = open(bestCFD_file, "a")
        file.write("".join(lists_of_targets_list[0]))
        file.close()

        file = open(bestMMBUL_file, "a")
        file.write("".join(lists_of_targets_list[1]))
        file.close()

        file = open(bestCRISTA_file, "a")
        file.write("".join(lists_of_targets_list[2]))
        file.close()

    # print(lists_of_targets_list[0])
    # snp_analysis_run = f"./new_simple_analysis.py {os.path.join(ref_folder,chr+'.fa')} {os.path.join(dictionaries_folder,'dictionaries_'+vcf_data,'my_dict_' + chr + '.json')} {os.path.join(output_folder,chr+'_process_before_simple_analysis.txt')} {pam_file} {os.path.join(output_folder,output_folder_name)} {mm}"
    # code = subprocess.run(snp_analysis_run, shell=True, capture_output=True)

    # write_to_verbose(code.stdout.decode("utf-8"))
    # if code.returncode != 0:
    #     write_to_error("simple analysis failed")
    #     write_to_error(code.stderr.decode("utf-8"))
    #     sys.exit(1)

    # os.chdir(current_working_directory)
    write_to_log(f"Post Process\tEnd\t" + str(datetime.datetime.now()))
    write_to_verbose(f"Post Process END")


def post_process_indels(target_file, vcf_data, ref_only=False):
    write_to_verbose(f"Starting post process indels")
    write_to_verbose(f"target_file is: {target_file}")

    write_to_log(f"Post Process Indels\tStart\t" + str(datetime.datetime.now()))

    target_df = pd.read_csv(os.path.join(output_folder, target_file), sep="\t")
    # write_to_verbose(f"fake_chr_list is: {fake_chr_list}")

    for chr in fake_chr_list:
        # fake_chr = "fake" + chr
        target_df_chr = target_df.loc[target_df["Chromosome"] == chr]
        target_df_chr["PAM_gen"] = "n"
        target_df_chr["Var_uniq"] = "n"
        target_df_chr["Samples"] = "n"
        target_df_chr["Annotation_type"] = "n"
        target_df_chr["Real_Guide"] = target_df_chr["crRNA"].str.replace("-", "")
        target_df_chr["rsID"] = "n"
        target_df_chr["AF"] = "n"
        target_df_chr["SNP_position"] = "n"
        target_df_chr.to_csv(
            os.path.join(output_folder, chr + "_process_before_simple_analysis.txt"),
            sep="\t",
            index=False,
        )

        os.chdir(processes_dir)
        # indel_analysis_run = "./analisi_indels_NNN.py"
        indel_analysis_run = f"./analisi_indels_NNN.py {os.path.join(output_folder,'.empty.txt')} {os.path.join(output_folder,chr+'_process_before_simple_analysis.txt')} {os.path.join(output_folder,output_folder_name)} {os.path.join(dictionaries_folder,'log_indels_'+vcf_data)} {pam_file} {mm} {os.path.join(ref_folder,chr.replace('fake','')+'.fa')} {guide_file} {bDNA} {bRNA}"
        code = subprocess.run(indel_analysis_run, shell=True, capture_output=True)

        write_to_verbose(code.stdout.decode("utf-8"))
        if code.returncode != 0:
            write_to_error("indel analysis failed")
            write_to_error(code.stderr.decode("utf-8"))
            sys.exit(1)

    write_to_verbose(f"Post process indel ended correctly")
    os.chdir(current_working_directory)

    write_to_log(f"Post Process Indels\tEnd\t" + str(datetime.datetime.now()))

    return 0


def fix_columns(output_folder_name):
    write_to_verbose(f"Starting fix columns in best files")

    os.chdir(processes_dir)
    adjust_col_run = f"./adjust_cols.py {os.path.join(output_folder,output_folder_name+'.bestCFD.txt')}"
    code = subprocess.run(adjust_col_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("adjust bestCFD failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    adjust_col_run = f"./adjust_cols.py {os.path.join(output_folder,output_folder_name+'.bestCRISTA.txt')}"
    code = subprocess.run(adjust_col_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("adjust bestCRISTA failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    adjust_col_run = f"./adjust_cols.py {os.path.join(output_folder,output_folder_name+'.bestmmblg.txt')}"
    code = subprocess.run(adjust_col_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("adjust bestMMBUL failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    os.chdir(current_working_directory)
    return 0


def remove_bad_indels(output_folder_name):
    write_to_verbose(f"Starting removing bad indels in best files")

    os.chdir(processes_dir)
    remove_bad_indels_run = f"./remove_bad_indel_targets.py {os.path.join(output_folder,output_folder_name+'.bestCFD.txt')}"
    code = subprocess.run(remove_bad_indels_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("bad indels removing bestCFD failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    remove_bad_indels_run = f"./remove_bad_indel_targets.py {os.path.join(output_folder,output_folder_name+'.bestCRISTA.txt')}"
    code = subprocess.run(remove_bad_indels_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("bad indels removing bestCRISTA failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    remove_bad_indels_run = f"./remove_bad_indel_targets.py {os.path.join(output_folder,output_folder_name+'.bestmmblg.txt')}"
    code = subprocess.run(remove_bad_indels_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("bad indels removing bestMMBUL failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    os.chdir(current_working_directory)
    write_to_verbose(f"Removing bad indels in best files ended correctly")
    return 0


def merge_results(output_folder_name):
    ##positional arguments to start merge results
    chrom = 5  # column for chromosome
    position = 7  # column for cluster_position
    total = 11  # column for total (mm+bul)
    true_guide = 16  # column for true guide (original guide without bulges)
    snp_info = 19  # column for snp info
    cfd = 21  # column for cfd score

    # sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
    # tail -n +2 $final_res.bestCRISTA.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCRISTA.txt
    # #sort using guide_seq,chr,cluster_pos,total(mm+bul)
    # head -1 $final_res.bestmmblg.txt >$final_res.tmp
    # tail -n +2 $final_res.bestmmblg.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestmmblg.txt

    write_to_verbose(f"Start merging results in best files")
    os.chdir(processes_dir)
    to_merge_df = pd.read_csv(
        os.path.join(output_folder, output_folder_name + ".bestCFD.txt"), sep="\t"
    )
    to_merge_df.sort_values(
        ["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    to_merge_df.to_csv(
        os.path.join(output_folder, output_folder_name + ".bestCFD_sorted.txt"),
        sep="\t",
        index=False,
    )
    merge_resuts_run = f"./remove_contiguous_samples_cfd.py {os.path.join(output_folder,output_folder_name+'.bestCFD_sorted.txt')} {os.path.join(output_folder,output_folder_name+'.bestCFD.trimmed.txt')} {merge_t} {chrom} {position} {total} {true_guide} {snp_info} {cfd} 'score'"
    code = subprocess.run(merge_resuts_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("merge bestCFD failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    to_merge_df = pd.read_csv(
        os.path.join(output_folder, output_folder_name + ".bestCRISTA.txt"), sep="\t"
    )
    to_merge_df.sort_values(
        ["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    to_merge_df.to_csv(
        os.path.join(output_folder, output_folder_name + ".bestCRISTA_sorted.txt"),
        sep="\t",
        index=False,
    )
    merge_resuts_run = f"./remove_contiguous_samples_cfd.py {os.path.join(output_folder,output_folder_name+'.bestCRISTA_sorted.txt')} {os.path.join(output_folder,output_folder_name+'.bestCRISTA.trimmed.txt')} {merge_t} {chrom} {position} {total} {true_guide} {snp_info} {cfd} 'score'"
    code = subprocess.run(merge_resuts_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("merge bestCFD failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    to_merge_df = pd.read_csv(
        os.path.join(output_folder, output_folder_name + ".bestmmblg.txt"), sep="\t"
    )
    to_merge_df.sort_values(
        ["Real_Guide", "Chromosome", "Cluster_Position"],
        ascending=[True, True, True],
        inplace=True,
    )
    to_merge_df.to_csv(
        os.path.join(output_folder, output_folder_name + ".bestmmblg_sorted.txt"),
        sep="\t",
        index=False,
    )
    merge_resuts_run = f"./remove_contiguous_samples_cfd.py {os.path.join(output_folder,output_folder_name+'.bestmmblg_sorted.txt')} {os.path.join(output_folder,output_folder_name+'.bestmmblg.trimmed.txt')} {merge_t} {chrom} {position} {total} {true_guide} {snp_info} {cfd} 'total'"
    code = subprocess.run(merge_resuts_run, shell=True, capture_output=True)
    write_to_verbose(code.stdout.decode("utf-8"))
    if code.returncode != 0:
        write_to_error("merge bestCFD failed")
        write_to_error(code.stderr.decode("utf-8"))
        sys.exit(1)

    write_to_verbose(f"merging results in best files completed")
    os.chdir(current_working_directory)
    return 0


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
post_process(target_file, "", ref_only=True)

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
        target_file.replace(ref_name, ref_name + "_" + vcf_data), vcf_data, False
    )
    # post_process_indels(
    #     target_file.replace(ref_name, ref_name + "_" + vcf_data + "_INDELS"),
    #     vcf_data,
    #     False,
    # )


##fix columns in best files
# fix_columns(output_folder_name)
# remove_bad_indels(output_folder_name)
# merge_results(output_folder_name)
