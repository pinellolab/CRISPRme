#!/usr/bin/env python

from operator import not_, truediv
import operator
from posixpath import expanduser
from intervaltree import IntervalTree
import sys
import time
import glob
import subprocess
import os


def rev_comp(a):
    # perform reverse complement on a string
    if a == "A" or a == "a":
        return "T"
    if a == "T" or a == "t":
        return "A"
    if a == "C" or a == "c":
        return "G"
    return "C"


# count A and C presence in user window
def base_counting(target_sequence: str, start: int, end: int) -> list:
    countA = 0
    countC = 0
    countG = 0
    countT = 0

    for elem in target_sequence[start:end]:
        # skip if gap
        if elem == "-":
            continue
        # check if elem is any nucleotide
        if elem.upper() == "A":
            countA += 1
            continue
        if elem.upper() == "C":
            countC += 1
            continue
        if elem.upper() == "G":
            countG += 1
            continue
        if elem.upper() == "T":
            countT += 1
            continue

    return [countA, countC, countG, countT]


# process seed analisys
def seed_processing(seed_ref: str, seed_alt: str, non_seed_ref: str, non_seed_alt: str):
    # count mm and bulges for seed and non-seed target
    count_seed_ref = 0
    count_seed_alt = 0
    count_non_seed_ref = 0
    count_non_seed_alt = 0

    for elem in seed_ref:
        if elem.islower():
            count_seed_ref += 1
            continue
        if elem == "-":
            count_seed_ref += 1
    for elem in non_seed_ref:
        if elem.islower():
            count_non_seed_ref += 1
            continue
        if elem == "-":
            count_non_seed_ref += 1
    for elem in seed_alt:
        if elem.islower():
            count_seed_alt += 1
            continue
        if elem == "-":
            count_seed_alt += 1
    for elem in non_seed_alt:
        if elem.islower():
            count_non_seed_alt += 1
            continue
        if elem == "-":
            count_non_seed_alt += 1
    return [count_seed_ref, count_seed_alt, count_non_seed_ref, count_non_seed_alt]


def createBedforMultiAlternative(variantList, samples):
    # 5096 conta degli alleli di 1000G
    # bcftools view -H -r 1:222531770 -s HG00001 originalVCFs/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
    # list containing snp position to search in the vcf
    toSearchList = list()
    # bed file with vcf data extracted
    vcfBedFile = outputDir + originFileName + ".temp.bed"
    # cycle over variant in list to extract the correct ones from vcf
    for elem in variantList:
        split = elem.strip().split("_")
        chrom = split[0]
        pos = split[1]
        vcfFile = ""
        for vcf in vcfList:
            if chrom + "." in vcf:
                vcfFile = vcf
        vcfChr = chrom.replace("chr", "")
        toSearchList.append(vcfChr + ":" + pos)
    # string with all the snps to search in the vcf file
    toSearchString = ",".join(toSearchList).strip()
    # subprocces run to call bcftools and extract the desidered snps
    subprocess.run(
        [
            "bcftools",
            "view",
            "-H",
            "-r",
            toSearchString,
            "-s",
            samples,
            vcfFile,
            "-o",
            str(vcfBedFile),
        ]
    )
    haplotypeDict = dict()
    allele1 = list()
    allele2 = list()
    haplotype = 0
    # read line from vcf exctraction to obtain the alleles
    for line in open(vcfBedFile, "r"):
        if "INDEL" in line:
            continue
        split = line.strip().split("\t")
        tempHaploList = list()
        for elem in split:
            if "|" in elem:
                tempHaploList.append(elem)
        haplotypeDict[split[1]] = tempHaploList
        allele1 = [1] * len(tempHaploList)
        allele2 = [1] * len(tempHaploList)
    # count the correct haplotype frequence by and with sum operation
    # e.g. 0|1 and 1|1 will return 0|1 that will be counted as 1 for that sample haplotype
    for snp in haplotypeDict:
        for count, sample in enumerate(haplotypeDict[snp]):
            allele1_snp = sample.strip().split("|")[0]
            allele2_snp = sample.strip().split("|")[1]
            allele1[count] = int(allele1[count]) and int(allele1_snp)
            allele2[count] = int(allele2[count]) and int(allele2_snp)
    # final sum to generete haplotype frequency of the specific target
    for count, allele in enumerate(allele1):
        haplotype += int(allele1[count]) + int(allele2[count])
    # return to save the haplo freq divided by the number of alleles in 1000G (i.e. 2548samples*2alleles)
    return float(haplotype) / 5096


print("READING INPUT FILES")

# READ INPUT FILES
crispritzResultFile = sys.argv[1]  # file with CRISPRitz results
empiricalResults = sys.argv[2]  # file with empirical results
annotationFile = sys.argv[3]  # annotation file used to find nearest gene
guideFile = sys.argv[4]  # real guide used in the search
outputDir = sys.argv[5]  # output file name
check = sys.argv[6].upper()  # check if user input a GENCODE annotation file
genomeRelease = str(sys.argv[7]).strip()  # genome used in the search phase
# directory of vcf to perform haplotype count with more than one SVs in single target
vcfFileDirectory = sys.argv[8]
base_start = sys.argv[9]
base_end = sys.argv[10]
base_set = sys.argv[11]

colnames = sys.argv[12]  # annotation colnames

# OPEN INPUT FILES AND PREPARE OUTPUT FILE
# crispritz results file open
inCrispritzResults = open(crispritzResultFile, "r")
# empirical-seq data open
inEmpiricalResults = open(empiricalResults, "r").readlines()
# annotation data open
inAnnotationFile = open(annotationFile, "r")
# guide file
# realguide = open(guideFile, 'r').readlines()
# open outputDir to write results
originFileName = crispritzResultFile.split("/")
originFileName = originFileName[len(originFileName) - 1]
outFile = open(outputDir + originFileName + ".integrated_results.tsv", "w")
outFile_name = outputDir + originFileName + ".integrated_results.tsv"

# list file in vcf directory
checkVCF = False
if "vuota" not in vcfFileDirectory:
    checkVCF = True
    vcfList = glob.glob(vcfFileDirectory + "/*.gz")
    redirectFile = open(outputDir + originFileName + ".redirectFile.out", "w")
    for vcfFile in vcfList:
        # create tabix index if not already done
        subprocess.run(["tabix", str(vcfFile)], stderr=redirectFile)
    redirectFile.close()


empiricalTree = IntervalTree()
empiricalList = []
genomeDict = {}
empiricalDict = {}
valueDict = {}
# check if a personal annotation is used in the search and preserve the column
check_personal_existence = False

saveDict = {
    "Spacer+PAM": "NA",
    "Chromosome": "NA",
    "Start_coordinate_(highest_CFD)": "NA",
    "Strand_(highest_CFD)": "NA",
    "Aligned_spacer+PAM_(highest_CFD)": "NA",
    "Aligned_protospacer+PAM_REF_(highest_CFD)": "NA",
    "Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA",
    "PAM_(highest_CFD)": "NA",
    "Mismatches_(highest_CFD)": "NA",
    "Bulges_(highest_CFD)": "NA",
    "Mismatches+bulges_(highest_CFD)": "NA",
    "Seed_mismatches+bulges_REF_(highest_CFD)": "NA",
    "Non_seed_mismatches+bulges_REF_(highest_CFD)": "NA",
    "Seed_mismatches+bulges_ALT_(highest_CFD)": "NA",
    "Non_seed_mismatches+bulges_ALT_(highest_CFD)": "NA",
    "Bulge_type_(highest_CFD)": "NA",
    "REF/ALT_origin_(highest_CFD)": "NA",
    "PAM_creation_(highest_CFD)": "NA",
    "CFD_score_(highest_CFD)": "NA",
    "CFD_score_REF_(highest_CFD)": "NA",
    "CFD_score_ALT_(highest_CFD)": "NA",
    "CFD_risk_score_(highest_CFD)": "NA",
    "Variant_info_spacer+PAM_(highest_CFD)": "NA",
    "Variant_info_genome_(highest_CFD)": "NA",
    "Variant_MAF_(highest_CFD)": "NA",
    "Variant_rsID_(highest_CFD)": "NA",
    "Variant_samples_(highest_CFD)": "NA",
    "Not_found_in_REF": "NA",
    "Other_motifs": "NA",
    "Start_coordinate_(fewest_mm+b)": "NA",
    "Strand_(fewest_mm+b)": "NA",
    "Aligned_spacer+PAM_(fewest_mm+b)": "NA",
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)": "NA",
    "Aligned_protospacer+PAM_ALT_(fewest_mm+b)": "NA",
    "PAM_(fewest_mm+b)": "NA",
    "Mismatches_(fewest_mm+b)": "NA",
    "Bulges_(fewest_mm+b)": "NA",
    "Mismatches+bulges_(fewest_mm+b)": "NA",
    "Seed_mismatches+bulges_REF_(fewest_mm+b)": "NA",
    "Non_seed_mismatches+bulges_REF_(fewest_mm+b)": "NA",
    "Seed_mismatches+bulges_ALT_(fewest_mm+b)": "NA",
    "Non_seed_mismatches+bulges_ALT_(fewest_mm+b)": "NA",
    "Bulge_type_(fewest_mm+b)": "NA",
    "REF/ALT_origin_(fewest_mm+b)": "NA",
    "PAM_creation_(fewest_mm+b)": "NA",
    "CFD_score_(fewest_mm+b)": "NA",
    "CFD_score_REF_(fewest_mm+b)": "NA",
    "CFD_score_ALT_(fewest_mm+b)": "NA",
    "CFD_risk_score_(fewest_mm+b)": "NA",
    "Variant_info_spacer+PAM_(fewest_mm+b)": "NA",
    "Variant_info_genome_(fewest_mm+b)": "NA",
    "Variant_MAF_(fewest_mm+b)": "NA",
    "Variant_rsID_(fewest_mm+b)": "NA",
    "Variant_samples_(fewest_mm+b)": "NA",
    "Start_coordinate_(highest_CRISTA)": "NA",
    "Strand_(highest_CRISTA)": "NA",
    "Aligned_spacer+PAM_(highest_CRISTA)": "NA",
    "Aligned_protospacer+PAM_REF_(highest_CRISTA)": "NA",
    "Aligned_protospacer+PAM_ALT_(highest_CRISTA)": "NA",
    "PAM_(highest_CRISTA)": "NA",
    "Mismatches_(highest_CRISTA)": "NA",
    "Bulges_(highest_CRISTA)": "NA",
    "Mismatches+bulges_(highest_CRISTA)": "NA",
    "Seed_mismatches+bulges_REF_(highest_CRISTA)": "NA",
    "Non_seed_mismatches+bulges_REF_(highest_CRISTA)": "NA",
    "Seed_mismatches+bulges_ALT_(highest_CRISTA)": "NA",
    "Non_seed_mismatches+bulges_ALT_(highest_CRISTA)": "NA",
    "Bulge_type_(highest_CRISTA)": "NA",
    "REF/ALT_origin_(highest_CRISTA)": "NA",
    "PAM_creation_(highest_CRISTA)": "NA",
    "CRISTA_score_(highest_CRISTA)": "NA",
    "CRISTA_score_REF_(highest_CRISTA)": "NA",
    "CRISTA_score_ALT_(highest_CRISTA)": "NA",
    "CRISTA_risk_score_(highest_CRISTA)": "NA",
    "Variant_info_spacer+PAM_(highest_CRISTA)": "NA",
    "Variant_info_genome_(highest_CRISTA)": "NA",
    "Variant_MAF_(highest_CRISTA)": "NA",
    "Variant_rsID_(highest_CRISTA)": "NA",
    "Variant_samples_(highest_CRISTA)": "NA",
    "Annotation_GENCODE": "NA",
    "Annotation_closest_gene_name": "NA",
    "Annotation_closest_gene_ID": "NA",
    "Annotation_closest_gene_distance_(kb)": "NA",
}
#     "Annotation_ENCODE": "NA",
#     "Annotation_DHS": "NA",
#     "Annotation_personal": "NA",
#     "Susceptible_to_ABE": "NA",
#     "Susceptible_to_CBE": "NA",
#     "Susceptible_to_GBE": "NA",
#     "Susceptible_to_TBE": "NA",
# }
for cname in colnames.split(","):
    saveDict[f"Annotation_{cname}"] = "NA"

# add base editing columns
bediting_cnames = ["Susceptible_to_ABE", "Susceptible_to_CBE", "Susceptible_to_GBE", "Susceptible_to_TBE"]
for cname in bediting_cnames:
    saveDict[cname] = "NA"

print(list(saveDict.keys()))

start_time = time.time()

print("CREATING INTERVAL TREES")

for count, line in enumerate(inEmpiricalResults):
    empList = line.strip().split("\t")
    empList = [elem.strip() for elem in empList]
    empList.append(count)
    # adding empirical data to the tree
    empiricalTree[int(empList[1]) : int(empList[2])] = empList
    # generate temp empirical dict to save empirical information per row
    empiricalDict[str(empList[4])] = 50
    # value of empirical row, keeping info about mm+bul for each empirical origin (seq data, in-silico, vitro, ecc)
    valueDict[str(empList[5])] = "NA"
    # update save dict with user-defined names from empirical data
    saveDict[str(empList[5])] = "NA"
    # newkey = str(empList[5])+'_mm+bul'
    # saveDict[newkey] = 'NA'

# writing header in file
save = ""
save += "\t".join(list(saveDict.keys()))
save += "\n"
outFile.write(save)

print("INTEGRATING RESULTS")

if "#" in inCrispritzResults.readline():
    print("SKIP HEADER")
else:
    inCrispritzResults.seek(0)

for nline, line in enumerate(inCrispritzResults):
    target = line.strip().split("\t")
    # file annotation reported after gencode association with gene
    annotationLine = inAnnotationFile.readline().strip().split("\t")
    # reset dict to save new line
    for key in saveDict:
        saveDict[key] = "NA"
    # reset empiricaldict
    for key in empiricalDict:
        empiricalDict[key] = 50
        valueDict[key] = "NA"

    if "NA" not in annotationLine and check == "TRUE":
        for elem in annotationLine:
            if "gene_id" in elem:
                temp = elem.strip().split(";")
                for name in temp:
                    if "gene_id" in name:
                        saveDict["Annotation_closest_gene_ID"] = name.strip().split(
                            "="
                        )[1]
                    if "gene_name" in name:
                        saveDict["Annotation_closest_gene_name"] = name.strip().split(
                            "="
                        )[1]
        saveDict["Annotation_closest_gene_distance_(kb)"] = str(
            float(annotationLine[len(annotationLine) - 1]) / 1000
        )
        if float(annotationLine[len(annotationLine) - 1]) != 0:
            saveDict["Annotation_GENCODE"] = "intergenic"
        else:
            saveDict["Annotation_GENCODE"] = annotationLine[11].strip()

    variantList = ["NA"]
    if str(target[18]) != "NA":  # check if target has variants reported (CFD)
        variantList = str(target[18]).strip().split(",")
        # generate variant nucleotide reverse complemented always in positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split("_")
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if "+" in str(target[7]):  # strand
                # var pos is equal to pos_of_variant-real_position
                var_pos = int(split[1]) - int(target[5])  # real_position
                if var_pos < 1:
                    var_pos = 1
                variantList[count] = str(var_pos) + str(split[2]) + ">" + str(split[3])
            else:
                firstcomp = ""
                secondcomp = ""
                var_pos = int(split[1]) - int(target[5])
                if var_pos < 1:
                    var_pos = len(target[1])
                else:
                    # if var pos non negative, count the position in reverse strand so starting from end of the sequence that is been reversed complemented
                    if split_second_len < split_one_len:
                        correction = split_one_len - split_second_len
                    var_pos = abs(
                        int(split[1]) - int(target[5]) - 1 - len(target[1]) + correction
                    )
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = (
                    str(var_pos)
                    + "".join(reversed(firstcomp))
                    + ">"
                    + "".join(reversed(secondcomp))
                )
    variantList_highest_cfd = variantList

    variantList = ["NA"]
    if str(target[42]) != "NA":  # check if target has variants reported (MMBUL)
        variantList = str(target[42]).strip().split(",")
        var_pos = []
        # generate variant position corrected to be in the positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split("_")
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if "+" in str(target[31]):
                # var pos is equal to pos_of_variant-real_position
                var_pos = int(split[1]) - int(target[29])
                if var_pos < 1:
                    var_pos = 1
                variantList[count] = str(var_pos) + str(split[2]) + ">" + str(split[3])
            else:
                firstcomp = ""
                secondcomp = ""
                var_pos = int(split[1]) - int(target[29])
                if var_pos < 1:
                    var_pos = len(target[25])  # aligned_sgRNA
                else:
                    # if var pos non negative, count the position in reverse strand so starting from end of the sequence that is been reversed complemented
                    if split_second_len < split_one_len:
                        correction = split_one_len - split_second_len
                    var_pos = abs(
                        int(split[1])
                        - int(target[29])
                        - 1
                        - len(target[25])
                        + correction
                    )
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = (
                    str(var_pos)
                    + "".join(reversed(firstcomp))
                    + ">"
                    + "".join(reversed(secondcomp))
                )
    variantList_fewest_mm_b = variantList

    variantList = ["NA"]
    if str(target[66]) != "NA":  # check if target has variants reported (CRISTA)
        variantList = str(target[66]).strip().split(",")
        var_pos = []
        # generate variant position corrected to be in the positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split("_")
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if "+" in str(target[55]):
                # var pos is equal to pos_of_variant-real_position
                var_pos = int(split[1]) - int(target[53])
                if var_pos < 1:
                    var_pos = 1
                variantList[count] = str(var_pos) + str(split[2]) + ">" + str(split[3])
            else:
                firstcomp = ""
                secondcomp = ""
                var_pos = int(split[1]) - int(target[53])
                if var_pos < 1:
                    var_pos = len(target[49])  # aligned_sgRNA
                else:
                    # if var pos non negative, count the position in reverse strand so starting from end of the sequence that is been reversed complemented
                    if split_second_len < split_one_len:
                        correction = split_one_len - split_second_len
                    var_pos = abs(
                        int(split[1])
                        - int(target[53])
                        - 1
                        - len(target[49])
                        + correction
                    )
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = (
                    str(var_pos)
                    + "".join(reversed(firstcomp))
                    + ">"
                    + "".join(reversed(secondcomp))
                )
    variantList_highest_crista = variantList

    saveDict["Spacer+PAM"] = str(target[15])
    saveDict["Strand_(highest_CFD)"] = str(target[7])
    saveDict["Chromosome"] = str(target[4])
    saveDict["Start_coordinate_(highest_CFD)"] = str(target[5])
    saveDict["Aligned_spacer+PAM_(highest_CFD)"] = str(target[1])
    saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"] = str(target[3])
    saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = str(target[2])
    saveDict["PAM_(highest_CFD)"] = "NA"
    saveDict["Mismatches_(highest_CFD)"] = str(target[8])
    saveDict["Bulges_(highest_CFD)"] = str(target[9])
    saveDict["Mismatches+bulges_(highest_CFD)"] = str(target[10])
    saveDict["Bulge_type_(highest_CFD)"] = str(target[0])
    saveDict["REF/ALT_origin_(highest_CFD)"] = (
        "ref" if str(target[13]) == "NA" else "alt"
    )

    # switch ref with alt alignment if target is ref
    if saveDict["REF/ALT_origin_(highest_CFD)"] == "ref":
        saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"] = str(target[2])
        saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = "NA"

    saveDict["PAM_creation_(highest_CFD)"] = str(target[11])
    saveDict["CFD_score_(highest_CFD)"] = (
        str(target[20]) if float(target[20]) > float(target[21]) else str(target[21])
    )
    saveDict["CFD_score_REF_(highest_CFD)"] = str(target[21])
    saveDict["CFD_score_ALT_(highest_CFD)"] = str(target[20])
    saveDict["CFD_risk_score_(highest_CFD)"] = str(target[22])
    saveDict["Variant_info_spacer+PAM_(highest_CFD)"] = ",".join(
        variantList_highest_cfd
    )
    saveDict["Variant_info_genome_(highest_CFD)"] = str(target[18])

    # remove 0 MAF
    maf_list = list()
    for elem in target[17].strip().split(","):
        if elem != "NA":
            if float(elem) == 0:
                maf_list.append(str(0.00001))
            else:
                maf_list.append(str(elem))
        else:
            maf_list.append("NA")

    saveDict["Variant_MAF_(highest_CFD)"] = ",".join(maf_list)

    saveDict["Variant_rsID_(highest_CFD)"] = (
        "NA" if str(target[16]) == "." else str(target[16])
    )
    saveDict["Variant_samples_(highest_CFD)"] = str(target[13])
    saveDict["Not_found_in_REF"] = "y" if str(target[12]) == "y" else "NA"
    saveDict["Other_motifs"] = str(target[19])
    # savedict for fewest_mm+b
    saveDict["Strand_(fewest_mm+b)"] = str(target[31])
    saveDict["Start_coordinate_(fewest_mm+b)"] = str(target[29])
    saveDict["Aligned_spacer+PAM_(fewest_mm+b)"] = str(target[25])
    saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"] = str(target[27])
    saveDict["Aligned_protospacer+PAM_ALT_(fewest_mm+b)"] = str(target[26])
    saveDict["PAM_(fewest_mm+b)"] = "NA"
    saveDict["Mismatches_(fewest_mm+b)"] = str(target[32])
    saveDict["Bulges_(fewest_mm+b)"] = str(target[33])
    saveDict["Mismatches+bulges_(fewest_mm+b)"] = str(target[34])
    saveDict["Bulge_type_(fewest_mm+b)"] = str(target[24])
    saveDict["REF/ALT_origin_(fewest_mm+b)"] = (
        "ref" if str(target[37]) == "NA" else "alt"
    )

    # switch ref with alt alignment if target is ref
    if saveDict["REF/ALT_origin_(fewest_mm+b)"] == "ref":
        saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"] = str(target[26])
        saveDict["Aligned_protospacer+PAM_ALT_(fewest_mm+b)"] = "NA"

    saveDict["PAM_creation_(fewest_mm+b)"] = str(target[35])
    saveDict["CFD_score_(fewest_mm+b)"] = (
        str(target[44]) if float(target[44]) > float(target[45]) else str(target[45])
    )
    saveDict["CFD_score_REF_(fewest_mm+b)"] = str(target[45])
    saveDict["CFD_score_ALT_(fewest_mm+b)"] = str(target[44])
    saveDict["CFD_risk_score_(fewest_mm+b)"] = str(target[46])
    saveDict["Variant_info_spacer+PAM_(fewest_mm+b)"] = ",".join(
        variantList_fewest_mm_b
    )
    saveDict["Variant_info_genome_(fewest_mm+b)"] = str(target[42])

    maf_list = list()
    for elem in target[41].strip().split(","):
        if elem != "NA":
            if float(elem) == 0:
                maf_list.append(str(0.00001))
            else:
                maf_list.append(str(elem))
        else:
            maf_list.append("NA")

    saveDict["Variant_MAF_(fewest_mm+b)"] = ",".join(maf_list)

    saveDict["Variant_rsID_(fewest_mm+b)"] = (
        "NA" if str(target[40]) == "." else str(target[40])
    )
    saveDict["Variant_samples_(fewest_mm+b)"] = str(target[37])
    # savedict for highestCRISTA
    saveDict["Strand_(highest_CRISTA)"] = str(target[55])
    saveDict["Start_coordinate_(highest_CRISTA)"] = str(target[53])
    saveDict["Aligned_spacer+PAM_(highest_CRISTA)"] = str(target[49])
    saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"] = str(target[51])
    saveDict["Aligned_protospacer+PAM_ALT_(highest_CRISTA)"] = str(target[50])
    saveDict["PAM_(highest_CRISTA)"] = "NA"
    saveDict["Mismatches_(highest_CRISTA)"] = str(target[56])
    saveDict["Bulges_(highest_CRISTA)"] = str(target[57])
    saveDict["Mismatches+bulges_(highest_CRISTA)"] = str(target[58])
    saveDict["Bulge_type_(highest_CRISTA)"] = str(target[48])
    saveDict["REF/ALT_origin_(highest_CRISTA)"] = (
        "ref" if str(target[61]) == "NA" else "alt"
    )

    # switch ref with alt alignment if target is ref
    if saveDict["REF/ALT_origin_(highest_CRISTA)"] == "ref":
        saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"] = str(target[50])
        saveDict["Aligned_protospacer+PAM_ALT_(highest_CRISTA)"] = "NA"

    saveDict["PAM_creation_(highest_CRISTA)"] = str(target[59])
    saveDict["CRISTA_score_(highest_CRISTA)"] = (
        str(target[68]) if float(target[68]) > float(target[69]) else str(target[69])
    )
    saveDict["CRISTA_score_REF_(highest_CRISTA)"] = str(target[69])
    saveDict["CRISTA_score_ALT_(highest_CRISTA)"] = str(target[68])
    saveDict["CRISTA_risk_score_(highest_CRISTA)"] = str(target[70])
    saveDict["Variant_info_spacer+PAM_(highest_CRISTA)"] = ",".join(
        variantList_highest_crista
    )
    saveDict["Variant_info_genome_(highest_CRISTA)"] = str(target[66])

    maf_list = list()
    for elem in target[65].strip().split(","):
        if elem != "NA":
            if float(elem) == 0:
                maf_list.append(str(0.00001))
            else:
                maf_list.append(str(elem))
        else:
            maf_list.append("NA")

    saveDict["Variant_MAF_(highest_CRISTA)"] = ",".join(maf_list)

    saveDict["Variant_rsID_(highest_CRISTA)"] = (
        "NA" if str(target[64]) == "." else str(target[64])
    )
    saveDict["Variant_samples_(highest_CRISTA)"] = str(target[61])

    # change_alt_ref_highest_cfd = False
    # if saveDict['Not_found_in_REF'] == 'NA' and saveDict['REF/ALT_origin_(highest_CFD)'] == 'alt' and saveDict['CFD_score_REF_(highest_CFD)'] != '-1.0' and saveDict['CFD_score_REF_(highest_CFD)'] == saveDict['CFD_score_ALT_(highest_CFD)']:
    #     change_alt_ref_highest_cfd = True

    # change_alt_ref_highest_crista = False
    # if saveDict['Not_found_in_REF'] == 'NA' and saveDict['REF/ALT_origin_(highest_CRISTA)'] == 'alt' and saveDict['CRISTA_score_REF_(highest_CRISTA)'] != '-1.0' and saveDict['CRISTA_score_REF_(highest_CRISTA)'] == saveDict['CRISTA_score_ALT_(highest_CRISTA)']:
    #     change_alt_ref_highest_crista = True

    # if change_alt_ref_highest_cfd:
    #     saveDict['Aligned_protospacer+PAM_ALT_(highest_CFD)'] = 'NA'
    #     saveDict['REF/ALT_origin_(highest_CFD)'] = 'ref'
    #     saveDict['Variant_info_spacer+PAM_(highest_CFD)'] = 'NA'
    #     saveDict['Variant_info_genome_(highest_CFD)'] = 'NA'
    #     saveDict['Variant_MAF_(highest_CFD)'] = 'NA'
    #     saveDict['Variant_rsID_(highest_CFD)'] = 'NA'
    #     saveDict['Variant_samples_(highest_CFD)'] = 'NA'
    #     saveDict['PAM_creation_(highest_CFD)'] = 'NA'

    # if change_alt_ref_highest_crista:
    #     saveDict['Aligned_protospacer+PAM_ALT_(highest_CRISTA)'] = 'NA'
    #     saveDict['REF/ALT_origin_(highest_CRISTA)'] = 'ref'
    #     saveDict['Variant_info_spacer+PAM_(highest_CRISTA)'] = 'NA'
    #     saveDict['Variant_info_genome_(highest_CRISTA)'] = 'NA'
    #     saveDict['Variant_MAF_(highest_CRISTA)'] = 'NA'
    #     saveDict['Variant_rsID_(highest_CRISTA)'] = 'NA'
    #     saveDict['Variant_samples_(highest_CRISTA)'] = 'NA'
    #     saveDict['PAM_creation_(highest_CRISTA)'] = 'NA'

    # check how long is the pam counting Ns in the guide
    count_N_in_guide = str(target[15]).count("N")
    # check if pam is at start of the sequence
    if str(target[15])[0] == "N":
        pam_at_start = True
    else:
        pam_at_start = False

    # count seed and non-seed for highestCFD
    if pam_at_start:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"][
            count_N_in_guide:
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"][
            count_N_in_guide:
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(highest_CFD)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(highest_CFD)"][count_N_in_guide:]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]
        non_seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        non_seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
    else:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"][
            : len(saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"])
            - count_N_in_guide
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"][
            : len(saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"])
            - count_N_in_guide
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(highest_CFD)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(highest_CFD)"][
                : len(saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"])
                - count_N_in_guide
            ]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
        non_seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        non_seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]

    if saveDict["REF/ALT_origin_(highest_CFD)"] == "ref":
        seed_alt == "NA"
        non_seed_alt == "NA"

    # count_seed_ref, count_seed_alt, count_non_seed_ref, count_non_seed_alt
    seed_list = seed_processing(seed_ref, seed_alt, non_seed_ref, non_seed_alt)

    saveDict["Seed_mismatches+bulges_REF_(highest_CFD)"] = str(
        int(seed_list[0]) + dna_bulge_count_seed
    )
    saveDict["Seed_mismatches+bulges_ALT_(highest_CFD)"] = str(
        int(seed_list[1]) + dna_bulge_count_seed
    )
    saveDict["Non_seed_mismatches+bulges_REF_(highest_CFD)"] = str(
        int(seed_list[2]) + dna_bulge_count_non_seed
    )
    saveDict["Non_seed_mismatches+bulges_ALT_(highest_CFD)"] = str(
        int(seed_list[3]) + dna_bulge_count_non_seed
    )

    # count seed and non-seed for fewest_mm+b
    if pam_at_start:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"][
            count_N_in_guide:
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(fewest_mm+b)"][
            count_N_in_guide:
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(fewest_mm+b)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(fewest_mm+b)"][count_N_in_guide:]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]
        non_seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        non_seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
    else:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"][
            : len(saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"])
            - count_N_in_guide
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(fewest_mm+b)"][
            : len(saveDict["Aligned_protospacer+PAM_ALT_(fewest_mm+b)"])
            - count_N_in_guide
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(fewest_mm+b)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(fewest_mm+b)"][
                : len(saveDict["Aligned_protospacer+PAM_REF_(fewest_mm+b)"])
                - count_N_in_guide
            ]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
        non_seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        non_seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]

    if saveDict["REF/ALT_origin_(fewest_mm+b)"] == "ref":
        seed_alt == "NA"
        non_seed_alt == "NA"

    # count_seed_ref, count_seed_alt, count_non_seed_ref, count_non_seed_alt
    seed_list = seed_processing(seed_ref, seed_alt, non_seed_ref, non_seed_alt)

    saveDict["Seed_mismatches+bulges_REF_(fewest_mm+b)"] = str(seed_list[0])
    saveDict["Seed_mismatches+bulges_ALT_(fewest_mm+b)"] = str(seed_list[1])
    saveDict["Non_seed_mismatches+bulges_REF_(fewest_mm+b)"] = str(seed_list[2])
    saveDict["Non_seed_mismatches+bulges_ALT_(fewest_mm+b)"] = str(seed_list[3])

    # count seed and non-seed for highest_CRISTA
    if pam_at_start:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"][
            count_N_in_guide:
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(highest_CRISTA)"][
            count_N_in_guide:
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(highest_CRISTA)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(highest_CRISTA)"][count_N_in_guide:]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]
        non_seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        non_seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
    else:
        real_target_ref = saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"][
            : len(saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"])
            - count_N_in_guide
        ]
        real_target_alt = saveDict["Aligned_protospacer+PAM_ALT_(highest_CRISTA)"][
            : len(saveDict["Aligned_protospacer+PAM_ALT_(highest_CRISTA)"])
            - count_N_in_guide
        ]
        dna_bulge_count_seed = 0
        dna_bulge_count_non_seed = 0
        if "DNA" in saveDict["Bulge_type_(highest_CRISTA)"]:
            dna_seq = saveDict["Aligned_spacer+PAM_(highest_CRISTA)"][
                : len(saveDict["Aligned_protospacer+PAM_REF_(highest_CRISTA)"])
                - count_N_in_guide
            ]
            for pos, nt in enumerate(dna_seq):
                if nt == "-" and pos >= int(len(real_target_ref) / 2):
                    dna_bulge_count_seed += 1
                elif nt == "-" and pos < int(len(real_target_ref) / 2):
                    dna_bulge_count_non_seed += 1
        seed_ref = real_target_ref[int(len(real_target_ref) / 2) :]
        seed_alt = real_target_alt[int(len(real_target_alt) / 2) :]
        non_seed_ref = real_target_ref[: int(len(real_target_ref) / 2)]
        non_seed_alt = real_target_alt[: int(len(real_target_alt) / 2)]

    if saveDict["REF/ALT_origin_(highest_CRISTA)"] == "ref":
        seed_alt == "NA"
        non_seed_alt == "NA"

    # count_seed_ref, count_seed_alt, count_non_seed_ref, count_non_seed_alt
    seed_list = seed_processing(seed_ref, seed_alt, non_seed_ref, non_seed_alt)

    saveDict["Seed_mismatches+bulges_REF_(highest_CRISTA)"] = str(
        int(seed_list[0]) + dna_bulge_count_seed
    )
    saveDict["Seed_mismatches+bulges_ALT_(highest_CRISTA)"] = str(
        int(seed_list[1]) + dna_bulge_count_seed
    )
    saveDict["Non_seed_mismatches+bulges_REF_(highest_CRISTA)"] = str(
        int(seed_list[2]) + dna_bulge_count_non_seed
    )
    saveDict["Non_seed_mismatches+bulges_ALT_(highest_CRISTA)"] = str(
        int(seed_list[3]) + dna_bulge_count_non_seed
    )

    # reset count seed and non-seed for ALT if target is REF
    if saveDict["REF/ALT_origin_(highest_CFD)"] == "ref":
        saveDict["Seed_mismatches+bulges_ALT_(highest_CFD)"] = "0"
        saveDict["Non_seed_mismatches+bulges_ALT_(highest_CFD)"] = "0"

    if saveDict["REF/ALT_origin_(fewest_mm+b)"] == "ref":
        saveDict["Seed_mismatches+bulges_ALT_(fewest_mm+b)"] = "0"
        saveDict["Non_seed_mismatches+bulges_ALT_(fewest_mm+b)"] = "0"

    if saveDict["REF/ALT_origin_(highest_CRISTA)"] == "ref":
        saveDict["Seed_mismatches+bulges_ALT_(highest_CRISTA)"] = "0"
        saveDict["Non_seed_mismatches+bulges_ALT_(highest_CRISTA)"] = "0"

    # extract PAM sequence
    # highestCFD
    if saveDict["REF/ALT_origin_(highest_CFD)"] == "ref":
        if pam_at_start:  # save pam sequence extracting directly from the ref sequence
            saveDict["PAM_(highest_CFD)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(highest_CFD)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(highest_CFD)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(highest_CFD)"
            ][-count_N_in_guide:]
    else:
        if pam_at_start:  # save pam sequence extracting directly from the var sequence
            saveDict["PAM_(highest_CFD)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(highest_CFD)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(highest_CFD)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(highest_CFD)"
            ][-count_N_in_guide:]

    # fewest_mm+b
    if saveDict["REF/ALT_origin_(fewest_mm+b)"] == "ref":
        if pam_at_start:  # save pam sequence extracting directly from the ref sequence
            saveDict["PAM_(fewest_mm+b)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(fewest_mm+b)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(fewest_mm+b)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(fewest_mm+b)"
            ][-count_N_in_guide:]
    else:
        if pam_at_start:  # save pam sequence extracting directly from the var sequence
            saveDict["PAM_(fewest_mm+b)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(fewest_mm+b)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(fewest_mm+b)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(fewest_mm+b)"
            ][-count_N_in_guide:]

    # highest_CRISTA
    if saveDict["REF/ALT_origin_(highest_CRISTA)"] == "ref":
        if pam_at_start:  # save pam sequence extracting directly from the ref sequence
            saveDict["PAM_(highest_CRISTA)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(highest_CRISTA)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(highest_CRISTA)"] = saveDict[
                "Aligned_protospacer+PAM_REF_(highest_CRISTA)"
            ][-count_N_in_guide:]
    else:
        if pam_at_start:  # save pam sequence extracting directly from the var sequence
            saveDict["PAM_(highest_CRISTA)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(highest_CRISTA)"
            ][:count_N_in_guide]
        else:
            saveDict["PAM_(highest_CRISTA)"] = saveDict[
                "Aligned_protospacer+PAM_ALT_(highest_CRISTA)"
            ][-count_N_in_guide:]

    # annotate with empirical and convert _personal and _gencode annotation to better visualization
    annotation_subdict = {f"Annotation_{cname}": set() for cname in colnames.split(",")}
    annotationList = target[14].split(",")
    if len(annotationList) != 1 and annotationList[0] != "NA":
        for ann in annotationList:
            cname = ann.split("__")[-1]
            annotation_subdict[f"Annotation_{cname}"].add(ann.replace(f"__{cname}", ""))
    for cname, ann in annotation_subdict.items():
        saveDict[cname] = ",".join(ann) if ann else "NA"

    # for elem in annotationList:
    #     if "_personal" in elem:
    #         personal_annotations.add(elem.replace("_personal", ""))
    #     elif "_gencode" in elem:
    #         gencode_annotations.add(elem.replace("_gencode", ""))
    #     elif "_DHS" in elem:
    #         DHS_annotations.add(elem.replace("_DHS", ""))
    #     else:
    #         encode_annotations.add(elem)

    # if len(DHS_annotations):
    #     saveDict["Annotation_DHS"] = ",".join(DHS_annotations)

    # if len(personal_annotations):
    #     saveDict["Annotation_personal"] = ",".join(personal_annotations)
    #     check_personal_existence = True

    # if len(encode_annotations):
    #     saveDict["Annotation_ENCODE"] = ",".join(encode_annotations)

    # if len(gencode_annotations):
    #     saveDict["Annotation_GENCODE"] = ",".join(gencode_annotations)

    foundEmpirical = sorted(empiricalTree[int(target[6]) - 4 : int(target[6]) + 4])

    for found in range(0, len(foundEmpirical)):
        empirical = foundEmpirical[found].data
        if str(saveDict["Chromosome"]) == str(empirical[0]):
            empiricalList.append(empirical[7])
            valueDict[str(empirical[4])] = empirical[5]
            empiricalDict[str(empirical[4])] = int(empirical[3])

    for key in empiricalDict:
        if int(empiricalDict[key]) < 50:
            saveDict[key] = str(valueDict[key])
            newkey = str(key) + "_mm+bul"
            saveDict[newkey] = empiricalDict[key]

    # fixed value for base editor analysis
    start = 1
    end = 0
    try:
        # check if base window is int
        base_start = int(base_start)
        base_end = int(base_end)
    except:
        # use default value
        base_start = 1
        base_end = 0
    if base_start > 0 and base_start <= len(saveDict["Spacer+PAM"]):
        start = base_start
    if base_end > 0 and base_end <= len(saveDict["Spacer+PAM"]):
        end = base_end

    # print('basestart:',base_start,' baseend:',base_end)
    if saveDict["REF/ALT_origin_(highest_CFD)"] == "alt":
        value = base_counting(
            saveDict["Aligned_protospacer+PAM_ALT_(highest_CFD)"], start - 1, end
        )
        # check count A
        if value[0]:
            saveDict["Susceptible_to_ABE"] = "y"
        # check count C
        if value[1]:
            saveDict["Susceptible_to_CBE"] = "y"
        # check count G
        if value[2]:
            saveDict["Susceptible_to_GBE"] = "y"
        # check count T
        if value[3]:
            saveDict["Susceptible_to_TBE"] = "y"
    else:
        value = base_counting(
            saveDict["Aligned_protospacer+PAM_REF_(highest_CFD)"], start - 1, end
        )
        # check count A
        if value[0]:
            saveDict["Susceptible_to_ABE"] = "y"
        # check count C
        if value[1]:
            saveDict["Susceptible_to_CBE"] = "y"
        # check count G
        if value[2]:
            saveDict["Susceptible_to_GBE"] = "y"
        # check count T
        if value[3]:
            saveDict["Susceptible_to_TBE"] = "y"

    save = "\t".join(list(saveDict.values()))
    save += "\n"
    outFile.write(save)

# close integrated file
outFile.close()


# remove unused columns from final integrated file
columns_of_dict = list(saveDict.keys())
columns_to_remove = list()
# added 1 to account for bash counting starting from 1
# columns_to_remove.append(columns_of_dict.index("Annotation_personal") + 1)
columns_to_remove.append(columns_of_dict.index("Susceptible_to_ABE") + 1)
columns_to_remove.append(columns_of_dict.index("Susceptible_to_CBE") + 1)
columns_to_remove.append(columns_of_dict.index("Susceptible_to_GBE") + 1)
columns_to_remove.append(columns_of_dict.index("Susceptible_to_TBE") + 1)

# check columns to keep
# if check_personal_existence:
#     columns_to_remove.remove(columns_of_dict.index("Annotation_personal") + 1)

if "none" not in base_set:
    for base in base_set.strip().split(","):
        base = base.upper()
        columns_to_remove.remove(
            columns_of_dict.index("Susceptible_to_" + base + "BE") + 1
        )
# convert int list to str list for joining in cut command
output = [str(x) for x in columns_to_remove]
string_to_remove = ",".join(output)
# call cut in terminal to remove unused columns
os.system(f"cut -f {string_to_remove} --complement {outFile_name} > {outFile_name}.tmp")
os.system(f"mv {outFile_name}.tmp {outFile_name}")

print("CHECKING MISSING RESULTS")

notFoundFile = open(outputDir + originFileName + ".empirical_not_found.tsv", "w")

for count, line in enumerate(inEmpiricalResults):
    if count not in empiricalList:
        notFoundFile.write(line)

# close notfound file
notFoundFile.close()

print("INTEGRATION COMPLETED IN: %s seconds" % (time.time() - start_time))
