#!/usr/bin/env python

from operator import truediv
from posixpath import expanduser
from intervaltree import IntervalTree
import sys
import time
import glob
import subprocess
import os

print('READING INPUT FILES')

# READ INPUT FILES
crispritzResultFile = sys.argv[1]  # file with CRISPRitz results
empiricalResults = sys.argv[2]  # file with empirical results
# annotationFile = sys.argv[3]  # annotation file used to find nearest gene
# guideFile = sys.argv[4]  # real guide used in the search
outputDir = sys.argv[3]  # output file name
# check = sys.argv[6].upper()  # check if user input a GENCODE annotation file
# genomeRelease = str(sys.argv[7]).strip()  # genome used in the search phase
# directory of vcf to perform haplotype count with more than one SVs in single target
# vcfFileDirectory = sys.argv[8]

# OPEN INPUT FILES AND PREPARE OUTPUT FILE
# crispritz results file open
inCrispritzResults = open(crispritzResultFile, 'r')
# empirical-seq data open
inEmpiricalResults = open(empiricalResults, 'r').readlines()
# annotation data open
# inAnnotationFile = open(annotationFile, 'r')
# guide file
# realguide = open(guideFile, 'r').readlines()
# open outputDir to write results
originFileName = crispritzResultFile.split('/')[-1]
# originFileName = originFileName[len(originFileName)-1]
outFile = open(outputDir + '/' + originFileName +
               '.empirical_data.tsv', 'w')
# outFile_name = outputDir+originFileName+'.integrated_results.tsv'


empiricalTree = IntervalTree()
empiricalList = list()
empiricalDict = dict()

saveDict = {
    'Spacer+PAM': 'NA',
    'Chromosome': 'NA',
    'Start_coordinate_(highest_CFD)': 'NA',
    'Strand_(highest_CFD)': 'NA',
    'Aligned_spacer+PAM_(highest_CFD)': 'NA',
    'Aligned_protospacer+PAM_REF_(highest_CFD)': 'NA',
    'Aligned_protospacer+PAM_ALT_(highest_CFD)': 'NA',
    'PAM_(highest_CFD)': 'NA',
    'Mismatches_(highest_CFD)': 'NA',
    'Bulges_(highest_CFD)': 'NA',
    'Mismatches+bulges_(highest_CFD)': 'NA',
    'Bulge_type_(highest_CFD)': 'NA',
    'REF/ALT_origin_(highest_CFD)': 'NA',
    'PAM_creation_(highest_CFD)': 'NA',
    'CFD_score_(highest_CFD)': 'NA',
    'CFD_score_REF_(highest_CFD)': 'NA',
    'CFD_score_ALT_(highest_CFD)': 'NA',
    'CFD_risk_score_(highest_CFD)': 'NA',
    'Variant_info_spacer+PAM_(highest_CFD)': 'NA',
    'Variant_info_genome_(highest_CFD)': 'NA',
    'Variant_MAF_(highest_CFD)': 'NA',
    'Variant_rsID_(highest_CFD)': 'NA',
    'Variant_samples_(highest_CFD)': 'NA',
    'Not_found_in_REF': 'NA',
    'Other_motifs': 'NA',
    'Start_coordinate_(fewest_mm+b)': 'NA',
    'Strand_(fewest_mm+b)': 'NA',
    'Aligned_spacer+PAM_(fewest_mm+b)': 'NA',
    'Aligned_protospacer+PAM_REF_(fewest_mm+b)': 'NA',
    'Aligned_protospacer+PAM_ALT_(fewest_mm+b)': 'NA',
    'PAM_(fewest_mm+b)': 'NA',
    'Mismatches_(fewest_mm+b)': 'NA',
    'Bulges_(fewest_mm+b)': 'NA',
    'Mismatches+bulges_(fewest_mm+b)': 'NA',
    'Bulge_type_(fewest_mm+b)': 'NA',
    'REF/ALT_origin_(fewest_mm+b)': 'NA',
    'PAM_creation_(fewest_mm+b)': 'NA',
    'CFD_score_(fewest_mm+b)': 'NA',
    'CFD_score_REF_(fewest_mm+b)': 'NA',
    'CFD_score_ALT_(fewest_mm+b)': 'NA',
    'CFD_risk_score_(fewest_mm+b)': 'NA',
    'Variant_info_spacer+PAM_(fewest_mm+b)': 'NA',
    'Variant_info_genome_(fewest_mm+b)': 'NA',
    'Variant_MAF_(fewest_mm+b)': 'NA',
    'Variant_rsID_(fewest_mm+b)': 'NA',
    'Variant_samples_(fewest_mm+b)': 'NA',
    'Annotation_GENCODE': 'NA',
    'Annotation_closest_gene_name': 'NA',
    'Annotation_closest_gene_ID': 'NA',
    'Annotation_closest_gene_distance_(kb)': 'NA',
    'Annotation_ENCODE': 'NA',
    'Annotation_personal': 'NA'
}


start_time = time.time()

print('CREATING INTERVAL TREES')

for count, line in enumerate(inEmpiricalResults):
    empList = line.strip().split('\t')
    empList = [elem.strip() for elem in empList]
    # example row for empirical list
    # mand  mand        mand       mand mand       mand  optional
    # chr10	33753323	33753346	4	CIRCLEseq	OT1 aTtACAGcTGCaTTTATCACAGG
    empList.append(count)
    # adding empirical data to the tree
    empiricalTree[int(empList[1]):int(empList[2])] = empList
    # to save header
    saveDict[str(empList[4])] = 'NA'
    newkey = str(empList[4])+'_mm+bul'
    saveDict[newkey] = 'NA'
    # to save data of empirical
    empiricalDict[str(empList[4])] = 'NA'
    empiricalDict[newkey] = 'NA'

# writing header in file
save = ''
save += '\t'.join(list(saveDict.keys()))
save += '\n'
outFile.write(save)

print('INTEGRATING RESULTS')

if 'Chromosome' in inCrispritzResults.readline():
    print('SKIP HEADER')
else:
    inCrispritzResults.seek(0)

for nline, line in enumerate(inCrispritzResults):
    target = line.strip().split('\t')

    # for key in saveDict:
    #     saveDict[key] = 'NA'

    for key in empiricalDict:
        empiricalDict[key] = 'NA'
    #     valueDict[key] = 'NA'
    # lowestEmpirical = 100

    # read chr from target line
    # saveDict['Chromosome'] = target[1]

    # search empirical target using in-silico target position with a window
    foundEmpirical = sorted(empiricalTree[int(target[2])-4:int(target[2])+4])

    # search in the list of found empirical to extract data with same chr and save them
    for found in range(0, len(foundEmpirical)):
        empirical = foundEmpirical[found].data
        if str(target[1]) == str(empirical[0]):
            empiricalList.append(empirical[-1])
            empiricalDict[str(empirical[4])] = str(empirical[5])
            empiricalDict[str(empirical[4])+'_mm+bul'] = str(empirical[3])
            # valueDict[str(empirical[4])] = empirical[5]
            # empiricalDict[str(empirical[4])] = int(empirical[3])

    # update the empirical dict with found targets
    # for key in empiricalDict:
    #     if int(empiricalDict[key]) < 50:
    #         saveDict[key] = str(valueDict[key])
    #         newkey = str(key)+'_mm+bul'
    #         saveDict[newkey] = empiricalDict[key]
            # if int(empiricalDict[key]) < lowestEmpirical:
            # saveDict['lowest_empirical'] = str(empiricalDict[key])

    # save = str()
    # for key in saveDict:
    #     save += str(saveDict[key])+'\t'
    # save += '\n'

    # save row of target with empirical data
    outFile.write('\t'.join(target) + '\t' +
                  '\t'.join(empiricalDict.values())+'\n')

print('CHECKING MISSING RESULTS')

notFoundFile = open(outputDir + '/' + originFileName +
                    '.empirical_not_found.tsv', 'w')

for count, line in enumerate(inEmpiricalResults):
    if count not in empiricalList:
        notFoundFile.write(line)

print("INTEGRATION COMPLETED IN: %s seconds" % (time.time() - start_time))
