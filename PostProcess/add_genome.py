#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 19:23:05 2020

@author: whitebreeze
"""
import os
from os.path import isfile, isdir, join  # for getting directories
import sys
from shutil import copy, rmtree

fileDir = sys.argv[1]
appDir = sys.argv[2]
genomeDir = sys.argv[3]
PAMFile = sys.argv[4]
annotationFile = sys.argv[5]
bMax = sys.argv[6]
VCFDir = sys.argv[7]
sampleFile = sys.argv[8]
enrichedName = sys.argv[9]

is_vcf_given = False
if VCFDir != "None":
    is_vcf_given = True

cleanName = ''.join([char for char in os.path.basename(
    os.path.normpath(genomeDir)) if char != "+"])


os.chdir(fileDir)
with open("Genomes/aggiunta_nuovo_genoma.txt", 'w') as log:
    log.write("0 Copy-files")

if not os.path.exists(fileDir+"/Genomes/"+cleanName+"_ref"):
    os.mkdir(fileDir+"/Genomes/"+cleanName+"_ref")

for item in os.listdir(genomeDir+"/"):
    if not os.path.exists(fileDir+"/Genomes/"+cleanName+"_ref/"+item):
        copy(genomeDir+"/"+item, fileDir+"/Genomes/"+cleanName+"_ref/"+item)

with open(PAMFile) as pam:
    line = pam.read().strip()
    pam = line.split(' ')[0]
    len_pam = int(line.split(' ')[1])
    guide_len = len(pam) - len_pam
    pos_beg = 0
    pos_end = None
    pam_begin = 0
    pam_end = len_pam * (-1)
    if len_pam < 0:
        guide_len = len(pam) + len_pam
        pam = pam[: (len_pam * (-1))]
        len_pam = len_pam * (-1)
        pos_beg = len_pam
        pos_end = None
        pam_begin = 0
        pam_end = len_pam
        pam_at_beginning = True
    else:
        pam = pam[(len_pam * (-1)):]
        pos_beg = 0
        pos_end = len_pam * (-1)
        pam_begin = len_pam * (-1)
        pam_end = None

# Get list of pam in pam directory
onlypam = [f for f in os.listdir(
    fileDir + 'pam') if isfile(join(fileDir + 'pam', f))]
removePamFile = False
for p in onlypam:
    if pam in p:
        removePamFile = True
        # Copy pam file into pam/tmpPam.txt
        with open(fileDir+"/pam/tmpPam.txt", "w") as pam_file:
            pam_file.write(pam+" "+str(len_pam))
        PAMFile = 'tmpPam.txt'
        break
else:  # The input pam is not available in the 'pam' folder
    with open(fileDir+"/pam/"+os.path.basename(PAMFile), "w") as pam_file:
        pam_file.write(pam+" "+str(len_pam))

if not os.path.exists(fileDir+"/annotations/"+cleanName+"_ref.Annotations.bed"):
    copy(annotationFile, fileDir+"/annotations/" +
         cleanName+"_ref.Annotations.bed")

if is_vcf_given and not os.path.exists(fileDir+"/samplesID/samples_"+cleanName+"_ref+"+cleanName+"_"+enrichedName+".txt"):
    copy(sampleFile, fileDir+"/samplesID/samples_" +
         cleanName+"_ref+"+cleanName+"_"+enrichedName+".txt")

if is_vcf_given:
    with open(fileDir+"Genomes/aggiunta_nuovo_genoma.txt", 'w') as log:
        log.write("1 Adding-variants")

    print("#####################################")
    print("Creating enriched genome")
    print("#####################################")
    os.system("crispritz.py add-variants "+VCFDir+"/ "+genomeDir+"/")
    if not os.path.exists("Genomes/"+cleanName+"_ref+"+cleanName+"_"+enrichedName):
        os.mkdir(fileDir+"/Genomes/"+cleanName +
                 "_ref+"+cleanName+"_"+enrichedName)
    for item in os.listdir(fileDir+"/variants_genome/SNPs_genome/"+cleanName+"_enriched/"):
        copy(fileDir+"/variants_genome/SNPs_genome/"+cleanName+"_enriched/"+item, fileDir+"/Genomes/" +
             cleanName+"_ref+"+cleanName+"_"+enrichedName+"/"+item)
    rmtree("variants_genome/")

with open(fileDir+"Genomes/aggiunta_nuovo_genoma.txt", 'w') as log:
    log.write("2 Indexing")
print("#####################################")
print("Creating indexes for reference genome")
print("#####################################")
os.system("crispritz.py index-genome " +
          cleanName+"_ref " +
          fileDir+"/Genomes/"+cleanName+"_ref/ " +
          fileDir+"/pam/"+os.path.basename(PAMFile)+" -bMax "+str(bMax))

if is_vcf_given:
    print("#####################################")
    print("Creating indexes for enriched genome")
    print("#####################################")
    os.system("crispritz.py index-genome "+cleanName+"_ref+"+cleanName+"_"+enrichedName +
              " "+fileDir+"/Genomes/"+cleanName+"_ref+"+cleanName+"_"+enrichedName +
              "/ "+fileDir+"/pam/"+os.path.basename(PAMFile)+" -bMax "+str(bMax))

    with open(fileDir+"Genomes/aggiunta_nuovo_genoma.txt", 'w') as log:
        log.write("3 Creating-dictionaries")
    print("#####################################")
    print("Creating dictionaries")
    print("#####################################")
    os.chdir(appDir+"/PostProcess/")
    if not os.path.exists(fileDir+"/dictionaries/"+cleanName+"_ref+"+cleanName+"_"+enrichedName):
        os.mkdir(fileDir+"/dictionaries/"+cleanName +
                 "_ref+"+cleanName+"_"+enrichedName)
    for item in os.listdir(VCFDir+"/"):
        print("Creating dictionary for file "+item)
        os.system("python creazione_dizionari.py "+VCFDir+"/"+item+" " +
                  fileDir+"/dictionaries/"+cleanName+"_ref+"+cleanName +
                  "_"+enrichedName+"/"+item)

if removePamFile:
    os.remove(fileDir + '/pam/' + PAMFile)

with open(fileDir+"Genomes/aggiunta_nuovo_genoma.txt", 'w') as log:
    log.write("4 Done")

"""
with open(fileDir+"/logGenomes/"+os.path.basename(os.path.normpath(genomeDir))+"+"+
               enrichedName+".log","w") as logFile:
    logFile.writelines("Reference Genome\t"+genomeDir+"\n")
    logFile.writelines("VCF\t"+VCFDir+"\n")
    logFile.writelines("PAM\t"+PAMFile+"\n")
    logFile.writelines("Annotation\t"+annotationFile+"\n")
    logFile.writelines("Sample List\t"+sampleFile+"\n")
    logFile.writelines("# Bulges\t"+str(bMax)+"\n")
    logFile.writelines("Name Enriched\t"+enrichedName+"\n")
"""

print("#####################################")
print("Procedure Finished")
print("#####################################")
