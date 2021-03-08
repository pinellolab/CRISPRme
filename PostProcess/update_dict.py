#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 23:27:45 2020

@author: whitebreeze
"""
import os
import sys
import gzip
from change_dict import updateDictionary
from shutil import copy, rmtree, move


def getAssociatedPAM(genomeName, fileDir):
    pams, bMaxs = [], []
    for d in os.listdir(fileDir+"/genome_library/"):
        if genomeName in d:
            parts = d.split("_")
            pams.append(parts[0])
            bMaxs.append(parts[1])
    return pams, bMaxs


# TODO to prevent pam with no guide len in the name, put 25 Ns , check if positions are correct
def createTempPAM(fileDir, pam, pamFile):
    parts = pamFile.split("-")
    if parts[0] != pam:
        with open(fileDir+"/pam/tempPAM.txt", "w") as tempPAM:
            fullPAM = "".join(["N"]*int(parts[0][0:2])) + \
                pam + " " + str(len(pam))
            tempPAM.write(fullPAM)
    else:
        with open(fileDir+"/pam/tempPAM.txt", "w") as tempPAM:
            fullPAM = pam + "".join(["N"]*int(parts[1][0:2])
                                    ) + " -" + str(len(pam))
            tempPAM.write(fullPAM)


fileDir = sys.argv[1]
oldDicts = sys.argv[2]
VCFDir = sys.argv[3]
sampleFile = sys.argv[4]

os.chdir(fileDir)
with open("dictionaries/update_dizionari.txt", 'w') as log:
    log.write("0 Updating-Dictionaries")

print("#####################################")
print("Updating dictionaries")
print("#####################################")
dict_vcf = []
for i, VCFFile in enumerate(os.listdir(VCFDir)):
    with gzip.open(VCFDir+"/"+VCFFile, 'rb') as targets:
        # Skip vcf header
        for line in targets:
            line = line.decode('ascii')
            if ('#CHROM') in line:
                # Save this header for retrieving sample id
                column_vcf = line.strip().split('\t')
                break

        # Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        line = targets.readline()
        line = line.decode('ascii').strip().split('\t')
        chrom = line[0]
        # print(chrom)
        if "chr" not in chrom:
            chrom = "chr"+chrom
        for oldDict in os.listdir(oldDicts):
            if oldDict[8:len(oldDict)-5] == chrom:
                dict_vcf.append((oldDict, VCFFile))
                break
    if len(dict_vcf) != i+1:
        print("WARNING: No dictionaries associated to "+VCFFile)


for oldDict, VCFFile in dict_vcf:
    print("Updating "+oldDict)
    updateDictionary(oldDicts+"/"+oldDict, VCFDir+"/"+VCFFile)
    # os.system("python change_dict.py "+oldDicts+"/"+oldDict+" "+VCFDir+"/"+VCFFile)


with open("dictionaries/update_dizionari.txt", 'w') as log:
    log.write("1 Updating-Genome")
print("#####################################")
print("Updating enriched genome with new variants")
print("#####################################")
genomeEnr = os.path.basename(oldDicts).split('dictionary_')[-1]
print(genomeEnr)
for item in os.listdir(fileDir+"/Genomes/"+genomeEnr+"/"):
    chrom = item.split('.')[0]
    print(chrom)
    os.rename(fileDir+"/Genomes/"+genomeEnr+"/"+item,
              fileDir+"/Genomes/"+genomeEnr+"/"+chrom+".fa")
print(os.listdir(fileDir+"/Genomes/"+genomeEnr+"/"))
os.system("crispritz.py add-variants "+VCFDir +
          "/ "+fileDir+"/Genomes/"+genomeEnr+"/")
rmtree(fileDir+"/Genomes/"+genomeEnr+"/")
if not os.path.exists("Genomes/"+genomeEnr):
    os.mkdir(fileDir+"/Genomes/"+genomeEnr)
for item in os.listdir(fileDir+"/variants_genome/SNPs_genome/" + genomeEnr + '_enriched/'):
    copy(fileDir+"/variants_genome/SNPs_genome/" + genomeEnr +
         '_enriched/'+item, fileDir+"/Genomes/"+genomeEnr+"/"+item)
# move(fileDir+"/variants_genome/SNPs_genome/" + genomeEnr + '_enriched', fileDir+"/Genomes/"+genomeEnr) # problema sulla mia VM con Ubuntu
rmtree("variants_genome/")

with open("dictionaries/update_dizionari.txt", 'w') as log:
    log.write("2 Updating-Index")

print("#####################################")
print("Indexing enriched genome with new variants")
print("#####################################")
pams, bMaxs = getAssociatedPAM(genomeEnr, fileDir)
for pam, bMax in zip(pams, bMaxs):
    fullPAM = ""
    for d in os.listdir(fileDir+"/pam/"):
        if pam in d:
            fullPAM = d
            break
    if fullPAM == "":
        print("PAM NOT FOUND!")
        raise("PAM not found")
    # print(fullPAM)
    createTempPAM(fileDir, pam, fullPAM)
    os.system("crispritz.py index-genome " +
              genomeEnr+" " +
              fileDir+"/Genomes/"+genomeEnr+"/ " +
              fileDir+"/pam/tempPAM.txt"+" -bMax "+str(bMax))
    os.remove(fileDir+"/pam/tempPAM.txt")

with open("dictionaries/update_dizionari.txt", 'w') as log:
    log.write("3 Updating-Samples")

print("#####################################")
print("Adding samplesID")
print("#####################################")
with open(fileDir+'/samplesID/samples_'+genomeEnr+".txt", 'a') as oldSampleFile:
    with open(sampleFile, 'r') as newSampleFile:
        newSampleFile.readline()  # header
        for line in newSampleFile:
            oldSampleFile.write("\n"+line.strip())

with open("dictionaries/update_dizionari.txt", 'w') as log:
    log.write("4 Done")

print("#####################################")
print("Procedure Finished")
print("#####################################")
