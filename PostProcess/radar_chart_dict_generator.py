#!/usr/bin/env python

# Libraries
from operator import itemgetter
from os.path import isfile, join
from os import listdir
import os
import warnings
import glob
from itertools import islice
import sys
import numpy as np
import scipy.spatial.distance as sp
from math import pi
import pandas as pd
from matplotlib import patches as mpatches
from matplotlib import pyplot as plt
import math
import matplotlib
import time
import random
import multiprocessing
import json
import pysam

warnings.filterwarnings("ignore")
# matplotlib.use("TkAgg")
matplotlib.use("Agg")


plt.style.use("seaborn-poster")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

random.seed(a=None, version=2)
# final file containing all the results after post processing
inGuideFile = open(sys.argv[1], "r")  # guide file used during search
inFinalFile = open(sys.argv[2], "r")  # final result file from search
inSamplesIDFile = open(sys.argv[3], "r").readlines()  # sampleID file
inSamplesIDFile.pop(0)  # pop header from sampleID file
# annotation file used during search
annotation_fname = sys.argv[4]
inAnnotationsFile = None if annotation_fname == "vuoto.txt" else pysam.TabixFile(sys.argv[4], "r")
outDir = sys.argv[5]  # directory to output the figures
threads = int(sys.argv[6])  # number of concurrent execution of image creation
max_mm = int(sys.argv[7])
max_bulges = int(sys.argv[8])
# criteria to generate the plots (CFD,CRISTA,FEWEST)
selection_criteria = sys.argv[9]

web_server = False
population = False
# file_extension = 'pdf'
file_extension = "png"

if "-ws" in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = "png"

def motifDictCreation(guide):
    # create time stamp for motif dict
    motifDict = dict()
    for total in range(0, 15):
        motifDict[total] = dict()
        motifDict[total]["A"] = [0] * len(guide)
        motifDict[total]["C"] = [0] * len(guide)
        motifDict[total]["G"] = [0] * len(guide)
        motifDict[total]["T"] = [0] * len(guide)
        motifDict[total]["RNA"] = [0] * len(guide)
        motifDict[total]["DNA"] = [0] * len(guide)
    return motifDict


def guideDictCreation(annotationsSet):
    # create one stamp for guideDict
    guideDict = dict()
    for total in range(0, 15):
        guideDict[total] = dict()
        guideDict[total]["General"] = 0
        for annot in annotationsSet:
            guideDict[total][annot] = 0
    return guideDict


def fillDict(guide, guideDict, motifDict):
    # fill dictionary with info read from the final file
    inFinalFile.seek(0)
    if "#" in inFinalFile.readline():
        print("SKIP HEADER")
    else:
        inFinalFile.seek(0)

    for line in inFinalFile:
        split = line.strip().split("\t")
        if guide not in split[15]:
            continue
        alignedSequence = split[2]  # aligned target with the guide
        mismatch = int(split[8])  # mm extracted from target file
        bulge = int(split[9])  # bul extracted from target file
        total = mismatch + bulge

        # add target to TOTAL dict
        for over in range(total, 15):
            # count the target in each level, from total to infinite, this way we autoinclude
            # each target in the proper level and in all the following levels
            guideDict[over]["General"] += 1
            if split[14] != "NA":
                annotationsList = split[14].strip().split(",")
                # to avoid duplicate categories
                annotationsList = set(annotationsList)
                for annotation in annotationsList:
                    if "CTCF-bound" in annotation:
                        guideDict[over]["CTCF-only"] += 1
                        guideDict[over][annotation.split(";")[0]] += 1
                    elif "_gencode" in annotation:
                        guideDict[over][annotation.replace("_gencode", "")] += 1
            if "RNA,DNA" in split[0]:
                continue
            # find motif in X and RNA/DNA targets
            if "DNA" not in split[0]:
                for count, nucleotide in enumerate(alignedSequence):
                    if nucleotide.islower():
                        motifDict[over][nucleotide.upper()][count] += 1
                    elif nucleotide == "-":
                        motifDict[over][split[0]][count] += 1
                    if guide[count] == "N":
                        motifDict[over][nucleotide.upper()][count] += 1
            else:
                alignedGuide = split[1]
                for count, nucleotide in enumerate(alignedGuide[bulge:]):
                    if nucleotide == "-":
                        motifDict[over][split[0]][count] += 1
                for count, nucleotide in enumerate(alignedSequence[bulge:]):
                    if nucleotide.islower():
                        motifDict[over][nucleotide.upper()][count] += 1
                    if guide[0] != "N":
                        if guide[count] == "N":
                            motifDict[over][nucleotide.upper()][count] += 1
                # to correct reading of guides N's when upstream PAM
                for count, ennes in enumerate(guide):
                    if ennes == "N":
                        motifDict[over][alignedSequence[count].upper()][count] += 1
                    else:
                        break


# data containing populations and annotations
annotationsSet = set()
populationDict = dict()

# read all the annotations
if inAnnotationsFile is not None:
    for line in inAnnotationsFile.fetch():
        if (
            "vuoto.txt" in sys.argv[4]
        ):  # se vuoto.txt usato come annotazione, skippa lettura
            break
        annotations_list = line.strip().split("\t")[3].split(",")
        for annotation in annotations_list:
            if "_personal" not in annotation:
                annotationsSet.add(annotation.replace("_gencode", ""))

annotationsSet.add("CTCF-only")
annotationsSet = sorted(annotationsSet)

# read all the population, superpop and samples
for line in inSamplesIDFile:
    split = line.strip().split("\t")
    populationDict[split[0]] = [split[2], split[1]]

for guide in inGuideFile:
    guide = guide.strip()
    guideDict = guideDictCreation(annotationsSet)
    motifDict = motifDictCreation(guide)
    fillDict(guide, guideDict, motifDict)
    # for total in range(max_mm+max_bulges):
    #     generatePlot(guide, guideDict[total],
    #                  motifDict[total], total, 0, 'TOTAL')
    json.dump(
        guideDict, open(outDir + f"/.guide_dict_{guide}_{selection_criteria}.json", "w")
    )
    json.dump(
        motifDict, open(outDir + f"/.motif_dict_{guide}_{selection_criteria}.json", "w")
    )
