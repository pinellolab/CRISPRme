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

# matplotlib.use("TkAgg")
matplotlib.use('Agg')

warnings.filterwarnings("ignore")

plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
random.seed(a=None, version=2)
# final file containing all the results after post processing
inGuideFile = open(sys.argv[1], 'r')  # guide file used during search
inFinalFile = open(sys.argv[2], 'r')  # final result file from search
inSamplesIDFile = open(sys.argv[3], 'r').readlines()  # sampleID file
inSamplesIDFile.pop(0)  # pop header from sampleID file
# annotation file used during search
inAnnotationsFile = open(sys.argv[4], 'r')
outDir = sys.argv[5]  # directory to output the figures
threads = int(sys.argv[6])  # number of concurrent execution of image creation

web_server = False
population = False
# file_extension = 'pdf'
file_extension = 'png'

if '-ws' in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = 'png'


# def generatePlot(guide, guideDict, motifDict, mismatch, bulge, source):

#     # check if no targets are found for that combination source/totalcount and skip the execution
#     if guideDict['General'] == 0:
#         return
#     percentage_list = []
#     for elem in guideDict:
#         if float(guideDict['General']) != 0:
#             percentage_list.append(
#                 float(str(float(guideDict[elem])/float(guideDict['General']))[0:5]))
#         else:
#             percentage_list.append(float(0))

#     guideDataFrame = pd.DataFrame.from_dict(guideDict, orient='index')
#     guideDataFrame['Percentage'] = percentage_list
#     guideDataFrame.columns = ['Total', 'Percentage']
#     guideDataFrame = guideDataFrame.T
#     # number of variable
#     categories = list(guideDataFrame)[0:]
#     N = len(categories)

#     # We are going to plot the first line of the data frame.
#     # But we need to repeat the first value to close the circular graph:
#     values = guideDataFrame.loc['Percentage'].values.flatten().tolist()
#     values += values[:1]

#     # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
#     angles = [n / float(N) * 2 * pi for n in range(N)]
#     angles += angles[:1]

#     # Initialise the spider plot
#     ax = plt.subplot(2, 2, 1, polar=True)

#     for label, rot in zip(ax.get_xticklabels(), angles):
#         if (rot == 0):
#             label.set_horizontalalignment("center")
#         if (rot > 0):
#             label.set_horizontalalignment("left")
#         if (rot > 3):
#             label.set_horizontalalignment("center")
#         if (rot > 4):
#             label.set_horizontalalignment("right")

#     # Draw one axe per variable + add labels labels yet
#     plt.xticks(angles[:-1], categories, color='black', size=14)

#     # Draw ylabels
#     # # # Draw ylabels
#     ax.set_rlabel_position(0)
#     plt.yticks([0, 0.25, 0.50, 0.75], ["0", "0.25",
#                                        "0.50", "0.75"], color="black", size=12)
#     plt.ylim(0, 1)

#     # Fill area
#     ax.fill(angles, values, 'b', alpha=0.1)

#     # # # offset posizione y-axis
#     ax.set_theta_offset(pi / 2)
#     ax.set_theta_direction(-1)
#     # Plot data
#     ax.plot(angles, values, linewidth=1, linestyle='solid')

#     plt.subplot(2, 2, 2)
#     transpose_list = []
#     guideDataFrame = guideDataFrame.T
#     for elem in categories:
#         transpose_list.append(list(guideDataFrame.loc[elem]))
#         # transpose_list.append(
#         #     [guideDataFrame.loc[elem, 'Total'], guideDataFrame.loc[elem, 'Percentage']])

#     plt.axis('off')
#     table = plt.table(cellText=transpose_list, rowLabels=categories, colLabels=['Total', 'Percentage'],
#                       loc='best', colWidths=[0.25, 0.25])
#     table.auto_set_font_size(False)
#     table.set_fontsize(13)

#     totalMotif = [0]*len(guide)
#     for count in range(len(guide)):
#         for nuc in motifDict:
#             totalMotif[count] += motifDict[nuc][count]

#     maxmax = max(totalMotif)
#     for count in range(len(guide)):
#         for nuc in motifDict:
#             if maxmax != 0:
#                 motifDict[nuc][count] = float(
#                     motifDict[nuc][count]/float(maxmax))
#                 # motifDict[nuc][count] = float(
#                 #     str(float(motifDict[nuc][count])/float(maxmax))[0:5])

#     # ind = np.arange(0, len(guide), 1) + 0.15
#     ind = np.arange(0, len(guide), 1)
#     width = 0.7  # the width of the bars: can also be len(x) sequence

#     motif = plt.subplot(2, 1, 2, frameon=False)

#     A = np.array(motifDict['A'], dtype=float)
#     C = np.array(motifDict['C'], dtype=float)
#     G = np.array(motifDict['G'], dtype=float)
#     T = np.array(motifDict['T'], dtype=float)
#     RNA = np.array(motifDict['RNA'], dtype=float)
#     DNA = np.array(motifDict['DNA'], dtype=float)

#     p1 = plt.bar(ind, A, width, align='center')
#     p2 = plt.bar(ind, C, width, bottom=A, align='center')
#     p3 = plt.bar(ind, G, width, bottom=A+C, align='center')
#     p4 = plt.bar(ind, T, width, bottom=C+G+A, align='center')
#     p5 = plt.bar(ind, RNA, width, bottom=C+G+A+T, align='center')
#     p6 = plt.bar(ind, DNA, width, bottom=C+G+A+T+RNA,
#                  align='center')

#     # plt.xlim(0, len(guide))
#     # strArray = np.array([list(guide)])
#     plt.xticks(ticks=ind, labels=list(guide))

#     plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
#                ('A', 'C', 'G', 'T', 'bRNA', 'bDNA'), fontsize=13, loc='upper right', ncol=6)

#     # strArray = np.array([list(guide)])
#     # table = plt.table(cellText=strArray, loc='bottom',
#     #                   cellLoc='center', rowLoc='bottom')
#     # table.auto_set_font_size(False)
#     # table.set_fontsize(12)
#     # # table.scale(1, 1.6)
#     # table.xticks = ([])
#     # table.yticks = ([])
#     # plt.xticks

#     plt.suptitle(str(mismatch)+" Mismatches+ "+str(bulge)+" Bulge "+str(source),
#                  horizontalalignment='center', color='black', size=25)

#     plt.tight_layout()
#     plt.subplots_adjust(top=0.85, bottom=0.05, left=0.05,
#                         right=0.95, wspace=0.1)
#     plt.savefig(outDir+"/summary_single_guide_" + str(guide) + "_" + str(mismatch) +
#                 "."+str(bulge) + '_' + str(source) + "." + file_extension, format=file_extension)

#     plt.close('all')


def motifDictCreation(guide):
    # create time stamp for motif dict
    motifDict = dict()
    for mismatch in range(0, 11):
        motifDict[mismatch] = dict()
        for bulge in range(0, 6):
            motifDict[mismatch][bulge] = dict()
            motifDict[mismatch][bulge]['TOTAL'] = dict()
            motifDict[mismatch][bulge]['TOTAL']['A'] = [0]*len(guide)
            motifDict[mismatch][bulge]['TOTAL']['C'] = [0]*len(guide)
            motifDict[mismatch][bulge]['TOTAL']['G'] = [0]*len(guide)
            motifDict[mismatch][bulge]['TOTAL']['T'] = [0]*len(guide)
            motifDict[mismatch][bulge]['TOTAL']['RNA'] = [0]*len(guide)
            motifDict[mismatch][bulge]['TOTAL']['DNA'] = [0]*len(guide)
            for sample in populationDict:
                superpop = populationDict[sample][0]
                pop = populationDict[sample][1]
                motifDict[mismatch][bulge][sample] = dict()
                motifDict[mismatch][bulge][sample]['A'] = [0]*len(guide)
                motifDict[mismatch][bulge][sample]['C'] = [0]*len(guide)
                motifDict[mismatch][bulge][sample]['G'] = [0]*len(guide)
                motifDict[mismatch][bulge][sample]['T'] = [0]*len(guide)
                motifDict[mismatch][bulge][sample]['RNA'] = [0]*len(guide)
                motifDict[mismatch][bulge][sample]['DNA'] = [0]*len(guide)
                if superpop not in motifDict[mismatch][bulge]:
                    motifDict[mismatch][bulge][superpop] = dict()
                    motifDict[mismatch][bulge][superpop]['A'] = [0]*len(guide)
                    motifDict[mismatch][bulge][superpop]['C'] = [0]*len(guide)
                    motifDict[mismatch][bulge][superpop]['G'] = [0]*len(guide)
                    motifDict[mismatch][bulge][superpop]['T'] = [0]*len(guide)
                    motifDict[mismatch][bulge][superpop]['RNA'] = [
                        0]*len(guide)
                    motifDict[mismatch][bulge][superpop]['DNA'] = [
                        0]*len(guide)
                if pop not in motifDict[mismatch][bulge]:
                    motifDict[mismatch][bulge][pop] = dict()
                    motifDict[mismatch][bulge][pop]['A'] = [0]*len(guide)
                    motifDict[mismatch][bulge][pop]['C'] = [0]*len(guide)
                    motifDict[mismatch][bulge][pop]['G'] = [0]*len(guide)
                    motifDict[mismatch][bulge][pop]['T'] = [0]*len(guide)
                    motifDict[mismatch][bulge][pop]['RNA'] = [0]*len(guide)
                    motifDict[mismatch][bulge][pop]['DNA'] = [0]*len(guide)
    return motifDict


def guideDictCreation():
    # create one stamp for guideDict
    guideDict = dict()
    for mismatch in range(0, 11):
        guideDict[mismatch] = dict()
        for bulge in range(0, 6):
            guideDict[mismatch][bulge] = dict()
            guideDict[mismatch][bulge]['TOTAL'] = dict()
            guideDict[mismatch][bulge]['TOTAL']['General'] = 0
            for annot in annotationsSet:
                guideDict[mismatch][bulge]['TOTAL'][annot] = 0
            for sample in populationDict:
                superpop = populationDict[sample][0]
                pop = populationDict[sample][1]
                guideDict[mismatch][bulge][sample] = dict()
                guideDict[mismatch][bulge][sample]['General'] = 0
                if superpop not in guideDict[mismatch][bulge]:
                    guideDict[mismatch][bulge][superpop] = dict()
                    guideDict[mismatch][bulge][superpop]['General'] = 0
                if pop not in guideDict[mismatch][bulge]:
                    guideDict[mismatch][bulge][pop] = dict()
                    guideDict[mismatch][bulge][pop]['General'] = 0
                for annot in annotationsSet:
                    guideDict[mismatch][bulge][superpop][annot] = 0
                    guideDict[mismatch][bulge][pop][annot] = 0
                    guideDict[mismatch][bulge][sample][annot] = 0
    return guideDict


def fillDict(guide, guideDict, motifDict):
    # fill dictionary with info read from the final file
    inFinalFile.seek(0)
    if '#' in inFinalFile.readline():
        print('SKIP HEADER')
    else:
        inFinalFile.seek(0)

    for line in inFinalFile:
        # print(line)
        split = line.strip().split('\t')
        if guide not in split[15]:
            continue
        alignedSequence = split[2]  # aligned target with the guide
        mismatch = int(split[8])  # mm extracted from target file
        bulge = int(split[9])  # bul extracted from target file

        # add target to TOTAL dict
        # print(mismatch, bulge, guideDict[mismatch][bulge]['TOTAL']['General'])
        guideDict[mismatch][bulge]['TOTAL']['General'] += 1
        # find annotations and add data to the TOTAL dict
        if split[14] != 'n':
            annotationsList = split[14].strip().split(',')
            for annotation in annotationsList:
                guideDict[mismatch][bulge]['TOTAL'][annotation] += 1
        # find motif in X and RNA/DNA targets
        if 'DNA' not in split[0]:
            for count, nucleotide in enumerate(alignedSequence):
                if nucleotide.islower():
                    motifDict[mismatch][bulge]['TOTAL'][nucleotide.upper()][count] += 1
                elif nucleotide == '-':
                    motifDict[mismatch][bulge]['TOTAL'][split[0]][count] += 1
                if guide[count] == 'N':
                    motifDict[mismatch][bulge]['TOTAL'][nucleotide.upper()][count] += 1
        else:
            alignedGuide = split[1]
            for count, nucleotide in enumerate(alignedGuide[bulge:]):
                if nucleotide == '-':
                    motifDict[mismatch][bulge]['TOTAL'][split[0]][count] += 1
            for count, nucleotide in enumerate(alignedSequence[bulge:]):
                if nucleotide.islower():
                    motifDict[mismatch][bulge]['TOTAL'][nucleotide.upper()][count] += 1
                if guide[0] != 'N':
                    if guide[count] == 'N':
                        motifDict[mismatch][bulge]['TOTAL'][nucleotide.upper()][count] += 1
            for count, ennes in enumerate(guide): #to correct reading of guides N's when upstream PAM
                if ennes == 'N':
                    motifDict[mismatch][bulge]['TOTAL'][alignedSequence[count].upper()][count] += 1
                else:
                    break
        # check if samples are available and take them if yes
        samplePresent = False
        if 'n' not in split[13] and 'NO_SAMPLES' not in split[13]:
            samplePresent = split[13].strip().split(',')
        if samplePresent != False:  # if samples are present analyze them
            alreadySampled = set()
            for sample in samplePresent:
                superpop = populationDict[sample][0]
                pop = populationDict[sample][1]
                if superpop not in alreadySampled:  # check if superpop is already been updated for the current target to avoid duplicate count
                    alreadySampled.add(superpop)
                    guideDict[mismatch][bulge][superpop]['General'] += 1
                    # find motif in X and RNA/DNA targets
                    if 'DNA' not in split[0]:
                        for count, nucleotide in enumerate(alignedSequence):
                            if nucleotide.islower():
                                motifDict[mismatch][bulge][superpop][nucleotide.upper()][count] += 1
                            elif nucleotide == '-':
                                motifDict[mismatch][bulge][superpop][split[0]][count] += 1
                            if guide[count] == 'N':
                                motifDict[mismatch][bulge][superpop][nucleotide.upper()][count] += 1
                    else:
                        alignedGuide = split[1]
                        for count, nucleotide in enumerate(alignedGuide[bulge:]):
                            if nucleotide == '-':
                                motifDict[mismatch][bulge][superpop][split[0]][count] += 1
                        for count, nucleotide in enumerate(alignedSequence[bulge:]):
                            if nucleotide.islower():
                                motifDict[mismatch][bulge][superpop][nucleotide.upper()][count] += 1
                            if guide[0] != 'N':
                                if guide[count] == 'N':
                                    motifDict[mismatch][bulge][superpop][nucleotide.upper()][count] += 1
                        for count, ennes in enumerate(guide): #to correct reading of guides N's when upstream PAM
                            if ennes == 'N':
                                motifDict[mismatch][bulge][superpop][alignedSequence[count].upper()][count] += 1
                            else:
                                break
                    if split[14] != 'n':
                        annotationsList = split[14].strip().split(',')
                        for annotation in annotationsList:
                            guideDict[mismatch][bulge][superpop][annotation] += 1
                if pop not in alreadySampled:  # check if pop is already been updated for the current target to avoid duplicate count
                    alreadySampled.add(pop)
                    guideDict[mismatch][bulge][pop]['General'] += 1
                    # find motif in X and RNA/DNA targets
                    if 'DNA' not in split[0]:
                        for count, nucleotide in enumerate(alignedSequence):
                            if nucleotide.islower():
                                motifDict[mismatch][bulge][pop][nucleotide.upper()][count] += 1
                            elif nucleotide == '-':
                                motifDict[mismatch][bulge][pop][split[0]][count] += 1
                            if guide[count] == 'N':
                                motifDict[mismatch][bulge][pop][nucleotide.upper()][count] += 1
                    else:
                        alignedGuide = split[1]
                        for count, nucleotide in enumerate(alignedGuide[bulge:]):
                            if nucleotide == '-':
                                motifDict[mismatch][bulge][pop][split[0]][count] += 1
                        for count, nucleotide in enumerate(alignedSequence[bulge:]):
                            if nucleotide.islower():
                                motifDict[mismatch][bulge][pop][nucleotide.upper()][count] += 1
                            if guide[0] != 'N':
                                if guide[count] == 'N':
                                    motifDict[mismatch][bulge][pop][nucleotide.upper()][count] += 1
                        for count, ennes in enumerate(guide): #to correct reading of guides N's when upstream PAM
                            if ennes == 'N':
                                motifDict[mismatch][bulge][pop][alignedSequence[count].upper()][count] += 1
                            else:
                                break
                    if split[14] != 'n':
                        annotationsList = split[14].strip().split(',')
                        for annotation in annotationsList:
                            guideDict[mismatch][bulge][pop][annotation] += 1
                # add target to single sample, in this case duplicate are allowed
                guideDict[mismatch][bulge][sample]['General'] += 1
                # find motif in X and RNA/DNA targets
                if 'DNA' not in split[0]:
                    for count, nucleotide in enumerate(alignedSequence):
                        if nucleotide.islower():
                            motifDict[mismatch][bulge][sample][nucleotide.upper()][count] += 1
                        elif nucleotide == '-':
                            motifDict[mismatch][bulge][sample][split[0]][count] += 1
                        if guide[count] == 'N':
                            motifDict[mismatch][bulge][sample][nucleotide.upper()][count] += 1
                else:
                    alignedGuide = split[1]
                    for count, nucleotide in enumerate(alignedGuide[bulge:]):
                        if nucleotide == '-':
                            motifDict[mismatch][bulge][sample][split[0]][count] += 1
                    for count, nucleotide in enumerate(alignedSequence[bulge:]):
                        if nucleotide.islower():
                            motifDict[mismatch][bulge][sample][nucleotide.upper()][count] += 1
                        if guide[0] != 'N':
                            if guide[count] == 'N':
                                motifDict[mismatch][bulge][pop][nucleotide.upper()][count] += 1
                    for count, ennes in enumerate(guide): #to correct reading of guides N's when upstream PAM
                        if ennes == 'N':
                            motifDict[mismatch][bulge][pop][alignedSequence[count].upper()][count] += 1
                        else:
                            break
                if split[14] != 'n':
                    annotationsList = split[14].strip().split(',')
                    for annotation in annotationsList:
                        guideDict[mismatch][bulge][sample][annotation] += 1


# data containing populations and annotations
annotationsSet = set()
populationDict = dict()

# read all the annotations
for line in inAnnotationsFile:
    for elem in line.strip().split('\t')[3].strip().split(','):
        annotationsSet.add(elem)
annotationsSet = sorted(annotationsSet)

# read all the population, superpop and samples
for line in inSamplesIDFile:
    split = line.strip().split('\t')
    populationDict[split[0]] = [split[2], split[1]]

for guide in inGuideFile:
    # guide = guide.strip().replace('N', '')
    guide = guide.strip()
    guideDict = guideDictCreation()
    motifDict = motifDictCreation(guide)
    fillDict(guide, guideDict, motifDict)
    json.dump(guideDict, open(outDir+f"/guide_dict_{guide}.json", 'w'))
    json.dump(motifDict, open(outDir+f"/motif_dict_{guide}.json", 'w'))
    # print('Starting image creation for guide: ', guide)
    # if __name__ == '__main__':
    #     pool = multiprocessing.Pool(threads)
    #     for mismatch in range(0, 8):
    #         for bulge in range(0, 3):
    #             # skip creation if no target in global category
    #             if guideDict[mismatch][bulge]['TOTAL']['General'] == 0:
    #                 continue
    #             for elem in guideDict[mismatch][bulge]:
    #                 pool.apply_async(generatePlot, args=(
    #                     guide, guideDict[mismatch][bulge][elem], motifDict[mismatch][bulge][elem], mismatch, bulge, elem))
    #     pool.close()
    #     pool.join()
