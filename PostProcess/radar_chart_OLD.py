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

# matplotlib.use("TkAgg")
matplotlib.use('Agg')

warnings.filterwarnings("ignore")

# argv 1 is Guide
# argv 2 is total value
# argv 3 is annotation file (new version)
# argv 4 is extended profile (can be 'no' for web server)
# argv 5 is second summary for barplot (optional, 'no' if not given in input)
# argv 6 is gecko Summary file
# argv 7 is for web server: create png instead of pguideDataFrame (-ws) (optional)
# argv 8 is for sample/pop/superpop report (-sample HG001/EUR/TSI) with this option, the annotation file is .sample_annotation.GUIDE.sample.txt
# TODO if argv8 is selected, from argv6 get the directory containing the summary of the specific sample for gecko -> line 618
# NOTE al momento lo script funziona con gecko.annotation.summary; la comparison con i sample la fa su annotation.summary. Quando l'analisi per sample sarà fatta, la comparison
# con i sampla dovrà essere fatta sul file specifico
plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams["figure.figsize"] = [100, 80]
# SIZE_GECKO = 123411  # NOTE modify if new gecko annotations are done
# SIZE_GECKO = 111671
random.seed(a=None, version=2)
# final file containing all the results after post processing
inGuideFile = open(sys.argv[1], 'r')
inFinalFile = open(sys.argv[2], 'r')
inSamplesIDFile = open(sys.argv[3], 'r').readlines()
inSamplesIDFile.pop(0)
inAnnotationsFile = open(sys.argv[4], 'r')
outDir = sys.argv[5]
threads = int(sys.argv[6])

web_server = False
population = False
# file_extension = 'pdf'
file_extension = 'png'

if '-ws' in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = 'png'


def generatePlot(guide, guideDict, motifDict, totalcount, source):

    # check if no targets are found for that combination source/totalcount and skip the execution
    if guideDict['General'] == 0:
        return
    percentage_list = []
    for elem in guideDict:
        if float(guideDict['General']) != 0:
            percentage_list.append(
                float(str(float(guideDict[elem])/float(guideDict['General']))[0:5]))
        else:
            percentage_list.append(float(0))

    guideDataFrame = pd.DataFrame.from_dict(guideDict, orient='index')
    guideDataFrame['Percentage'] = percentage_list
    guideDataFrame.columns = ['Total', 'Percentage']
    guideDataFrame = guideDataFrame.T
    # number of variable
    categories = list(guideDataFrame)[0:]
    N = len(categories)

    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values = guideDataFrame.loc['Percentage'].values.flatten().tolist()
    values += values[:1]

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(2, 2, 1, polar=True)

    for label, rot in zip(ax.get_xticklabels(), angles):
        if (rot == 0):
            label.set_horizontalalignment("center")
        if (rot > 0):
            label.set_horizontalalignment("left")
        if (rot > 3):
            label.set_horizontalalignment("center")
        if (rot > 4):
            label.set_horizontalalignment("right")

    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='grey', size=10)

    # Draw ylabels
    # # # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0, 0.25, 0.50, 0.75], ["0", "0.25",
                                       "0.50", "0.75"], color="black", size=12)
    plt.ylim(0, 1)

    # Fill area
    ax.fill(angles, values, 'b', alpha=0.1)

    # # # offset posizione y-axis
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid')

    plt.subplot(2, 2, 2)
    transpose_list = []
    for count, elem in enumerate(categories):
        transpose_list.append(
            [guideDataFrame.T.iloc[count]['Total'], guideDataFrame.T.iloc[count]['Percentage']])

    plt.axis('off')
    table = plt.table(cellText=transpose_list, rowLabels=categories, colLabels=['Total', 'Percentage'],
                      loc='center', colWidths=[0.35 for x in guideDataFrame.columns])

    totalMotif = [0]*len(guide)
    for count in range(len(guide)):
        for nuc in motifDict:
            totalMotif[count] += motifDict[nuc][count]

    maxmax = max(totalMotif)
    for count in range(len(guide)):
        for nuc in motifDict:
            if maxmax != 0:
                motifDict[nuc][count] = float(
                    str(float(motifDict[nuc][count])/float(maxmax))[0:5])

    ind = np.arange(0, len(guide), 1) + 0.15
    width = 0.7  # the width of the bars: can also be len(x) sequence

    motif = plt.subplot(2, 1, 2, frameon=False)

    A = np.array(motifDict['A'], dtype=float)
    C = np.array(motifDict['C'], dtype=float)
    G = np.array(motifDict['G'], dtype=float)
    T = np.array(motifDict['T'], dtype=float)
    RNA = np.array(motifDict['RNA'], dtype=float)

    p1 = plt.bar(ind, A, width, color='red', align='edge')
    p2 = plt.bar(ind, C, width, color='blue', bottom=A, align='edge')
    p3 = plt.bar(ind, G, width, color='green', bottom=A+C, align='edge')
    p4 = plt.bar(ind, T, width, color='purple', bottom=C+G+A, align='edge')
    p5 = plt.bar(ind, RNA, width, color='magenta',
                 bottom=C+G+A+T, align='edge')
    # p6 = plt.bar(ind, RNA, width, color='gold', bottom=C+G+A+T+DNA, align='edge')
    # p6 = plt.bar(ind, RNA, width, color='gold', bottom=C+G+A+T, align='edge')
    plt.xlim(0, len(guide))
    plt.xticks([])

    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]),
               ('A', 'C', 'G', 'T', 'RNA'), fontsize=18)

    strArray = np.array([list(guide)])
    table = plt.table(cellText=strArray, loc='bottom',
                      cellLoc='center', rowLoc='bottom')
    table.auto_set_font_size(False)
    table.set_fontsize(18)
    table.scale(1, 1.6)
    table.xticks = ([])
    table.yticks = ([])

    plt.suptitle(str(totalcount)+" Mismatches + Bulge",
                 horizontalalignment='center', color='black', size=25)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.07, left=0.09,
                        right=0.99, wspace=0.25)
    plt.savefig(outDir+"/summary_single_guide_" + str(guide) + "_"+str(totalcount) +
                "_total"+'_'+source + "." + file_extension, format=file_extension)

    plt.close('all')


def motifDictCreation(guide):
    # create time stamp for motif dict
    motifDict = dict()
    for totalcount in range(0, 10):
        motifDict[totalcount] = dict()
        motifDict[totalcount]['REF'] = dict()
        motifDict[totalcount]['REF']['A'] = [0]*len(guide)
        motifDict[totalcount]['REF']['C'] = [0]*len(guide)
        motifDict[totalcount]['REF']['G'] = [0]*len(guide)
        motifDict[totalcount]['REF']['T'] = [0]*len(guide)
        motifDict[totalcount]['REF']['RNA'] = [0]*len(guide)
        for sample in populationDict:
            superpop = populationDict[sample][0]
            pop = populationDict[sample][1]
            motifDict[totalcount][sample] = dict()
            motifDict[totalcount][sample]['A'] = [0]*len(guide)
            motifDict[totalcount][sample]['C'] = [0]*len(guide)
            motifDict[totalcount][sample]['G'] = [0]*len(guide)
            motifDict[totalcount][sample]['T'] = [0]*len(guide)
            motifDict[totalcount][sample]['RNA'] = [0]*len(guide)
            if superpop not in motifDict[totalcount]:
                motifDict[totalcount][superpop] = dict()
                motifDict[totalcount][superpop]['A'] = [0]*len(guide)
                motifDict[totalcount][superpop]['C'] = [0]*len(guide)
                motifDict[totalcount][superpop]['G'] = [0]*len(guide)
                motifDict[totalcount][superpop]['T'] = [0]*len(guide)
                motifDict[totalcount][superpop]['RNA'] = [0]*len(guide)
            if pop not in motifDict[totalcount]:
                motifDict[totalcount][pop] = dict()
                motifDict[totalcount][pop]['A'] = [0]*len(guide)
                motifDict[totalcount][pop]['C'] = [0]*len(guide)
                motifDict[totalcount][pop]['G'] = [0]*len(guide)
                motifDict[totalcount][pop]['T'] = [0]*len(guide)
                motifDict[totalcount][pop]['RNA'] = [0]*len(guide)
    return motifDict


def guideDictCreation():
    # create one stamp for guideDict
    guideDict = dict()
    for totalcount in range(0, 10):
        guideDict[totalcount] = dict()
        guideDict[totalcount]['REF'] = dict()
        guideDict[totalcount]['REF']['General'] = 0
        for annot in annotationsSet:
            guideDict[totalcount]['REF'][annot] = 0
        for sample in populationDict:
            superpop = populationDict[sample][0]
            pop = populationDict[sample][1]
            guideDict[totalcount][sample] = dict()
            guideDict[totalcount][sample]['General'] = 0
            if superpop not in guideDict[totalcount]:
                guideDict[totalcount][superpop] = dict()
                guideDict[totalcount][superpop]['General'] = 0
            if pop not in guideDict[totalcount]:
                guideDict[totalcount][pop] = dict()
                guideDict[totalcount][pop]['General'] = 0
            for annot in annotationsSet:
                # guideDict[totalcount]['REF'][annot] = 0
                guideDict[totalcount][superpop][annot] = 0
                guideDict[totalcount][pop][annot] = 0
                guideDict[totalcount][sample][annot] = 0
    return guideDict


def fillDict(guide, guideDict, motifDict):
    # fill dictionary with info read from the final file
    if '#' in inFinalFile.readline():
        print('SKIP HEADER')
    else:
        inFinalFile.seek(0)

    for line in inFinalFile:
        split = line.strip().split('\t')
        alignedSequence = split[2]  # aligned target with the guide
        totalcount = int(split[10])  # total= mm+bul

        # add target to reference dict
        guideDict[totalcount]['REF']['General'] += 1
        # find annotations and add data to the reference dict
        if split[14] != 'n':
            annotationsList = split[14].strip().split(',')
            for annotation in annotationsList:
                guideDict[totalcount]['REF'][annotation] += 1
        # find motif in X and RNA targets TODO: integrate DNA targets
        if 'DNA' not in split[0]:
            for count, nucleotide in enumerate(alignedSequence):
                if nucleotide.islower():
                    motifDict[totalcount]['REF'][nucleotide.upper()
                                                 ][count] += 1
                elif nucleotide == '-':
                    motifDict[totalcount]['REF'][split[0]][count] += 1

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
                    guideDict[totalcount][superpop]['General'] += 1
                    if 'DNA' not in split[0]:
                        for count, nucleotide in enumerate(alignedSequence):
                            if nucleotide.islower():
                                motifDict[totalcount][superpop][nucleotide.upper(
                                )][count] += 1
                            elif nucleotide == '-':
                                motifDict[totalcount][superpop][split[0]
                                                                ][count] += 1
                    if split[14] != 'n':
                        annotationsList = split[14].strip().split(',')
                        for annotation in annotationsList:
                            guideDict[totalcount][superpop][annotation] += 1
                if pop not in alreadySampled:  # check if pop is already been updated for the current target to avoid duplicate count
                    alreadySampled.add(pop)
                    guideDict[totalcount][pop]['General'] += 1
                    if 'DNA' not in split[0]:
                        for count, nucleotide in enumerate(alignedSequence):
                            if nucleotide.islower():
                                motifDict[totalcount][pop][nucleotide.upper()
                                                           ][count] += 1
                            elif nucleotide == '-':
                                motifDict[totalcount][pop][split[0]
                                                           ][count] += 1
                    if split[14] != 'n':
                        annotationsList = split[14].strip().split(',')
                        for annotation in annotationsList:
                            guideDict[totalcount][pop][annotation] += 1
                # add target to single sample, in this case duplicate are allowed
                guideDict[totalcount][sample]['General'] += 1
                if 'DNA' not in split[0]:
                    for count, nucleotide in enumerate(alignedSequence):
                        if nucleotide.islower():
                            motifDict[totalcount][sample][nucleotide.upper()
                                                          ][count] += 1
                        elif nucleotide == '-':
                            motifDict[totalcount][sample][split[0]][count] += 1
                if split[14] != 'n':
                    annotationsList = split[14].strip().split(',')
                    for annotation in annotationsList:
                        guideDict[totalcount][sample][annotation] += 1


def cycleOverCount(totalcount):
    # cycle over count on the dictionary to create plots
    for elem in guideDict[totalcount]:
        generatePlot(guide, guideDict[totalcount][elem],
                     motifDict[totalcount][elem], totalcount, elem)


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
    guide = guide.strip().replace('N', '')
    guideDict = guideDictCreation()
    motifDict = motifDictCreation(guide)
    fillDict(guide, guideDict, motifDict)
    print('Starting image creation for guide: ', guide)
    # keyList = guideDict[0].keys()
    # print(keyList)
    if __name__ == '__main__':
        pool = multiprocessing.Pool(threads)
        for totalcount in range(0, 9):
            for elem in guideDict[totalcount]:
                pool.apply_async(generatePlot, args=(
                    guide, guideDict[totalcount][elem], motifDict[totalcount][elem], totalcount, elem))
        pool.close()
        pool.join()
    # for totalcount in range(0,9):
    #     p = Pool(9)
    #     p.map(cycleOverCount, [0,1,2,3,4,5,6,7,8])
        # for elem in guideDict[totalcount]:
        # generatePlot(guide, guideDict[totalcount][elem],motifDict[totalcount][elem],totalcount,elem)
    # p = Pool(20)
    # p.map(cycleOverCount, keyList)
    # p.map(cycleOverCount, foe)

    # for totalcount in range(0,9):
        # for elem in guideDict[totalcount]:
        # generatePlot(guide, guideDict[totalcount][elem],motifDict[totalcount][elem],totalcount,elem)
