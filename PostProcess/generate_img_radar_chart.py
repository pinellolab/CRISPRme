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

from pandas.core.tools.numeric import to_numeric
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

guide = sys.argv[1]  # guide file used during search
try:
    guideDictFile = open(sys.argv[2])
    motifDictFile = open(sys.argv[3])
except:
    sys.exit()
mismatch = sys.argv[4]
bulge = sys.argv[5]
elem = sys.argv[6]
outDir = sys.argv[7]  # directory to output the figures
# threads = int(sys.argv[6])  # number of concurrent execution of image creation

web_server = False
population = False
# file_extension = 'pdf'
file_extension = 'png'

if '-ws' in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = 'png'


def generatePlot(guide, guideDict, motifDict, mismatch, bulge, source):

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
    # convert to int total column
    # print('prima', guideDataFrame)
    # guideDataFrame['Total'] = guideDataFrame['Total'].astype('int32')
    guideDataFrame = guideDataFrame.T
    # print('dopo', guideDataFrame)
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
    plt.xticks(angles[:-1], categories, color='black', size=14)

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
    guideDataFrame = guideDataFrame.T
    # guideDataFrame['Total'] = to_numeric(
    #     guideDataFrame['Total'], downcast='integer')
    # print('lamadonna', guideDataFrame)
    for elem in categories:
        transpose_list.append(list(guideDataFrame.loc[elem]))
        # transpose_list.append(
        #     [guideDataFrame.loc[elem, 'Total'], guideDataFrame.loc[elem, 'Percentage']])
    # print(transpose_list)
    templist = list()
    for couple in transpose_list:
        couple[0] = int(couple[0])
        templist.append(couple)
    # templist to convert into only the total column
    transpose_list = templist

    plt.axis('off')
    table = plt.table(cellText=transpose_list, rowLabels=categories, colLabels=['Total', 'Percentage'],
                      loc='best', colWidths=[0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(13)

    totalMotif = [0]*len(guide)
    for count in range(len(guide)):
        for nuc in motifDict:
            totalMotif[count] += motifDict[nuc][count]

    maxmax = max(totalMotif)
    for count in range(len(guide)):
        for nuc in motifDict:
            if maxmax != 0:
                motifDict[nuc][count] = float(
                    motifDict[nuc][count]/float(maxmax))
                # motifDict[nuc][count] = float(
                #     str(float(motifDict[nuc][count])/float(maxmax))[0:5])

    # ind = np.arange(0, len(guide), 1) + 0.15
    ind = np.arange(0, len(guide), 1)
    width = 0.7  # the width of the bars: can also be len(x) sequence

    motif = plt.subplot(2, 1, 2, frameon=False)

    A = np.array(motifDict['A'], dtype=float)
    C = np.array(motifDict['C'], dtype=float)
    G = np.array(motifDict['G'], dtype=float)
    T = np.array(motifDict['T'], dtype=float)
    RNA = np.array(motifDict['RNA'], dtype=float)
    DNA = np.array(motifDict['DNA'], dtype=float)

    p1 = plt.bar(ind, A, width, align='center')
    p2 = plt.bar(ind, C, width, bottom=A, align='center')
    p3 = plt.bar(ind, G, width, bottom=A+C, align='center')
    p4 = plt.bar(ind, T, width, bottom=C+G+A, align='center')
    p5 = plt.bar(ind, RNA, width, bottom=C+G+A+T, align='center')
    p6 = plt.bar(ind, DNA, width, bottom=C+G+A+T+RNA,
                 align='center')

    # plt.xlim(0, len(guide))
    # strArray = np.array([list(guide)])
    plt.xticks(ticks=ind, labels=list(guide))

    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
               ('A', 'C', 'G', 'T', 'bRNA', 'bDNA'), fontsize=13, loc='upper left', ncol=6)

    # strArray = np.array([list(guide)])
    # table = plt.table(cellText=strArray, loc='bottom',
    #                   cellLoc='center', rowLoc='bottom')
    # table.auto_set_font_size(False)
    # table.set_fontsize(12)
    # # table.scale(1, 1.6)
    # table.xticks = ([])
    # table.yticks = ([])
    # plt.xticks

    plt.suptitle(str(mismatch)+" Mismatches + "+str(bulge)+" Bulge "+str(source),
                 horizontalalignment='center', color='black', size=25)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85, bottom=0.05, left=0.05,
                        right=0.95, wspace=0.1)
    plt.savefig(outDir+"/summary_single_guide_" + str(guide) + "_" + str(mismatch) +
                "."+str(bulge) + '_' + str(source) + "." + file_extension, format=file_extension)

    plt.close('all')


guide = guide.strip()
guideDict = json.load(guideDictFile)
motifDict = json.load(motifDictFile)
if __name__ == '__main__':
    # skip creation if no target in global category
    if guideDict[mismatch][bulge]['TOTAL']['General'] == 0:
        sys.exit()

    generatePlot(guide, guideDict[mismatch][bulge][elem],
                 motifDict[mismatch][bulge][elem], mismatch, bulge, elem)
