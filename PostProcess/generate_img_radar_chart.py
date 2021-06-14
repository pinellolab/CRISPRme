#!/usr/bin/env python

# Libraries
from operator import itemgetter
from os.path import isfile, join
from os import listdir
import os
from posixpath import split
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

# afont = {'fontname': 'Helvetica'}

guide = sys.argv[1]  # guide file used during search
try:
    guideDictFile = open(sys.argv[2])
    motifDictFile = open(sys.argv[3])
except:
    sys.exit()
mismatch = int(sys.argv[4])
bulge = int(sys.argv[5])
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


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def generatePlot(guide, guideDict, motifDict, mismatch, bulge, source):

    # check if no targets are found for that combination source/totalcount and skip the execution
    if guideDict['General'] == 0:
        return

    total = mismatch+bulge  # count total of mismatch and bulge requested
    titlesize = 18
    fontsize = 17

    percentage_list = []
    for elem in guideDict:
        if float(guideDict['General']) != 0:
            percentage_list.append(
                round(float(guideDict[elem])*100/float(guideDict['General']), 2))
        else:
            percentage_list.append(round(float(0), 0))

    guideDataFrame = pd.DataFrame.from_dict(guideDict, orient='index')

    guideDataFrame['Percentage'] = percentage_list
    guideDataFrame.columns = ['Count', 'Percentage']

    try:
        # dataframe for table creation
        dataframe_table = guideDataFrame.loc[['General', 'three_prime_UTR', 'five_prime_UTR',
                                              'exon', 'CDS', 'gene', 'DNase-H3K4me3', 'CTCF-only', 'dELS', 'pELS', 'PLS']]
        # dataframe_table = guideDataFrame.loc[[
        #     'General', '', 'pELS', 'dELS', 'PLS', '', 'exon', 'gene', 'CDS', '', '']]
        dataframe_table.rename(
            index={'CTCF-only': 'CTCF', 'General': 'Total', 'three_prime_UTR': "3'UTR", 'five_prime_UTR': "5'UTR"}, inplace=True)
        dataframe_table.sort_values(
            by=['Percentage'], ascending=False, inplace=True)
    except:
        dataframe_table = guideDataFrame.loc[['General']]
        dataframe_table.rename(index={'General': 'Total'}, inplace=True)

    # print(dataframe_table)

    # dataframe for radar chart, drop total column
    if len(list(dataframe_table.T)[0:]) > 1:
        dataframe_radar_chart = dataframe_table.drop(['Total'])
    else:
        dataframe_radar_chart = dataframe_table

    dataframe_radar_chart = dataframe_radar_chart.T

    categories_radar_chart = list(dataframe_radar_chart)[0:]
    categories_table = list(dataframe_table.T)[0:]

    count_radar_chart_categories = len(categories_radar_chart)

    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values_radar_chart = dataframe_radar_chart.loc['Percentage'].values.flatten(
    ).tolist()
    values_radar_chart += values_radar_chart[:1]

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles_radar_chart = [n / float(count_radar_chart_categories) * 2 *
                          pi for n in range(count_radar_chart_categories)]
    angles_radar_chart += angles_radar_chart[:1]

    # Initialise the spider plot
    ax = plt.subplot(2, 2, 1, polar=True)

    for label, rot in zip(ax.get_xticklabels(), angles_radar_chart):
        if (rot == 0):
            label.set_horizontalalignment("center")
        if (rot > 0):
            label.set_horizontalalignment("left")
        if (rot > 3):
            label.set_horizontalalignment("center")
        if (rot > 4):
            label.set_horizontalalignment("right")

    # Draw one axe per variable + add labels labels yet
    # plt.xticks(angles[:-1], categories, color='black', size=fontsize)
    plt.xticks(angles_radar_chart[:-1], categories_radar_chart,
               color='black', size=fontsize)

    # Draw ylabels
    # # # Draw ylabels
    ax.set_rlabel_position(0)
    max_value_radar_chart = round(max(values_radar_chart))
    # round to upper 10 multiple
    while int(max_value_radar_chart % 10):
        max_value_radar_chart = max_value_radar_chart + 1

    radar_chart_yticks = [elem for elem in range(0, max_value_radar_chart, 10)]
    radar_chart_yticks_labels = [
        str(elem) for elem in range(0, max_value_radar_chart, 10)]
    # plt.yticks([0, 0.25, 0.50, 0.75], ["0", "0.25", "0.50", "0.75"], color="black", size=12)
    # plt.yticks([0, 25, 50, 75], ["0", "25", "50", "75"], color="black", size=fontsize-2)
    plt.yticks(radar_chart_yticks, radar_chart_yticks_labels,
               color="black", size=fontsize)
    # plt.ylim(0, 1)
    plt.ylim(0, max_value_radar_chart)

    # Fill area
    ax.fill(angles_radar_chart, values_radar_chart, 'b', alpha=0.1)

    # # # offset posizione y-axis
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    # Plot data
    ax.plot(angles_radar_chart, values_radar_chart,
            linewidth=1, linestyle='solid')

    plt.subplot(2, 2, 2)
    transpose_list = list()
    # dataframe_table = dataframe_table.T
    for elem in categories_table:
        transpose_list.append(list(dataframe_table.loc[elem]))
    templist = list()
    for couple in transpose_list:
        couple[0] = int(couple[0])
        templist.append(couple)
    # templist to convert into only the total column
    transpose_list = templist

    # print('transp list', transpose_list)

    plt.axis('off')
    table = plt.table(cellText=transpose_list, rowLabels=categories_table, colLabels=['Count', 'Percentage'],
                      loc='best', colWidths=[0.25, 0.35])
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1, 1.5)

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

    plt.xticks(ticks=ind, labels=list(guide), size=fontsize)
    plt.yticks(size=fontsize)

    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
               ('A', 'C', 'G', 'T', 'bRNA', 'bDNA'), fontsize=fontsize, loc='upper left', ncol=6)
    plt.title(
        'Mismatch and bulge distribution for targets with up to ' + str(total) + ' mismatches and/or bulges', color='black', size=titlesize)

    plt.suptitle('Targets with up to ' + str(total) + ' mismatches and/or bulges by ENCODE/GENCODE annotations',
                 horizontalalignment='center', color='black', size=titlesize)

    plt.tight_layout()

    plt.savefig(outDir+'/summary_single_guide_' + str(guide) + '_' + str(mismatch) +
                '.' + str(bulge) + '_' + str(source) + '.ENCODE+GENCODE.' + file_extension, format=file_extension)

    plt.close('all')


guide = guide.strip()
guideDict = json.load(guideDictFile)
motifDict = json.load(motifDictFile)
if __name__ == '__main__':
    total = mismatch+bulge
    total = str(total)
    # skip creation if no target in global category
    # if guideDict[mismatch][bulge]['TOTAL']['General'] == 0:
    if guideDict[total]['General'] == 0:
        sys.exit()
    # try:
        # generatePlot(guide, guideDict[mismatch][bulge][elem], motifDict[mismatch][bulge][elem], mismatch, bulge, elem)
    generatePlot(guide, guideDict[total],
                 motifDict[total], mismatch, bulge, elem)
    # except:
    #     sys.exit()
