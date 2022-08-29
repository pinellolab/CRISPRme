#!/usr/bin/env python

import sys
import time
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
# from upsetplot import generate_counts
from upsetplot import UpSet
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

inSamplesID = open(sys.argv[1], 'r')  # read samplesIDs of targets to analyze
inTargets = open(sys.argv[2], 'r')  # read targets
# low limit in CFD difference threshold (risk score) to include the target
CFDThreshold = float(sys.argv[3])
# low limit in CFD value to include the target
CFDabsLowest = float(sys.argv[4])

sampleDict = dict()  # dict for all samples
superPopDict = dict()  # dict for all superpopulations


def printUpsetPlot(targetsDF):
    figu = plt.figure()
    targetsDF = targetsDF.groupby(list(sampleDict)).size()
    # print(targetsDF)
    # upset = UpSet(targetsDF, facecolor='black',sort_by='cardinality',sort_categories_by=None, element_size=60)
    upset = UpSet(targetsDF, facecolor='black', element_size=60)
    upset.plot(fig=figu)
    plt.title(
        'Intersection plot for all superpopulations with diffCFD >= '+str(CFDThreshold)+' and absolute CFD value >='+str(CFDabsLowest))
    # plt.tight_layout()
    plt.savefig(sys.argv[2]+'_upsetplot_' +
                str(CFDThreshold)+'_'+str(CFDabsLowest)+'.pdf')


def printDensityPlot():
    # create figure and set axis
    plt.figure()
    for superPop in sampleDict:  # for each superpopulation
        andamenti = list()
        permutationList = list()
        for sample in superPopDict[superPop]:
            # append samples to list to permute
            permutationList.append(sample)
        print('DOING SUPERPOPULATION PLOT FOR: ', superPop)
        for permutation in range(0, 100):
            np.random.shuffle(permutationList)
            andamento = list()
            alreadyAddedTargets = set()
            for sample in permutationList:
                alreadyAddedTargets = alreadyAddedTargets.union(
                    superPopDict[superPop][sample])
                andamento.append(len(alreadyAddedTargets))
            andamenti.append(andamento)
        # read values to generate plot
        andamentiArray = np.array(andamenti)
        media = np.mean(andamentiArray, axis=0)
        standarddev = np.std(andamentiArray, axis=0)
        standarderr = standarddev/np.sqrt(len(list(sampleDict[superPop])))
        z_score = 1.96  # for confidence 95%
        lowerbound = media-(z_score*standarderr)
        upperbound = media+(z_score*standarderr)
        # allMedie.append(media)
        plt.plot(media, label=str(superPop))
        plt.fill_between(range(len(media)), lowerbound,
                         upperbound, alpha=0.10)

    plt.title('populations_with diffCFD >=' + str(CFDThreshold) +
              ' and CI '+str(95)+'%'+' and CFD score >='+str(CFDabsLowest))
    plt.xlabel('# Individuals')
    plt.ylabel('# Cumulative Targets')
    plt.legend()
    plt.tight_layout()
    plt.savefig(sys.argv[2]+'_allpop_with_diffCFD_'+str(CFDThreshold) +
                'and_CI_95_and_CFD_score_'+str(CFDabsLowest)+'.pdf')


def printSwarmPlot():
    totalList = list()
    for superPop in sampleDict:
        for pop in sampleDict[superPop]:
            for sample in sampleDict[superPop][pop]:
                for value in sampleDict[superPop][pop][sample][2]:
                    totalList.append([superPop, value])
    superPopDataFrame = pd.DataFrame(
        totalList, columns=['SuperPop', 'diffCFD'])
    plt.figure()

    ax = sns.violinplot(x='SuperPop', y='diffCFD', data=superPopDataFrame)
    for violin in ax.collections[::2]:
        violin.set_alpha(0.2)
    ax = sns.stripplot(x='SuperPop', y='diffCFD', data=superPopDataFrame)

    plt.title('Difference CFD distribution with diffCFD ' +
              str(CFDThreshold)+' and CFD>='+str(CFDabsLowest))
    plt.tight_layout()
    plt.savefig(sys.argv[2]+'_swarmplot_with_CFDdiff_'+str(CFDThreshold) +
                '_'+str(CFDabsLowest)+'.pdf')

    totalList = list()
    for superPop in sampleDict:
        for pop in sampleDict[superPop]:
            for sample in sampleDict[superPop][pop]:
                totalList.append(
                    [superPop, sampleDict[superPop][pop][sample][0]])
    superPopDataFrame = pd.DataFrame(
        totalList, columns=['SuperPop', 'Targets'])
    plt.figure()

    ax = sns.violinplot(x='SuperPop', y='Targets', data=superPopDataFrame)
    for violin in ax.collections[::2]:
        violin.set_alpha(0.2)
    ax = sns.stripplot(x='SuperPop', y='Targets', data=superPopDataFrame)

    plt.title('Difference CFD distribution with diffCFD ' +
              str(CFDThreshold)+' and CFD >='+str(CFDabsLowest))
    plt.tight_layout()
    plt.savefig(sys.argv[2]+'_swarmplot_with_targets_count_per_sample_' +
                str(CFDThreshold)+'_'+str(CFDabsLowest)+'.pdf')

    # for superPop in sampleDict:
    #     totalList = list()
    #     for pop in sampleDict[superPop]:
    #         for sample in sampleDict[superPop][pop]:
    #             totalList.append(
    #                 [pop[0:3], sampleDict[superPop][pop][sample][0]])
    #     superPopDataFrame = pd.DataFrame(
    #         totalList, columns=['Pop', 'Targets'])
    #     plt.figure()
    #     try:
    #         ax = sns.stripplot(x='Pop', y='Targets', data=superPopDataFrame)
    #     except:
    #         continue
    #     plt.title('Difference CFD distribution with diffCFD ' +
    #               str(CFDThreshold)+' and CFD >='+str(CFDabsLowest))
    #     plt.tight_layout()
    #     plt.savefig('swarmplot_with_targets_count_per_sample_' +
    #                 str(CFDThreshold)+'_'+str(CFDabsLowest)+'_in superpopulation_'+str(superPop)+'.pdf')

    # for superPop in sampleDict:
    #     totalList = list()
    #     for pop in sampleDict[superPop]:
    #         for sample in sampleDict[superPop][pop]:
    #             for value in sampleDict[superPop][pop][sample][2]:
    #                 totalList.append([pop[0:3], value])
    #     superPopDataFrame = pd.DataFrame(
    #         totalList, columns=['Pop', 'diffCFD'])
    #     plt.figure()
    #     try:
    #         ax = sns.stripplot(x='Pop', y='diffCFD', data=superPopDataFrame)
    #     except:
    #         continue
    #     plt.title('Difference CFD distribution with diffCFD ' +
    #               str(CFDThreshold)+' and CFD >='+str(CFDabsLowest))
    #     plt.tight_layout()
    #     plt.savefig('swarmplot_with_CFDdiff_'+str(CFDThreshold) +
    #                 '_'+str(CFDabsLowest)+'_in superpopulation_'+str(superPop)+'.pdf')


print('READING FILES AND DICTIONARY CREATION')
start = time.time()
inSamplesID.readline()
for line in inSamplesID:
    split = line.strip().split('\t')
    # for each sample in samplesID create a list to keep counter for targets and ID of target (line number) and superpop
    if split[2] not in sampleDict:
        sampleDict[split[2]] = dict()
    if split[1] not in sampleDict[split[2]]:
        sampleDict[split[2]][split[1]] = dict()
    # 0 is target count,set() is target line, list() is CFD diff HGDP_CFD-max(1000G_CFD,REF_CFD)
    sampleDict[split[2]][split[1]][split[0]] = [0, set(), list()]
    if split[2] not in superPopDict:
        superPopDict[split[2]] = dict()
    superPopDict[split[2]][split[0]] = set()


targetsDF = pd.DataFrame(columns=list(sampleDict))
targetDict = dict()

for lineNumber, target in enumerate(inTargets):
    split = target.strip().split('\t')
    if float(split[22]) < CFDThreshold:
        continue
    if float(split[20]) < CFDabsLowest:
        continue

    for elem in list(sampleDict):  # reset target dict
        targetDict[elem] = 0

    # extract all superpopulations from target
    superPop = split[-1].strip().split(',')
    # if len(superPop) > 6:
    #     print(target)

    for pop in superPop:  # for each superpopulation with the target
        # realSuperPop = pop.strip().split('_')[1]
        realSuperPop = pop
        # update targetdict to reflect the presence of the target in the superpopulation
        targetDict[realSuperPop] = 1
        # for each pop in superpop
        for pop in sampleDict[realSuperPop]:
            # for each sample in sample list in target line
            for sample in split[13].strip().split(','):
                try:
                    # count each target per superpopulation
                    superPopDict[realSuperPop][sample].add(lineNumber)
                    # count of each target per sample
                    sampleDict[realSuperPop][pop][sample][0] += 1
                    # collect the line number to identify target
                    sampleDict[realSuperPop][pop][sample][1].add(lineNumber)
                    # collect the CFD difference of target
                    sampleDict[realSuperPop][pop][sample][2].append(
                        float(split[-2]))
                except:
                    continue
    # do serie from dict of the target than append to dataframe
    targetSerie = pd.Series(data=targetDict)
    targetsDF = targetsDF.append(targetSerie, ignore_index=True)
print('DONE IN TIME ', time.time()-start)

print('CREATING SWARMPLOTs')
printSwarmPlot()
plt.close('all')


print('CREATING UPSET PLOT')
printUpsetPlot(targetsDF)
plt.close('all')

print('CREATING DENSITY PLOTs')
printDensityPlot()
plt.close('all')
