#!/usr/bin/env python


'''
Not used
'''
# Esempio di input:
#python 3 new_radar_chart GUIDE -mm 4 -profile profile.xls -extprofile extprofile.xls -files ... ... ... ...  -> creo radar chart
#python 3 new_radar_chart GUIDE -mm 4 -profile profile.xls -extprofile extprofile.xls -files ... ... ... ... -sumref ... -sumenr ... -> creo radar chart + barplot di quella guida
#python 3 new_radar_chart -mm 4  -sumref ... -sumenr ... -> creo barplot totale
# Input:
# guide, mm, profile, extprofile, annotation_file, summaryone, summarytwo, gecko profile, gecko exons, gecko introns, gecko promoter, gecko dnase, geckoctcf

#BUG 428 quando guida non ha trovato targets
# Libraries
import math
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import pandas as pd
from math import pi
import scipy.spatial.distance as sp
import numpy as np
import sys
from itertools import islice
import glob
import warnings
import os
from os import listdir
from os.path import isfile, join
warnings.filterwarnings("ignore")

plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

only_radar = False
radar_barplot = False
only_barplot = False

guide = sys.argv[1]
if '-' in sys.argv[2]:
    missmatch = sys.argv[2].split('-')
    lowermm = int(missmatch[0])
    uppermm = int(missmatch[1])
else:
    missmatch = int(sys.argv[2])
    lowermm = 0
    uppermm = missmatch

profile_file = sys.argv[3]
ext_profile_file = sys.argv[4]
count_dir = sys.argv[5]
summary_one = sys.argv[6]   
summary_two = sys.argv[7]   #is the same as count_dir

try:
    gecko_profile_file = sys.argv[8]
    gecko_exons = sys.argv[9]
    gecko_introns = sys.argv[10]
    gecko_promoter = sys.argv[11]
    gecko_dnase = sys.argv[12]
    gecko_ctcf = sys.argv[13]
except:
    pass

if guide == 'no':
    only_barplot = True
elif summary_one != 'no' and summary_two != 'no':
    radar_barplot = True
else:
    only_radar = True


# lettura file
with open(profile_file) as profile:
    num_guides = profile.read().strip().split('\n')[1:]
    num_guides = len (num_guides)
with open(count_dir) as annotation_file:
    content = annotation_file.read().strip()
    onlyfiles = content.split('-Summary_Total\n')[0].split('-')
    summary_two = content.split('-Summary_Total\n')[1]
    if guide == 'no':
        summary_two = summary_two.strip().split('-')[0].strip().split('\n')
       
    else:
        summary_two = summary_two.strip().split('-Summary_')[1:]
        summary_two = [s.strip() for s in summary_two if guide in s]
        if not summary_two:
            print('Warning! The selected guide has no annotated results')
            quit()
        summary_two = summary_two[0].split('\n')[1:]
        

if summary_one != 'no':
    with open(summary_one) as annotation_file:
        content = annotation_file.read().strip()
        summary_one = content.split('-Summary_Total\n')[1]
        if guide == 'no':
            summary_one = summary_one.strip().split('-')[0].strip().split('\n')
        else:
            summary_one = summary_one.strip().split('-Summary_')[1:]
            summary_one = [s.strip() for s in summary_one if guide in s]
            summary_one = summary_one[0].split('\n')[1:]

count_files = []
inGuidesProfileExtended = open(ext_profile_file, 'r')
inGuidesProfile = open(profile_file, 'r')
onlyfiles.sort()

for i in onlyfiles:
    if i == '':
        continue
    i = i.strip()
    element_name = i[:i.find('\n')] #i[i[:i.find('Count')].rfind('.') +1 : i.find('Count')]
    count_files.append([element_name, i.split('\n')[1:]])

check_annotation_name = []
for i in count_files:
    check_annotation_name.append(i[0].lower())

if check_annotation_name != ['ctcf', 'dnase', 'exon', 'intron', 'promoter'] and gecko_profile_file != "no":     #NOTE if gecko annotations are updated, update the list adding the new annotations in lexico order
    print('Warning! Option \'-gecko\' can be used on files annotated with \'ctcf\', \'dnase\', \'exon\', \'intron\', \'promoter\'\nRemoving \'-gecko\' option')
    gecko_profile_file = 'no'
#count_files is now ['Element', [first_guide, second_guide...]] ...

# if exonsCountFile != "no":
#     inExonsCountFile = open(exonsCountFile, "r")
# if intronsCountFile != "no":
#     inIntronsCountFile = open(intronsCountFile, "r")
# if promotersCountFile != "no":
#     inPromotersCountFile = open(promotersCountFile, "r")
# if DNAseCountFile != "no":
#     inDNAseCountFile = open(DNAseCountFile, "r")
# if DNAseCountFile != "no":
#     inCTCFCountFile = open(CTCFCountFile, "r")

# lists for data storing and analysis
guidesExtendedProfile = []

# global counting for annotation types

profileMissmatchGlobal = []
# exonsMissmatchGlobal = []
# intronsMissmatchGlobal = []
# promotersMissmatchGlobal = []
# DNAseMissmatchGlobal = []
# CTCFMissmatchGlobal = []

arraySummaryCountOne = []
arraySummaryCountTwo = []


if summary_one != "no" and summary_two != "no":
    inSummaryCountOne = summary_one #open(summary_one, "r")
    inSummaryCountTwo = summary_two #open(summary_two, "r")
    
    x_ticks_labels = []
    n_annotation = 0
    mms_total = 0      #len of the array exon   0   0   0   0   0   0   0, taken from [1:]
    for line in inSummaryCountOne:
        x = line.strip().split('\t')
        arraySummaryCountOne.append(tuple(x[1:]))
        x_ticks_labels.append(x[0])
        n_annotation = n_annotation + 1
        mms_total = len(x[1:])
    

    for line in inSummaryCountTwo:
        x = line.strip().split('\t')
        arraySummaryCountTwo.append(tuple(x[1:]))
    arraySummaryCountOne = np.array(arraySummaryCountOne, dtype=int)
    arraySummaryCountOne.shape = (n_annotation, mms_total)

    arraySummaryCountTwo = np.array(arraySummaryCountTwo, dtype=int)
    arraySummaryCountTwo.shape = (n_annotation, mms_total)

    percentageGain = []
    for row in range(0, n_annotation):
        for col in range(0, uppermm+1):
            res = max((arraySummaryCountOne[row, col]/arraySummaryCountTwo[row, col]),
                      (arraySummaryCountTwo[row, col]/arraySummaryCountOne[row, col]))
            percentageGain.append(res)
    percentageGain = np.array(percentageGain, dtype=float)
    percentageGain.shape = (n_annotation, uppermm+1)

    intergenicGainOne = (np.sum(arraySummaryCountOne, axis=0))
    intergenicGainTwo = (np.sum(arraySummaryCountTwo, axis=0))
    intergenicGainOne.shape = (1, mms_total)
    intergenicGainTwo.shape = (1, mms_total)
    intergenicGainOne = intergenicGainOne - arraySummaryCountOne[0, :]
    intergenicGainTwo = intergenicGainTwo - arraySummaryCountTwo[0, :]

    # p1 = plt.bar(1, percentageGain[0, uppermm-1])
    # p2 = plt.bar(2, percentageGain[1, uppermm-1])
    # p3 = plt.bar(3, percentageGain[2, uppermm-1])
    # p4 = plt.bar(4, percentageGain[3, uppermm-1])
    # p5 = plt.bar(5, percentageGain[4, uppermm-1])
    # p6 = plt.bar(6, percentageGain[5, uppermm-1])
    # p7 = plt.bar(7, max((intergenicGainOne[0, uppermm-1]/intergenicGainTwo[0, uppermm-1]),
    #                     (intergenicGainTwo[0, uppermm-1]/intergenicGainOne[0, uppermm-1])))

    ind = np.arange(0, n_annotation, 1)
    # print('math ceil', np.arange(0, max(arraySummaryCountTwo[:, uppermm]) + math.ceil(max(
    #     arraySummaryCountTwo[:, uppermm])/10), 
    #     math.ceil(max(arraySummaryCountTwo[:, uppermm])/5)))  #BUG if max(arraySummaryCountTwo[:, uppermm], then arange (0,0,0)
    no_result = False
    try:
        y_range = np.arange(0, max(arraySummaryCountTwo[:, uppermm]) + math.ceil(max(
            arraySummaryCountTwo[:, uppermm])/10), math.ceil(max(arraySummaryCountTwo[:, uppermm])/5))
    except:
        y_range = np.arange(0,1,1)
        no_result = True
    width = 0.5

    p1 = plt.bar(
        ind, arraySummaryCountOne[:, uppermm], width, color='#67a9cf', align='edge')
    p2 = plt.bar(ind, (arraySummaryCountTwo[:, uppermm]-arraySummaryCountOne[:, uppermm]),
                 width, bottom=arraySummaryCountOne[:, uppermm], color='#ef8a62', align='edge')

    plt.legend((p1[0], p2[0]), ('Reference Genome',
                                'Enriched Genome'), fontsize=30)

    # plt.xlim(0, len(string))
    # plt.set_ylim([1, 1.5])
    plt.title('Relative Increase Enriched/Reference Genome with ' +
              str(uppermm) + ' Mismatches', size=25)
    # plt.xlabel('Annotations')
    if no_result:
        plt.annotate('No targets found with ' + str(missmatch)  + ' mismatches', [1.35,0], size = 22) #NOTE with 0-mm print only the mm pdf; 1.35 modificare se cambia la str
    else:
        for k in range(0, n_annotation):
            plt.annotate('%.2fx' % percentageGain[k, uppermm], [
                        k+0.05, arraySummaryCountTwo[k, uppermm]+(max(arraySummaryCountTwo[:, uppermm])/100)], size=22)
        # plt.ylim([0, max(arraySummaryCountTwo[:, uppermm])+2000],size=25)

    
    plt.xticks(ind+0.25, x_ticks_labels, size=25)
    plt.yticks(y_range, size=22)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.06, left=0.1, right=0.99)
    if guide != 'no':
        plt.savefig("summary_histogram_" + str(guide) + '_' + str(uppermm) + 
                    "mm" + ".pdf", format="pdf")
    else:
        plt.savefig("summary_histogram_" + str(uppermm) + 
                    "mm" + ".pdf", format="pdf")
    if guide != 'no':
        plt.savefig("summary_histogram_" + str(guide) + '_' + str(uppermm) + 
                    "mm" + ".png", format="png")
    else:
        plt.savefig("summary_histogram_" + str(uppermm) + 
                    "mm" + ".png", format="png")


if guide != 'no':
    # reading extendend profile to obtain results over mismatches counts
    for line in inGuidesProfileExtended:
        if ">" + guide in line:
            # print(line)
            next(inGuidesProfileExtended)
            # line=inGuidesProfileExtended.readline()
            for ciao in range(0, uppermm+1):
                line = inGuidesProfileExtended.readline()
                count = 0
                x = line.split('\t')
                guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
                                              x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
                for line in inGuidesProfileExtended:
                    if count < 6:
                        line = line.rstrip()
                        x = line.split('\t')
                        #y = str(x[20]).split('\n')
                        guidesExtendedProfile.append((x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],
                                                      x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]))
                        count += 1
                    else:
                        break
            break
    arrayguidesExtendedProfile = np.array(guidesExtendedProfile, dtype=int)
    
    arrayguidesExtendedProfile.shape = (7*((uppermm-0)+1), 20)

    # reading profile file to obtain results for every mismatch count in the general profile
    guide_len = len(guide)
    next(inGuidesProfile)
    for line in inGuidesProfile:
        line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
            "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
        x = line.split('\t')
        profileMissmatchGlobal.append((x[guide_len+3], x[guide_len+4], x[guide_len+5], x[guide_len+6],
                                       x[guide_len+7], x[guide_len+8], x[guide_len+9], x[guide_len+10], x[guide_len+11], x[guide_len+12]))
        if str(x[0]) == guide:
            arrayprofileMissmatch = np.array((x[guide_len+3], x[guide_len+4], x[guide_len+5], x[guide_len+6],
                                              x[guide_len+7], x[guide_len+8], x[guide_len+9], x[guide_len+10], x[guide_len+11], x[guide_len+12]), dtype=int)

    arrayexonsMissmatch = []
    exonsMissmatchGlobal = []
    for pos,element_file in enumerate(count_files):
        exonsMissmatchGlobal.append([])
        arrayexonsMissmatch.append([])
        for line in element_file[1]:
            line = line.strip()
            #line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
            #    "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            exonsMissmatchGlobal[pos].append(tuple(x[1:]))
            if str(x[0]) == guide:
                arrayexonsMissmatch[pos].append(np.array( tuple(x[1:]), dtype=int))

    if(gecko_profile_file != "no"):   #NOTE if new annotations files are provided, add them in lexicograph order, following the structure for ctcf, dnase etc etc  
        inGuidesProfile = open(gecko_profile_file, "r")
        inExonsCountFile = open(gecko_exons, "r")
        inIntronsCountFile = open(gecko_introns, "r")
        inPromotersCountFile = open(gecko_promoter, "r")
        inDNAseCountFile = open(gecko_dnase, "r")
        inCTCFCountFile = open(gecko_ctcf, "r")

        profileMissmatchGlobal.clear()
        # exonsMissmatchGlobal.clear()
        # intronsMissmatchGlobal.clear()
        # promotersMissmatchGlobal.clear()
        # DNAseMissmatchGlobal.clear()
        # CTCFMissmatchGlobal.clear()
        exonsMissmatchGlobal = []
        # GECKO 
        # reading profile file to obtain results for every mismatch count in the general profile
        guide_len = len(guide)
        next(inGuidesProfile)
        for line in inGuidesProfile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            profileMissmatchGlobal.append((x[guide_len+3], x[guide_len+4], x[guide_len+5], x[guide_len+6],
                                        x[guide_len+7], x[guide_len+8], x[guide_len+9], x[guide_len+10], x[guide_len+11], x[guide_len+12]))

        # reading every count file to obtain results for the ecdf and percentile count for annotated type
        tmp = []
        for line in inCTCFCountFile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            tmp.append(
                (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))
        exonsMissmatchGlobal.append(tmp)
        tmp = []
        for line in inDNAseCountFile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            tmp.append(
                (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))
        exonsMissmatchGlobal.append(tmp)
        tmp = []
        for line in inExonsCountFile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            tmp.append(
                (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))
        exonsMissmatchGlobal.append(tmp)
        tmp = []
        for line in inIntronsCountFile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            tmp.append(
                (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))
        exonsMissmatchGlobal.append(tmp)
        tmp = []
        for line in inPromotersCountFile:
            line += "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t" + \
                "0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"
            x = line.split('\t')
            tmp.append(
                (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]))
        exonsMissmatchGlobal.append(tmp)
               
        
        
    if profile_file != "no":
        arrayexonsMissmatchGlobal = []
        arrayprofileMissmatchGlobal = np.array(profileMissmatchGlobal, dtype=int)
        for elem in exonsMissmatchGlobal:
            arrayexonsMissmatchGlobal.append(np.array(elem, dtype = int))
        # arrayexonsMissmatchGlobal = np.array(exonsMissmatchGlobal, dtype=int)
        # arrayintronsMissmatchGlobal = np.array(intronsMissmatchGlobal, dtype=int)
        # arraypromotersMissmatchGlobal = np.array(promotersMissmatchGlobal, dtype=int)
        # arrayDNAseMissmatchGlobal = np.array(DNAseMissmatchGlobal, dtype=int)
        # arrayCTCFMissmatchGlobal = np.array(CTCFMissmatchGlobal, dtype=int)
        distances = []
        test_distances = []
        Generaldistance = [[0 for x in range(0, uppermm+1)]
                        for y in range(np.size(arrayprofileMissmatchGlobal, 0))]
        for elem in exonsMissmatchGlobal:
            distances.append([[0 for x in range(0, uppermm+1)] for y in range(np.size(elem, 0))]  )
            test_distances.append([[0 for x in range(0, uppermm+1)] for y in range(np.size(elem, 0))]  )
        
        # CTCFdistance = [[0 for x in range(0, uppermm+1)]
        #                 for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        # Exonsdistance = [[0 for x in range(0, uppermm+1)]
        #                  for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        # Intronsdistance = [[0 for x in range(0, uppermm+1)]
        #                    for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        # Promotersdistance = [[0 for x in range(0, uppermm+1)]
        #                      for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        # Generaldistance = [[0 for x in range(0, uppermm+1)]
        #                    for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        # DNAsedistance = [[0 for x in range(0, uppermm+1)]
        #                  for y in range(np.size(arrayCTCFMissmatchGlobal, 0))]
        for i in range(0, uppermm+1):
            
            for pos, eleme in enumerate(arrayexonsMissmatchGlobal):
                arrayexonsMissmatchGlobal[pos][:,i] = np.sort( eleme[:,i], axis=None)

            # arrayCTCFMissmatchGlobal[:, i] = np.sort(
            #     arrayCTCFMissmatchGlobal[:, i], axis=None)
            arrayprofileMissmatchGlobal[:, i] = np.sort(
                arrayprofileMissmatchGlobal[:, i], axis=None)

            # arrayintronsMissmatchGlobal[:, i] = np.sort(
            #     arrayintronsMissmatchGlobal[:, i], axis=None)
            # arraypromotersMissmatchGlobal[:, i] = np.sort(
            #     arraypromotersMissmatchGlobal[:, i], axis=None)
            # arrayDNAseMissmatchGlobal[:, i] = np.sort(
            #     arrayDNAseMissmatchGlobal[:, i], axis=None)
            # arrayexonsMissmatchGlobal[:, i] = np.sort(
            #     arrayexonsMissmatchGlobal[:, i], axis=None)
        
            for k in range(0, np.size(arrayexonsMissmatchGlobal[0], 0)):
                for pos, elem in enumerate(arrayexonsMissmatchGlobal):
                    distances[pos][k][i] = abs (arrayexonsMissmatch[pos][0][i] - elem[k,i])
                
                
                for pos, elem in enumerate(arrayexonsMissmatchGlobal):
                    test_distances[2][k][i] = abs (arrayexonsMissmatch[2][0][i] - arrayexonsMissmatchGlobal[0][k,i])
                
                # CTCFdistance[k][i] = abs(
                #     arrayCTCFMissmatch[i]-arrayCTCFMissmatchGlobal[k, i])
                # Exonsdistance[k][i] = abs(
                #     arrayexonsMissmatch[i]-arrayexonsMissmatchGlobal[k, i])
                # Intronsdistance[k][i] = abs(
                #     arrayintronsMissmatch[i]-arrayintronsMissmatchGlobal[k, i])
                # Promotersdistance[k][i] = abs(
                #     arraypromotersMissmatch[i]-arraypromotersMissmatchGlobal[k, i])
                Generaldistance[k][i] = abs(
                    arrayprofileMissmatch[i]-arrayprofileMissmatchGlobal[k, i])
                # DNAsedistance[k][i] = abs(
                #     arrayDNAseMissmatch[i]-arrayDNAseMissmatchGlobal[k, i])
        array_distances = []
        for pos, elem in enumerate(distances):
            array_distances.append(np.array(elem, dtype=int))
        # arrayCTCFdistance = np.array(CTCFdistance, dtype=int)
        # arrayExonsFdistance = np.array(Exonsdistance, dtype=int)
        # arrayIntronsdistance = np.array(Intronsdistance, dtype=int)
        # arrayPromotersdistance = np.array(Promotersdistance, dtype=int)
        arrayProfiledistance = np.array(Generaldistance, dtype=int)
        # arrayDNAsedistance = np.array(DNAsedistance, dtype=int)

    if profile_file != "no":
        # SINGLE MISMATCH COUNT
        if len(sys.argv[2]) == 1:       #mm inserted is single number 
            # Set data
            data_for_df = {'group': ['A'], 'General': [np.argmin(arrayProfiledistance[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0)]}
            for pos, el in enumerate(array_distances):
                data_for_df[count_files[pos][0]] = [np.argmin(el[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0) ]
            # df = pd.DataFrame({
            #     'group': ['A'],
            #     'Exons': [np.argmin(arrayExonsFdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)],
            #     'General': [np.argmin(arrayProfiledistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)],
            #     'Introns': [np.argmin(arrayIntronsdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)],
            #     'Promoters': [np.argmin(arrayPromotersdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)],
            #     'DNAse': [np.argmin(arrayDNAsedistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)],
            #     'CTCF': [np.argmin(arrayCTCFdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0)]
            # })
            df = pd.DataFrame(data_for_df)
            # number of variable
            categories = list(df)[1:]
            N = len(categories)

            # We are going to plot the first line of the data frame.
            # But we need to repeat the first value to close the circular graph:
            values = df.loc[0].drop('group').values.flatten().tolist()
            values += values[:1]

            # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
            angles = [n / float(N) * 2 * pi for n in range(N)]
            angles += angles[:1]

            # Initialise the spider plot
            ax = plt.subplot(2, 2, 1, polar=True)
            # ax=plt.subplot(1, 1, 1, polar=True)
            # plt.title('RADAR CHART')

            #labels = ['CTCF', 'DNAse', 'Exons', 'General', 'Introns', 'Promoters']
            labels = list(df.columns.values[1:])
            # Draw one axe per variable + add labels labels yet
            plt.xticks(angles[:-1], labels, color='black', size=18)
            for label, rot in zip(ax.get_xticklabels(), angles):
                if (rot == 0):
                    label.set_horizontalalignment("center")
                if (rot > 0):
                    label.set_horizontalalignment("left")
                if (rot > 3):
                    label.set_horizontalalignment("center")
                if (rot > 4):
                    label.set_horizontalalignment("right")

            # offset posizione y-axis
            ax.set_theta_offset(pi / 2)
            ax.set_theta_direction(-1)

            # Draw ylabels
            ax.set_rlabel_position(0)
            plt.yticks([0, 0.25, 0.50, 0.75, 1], ["0", "0.25",
                                                "0.50", "0.75"], color="black", size=18)
            plt.ylim(0, 1)

            # Plot data
            ax.plot(angles, values, linewidth=1, linestyle='solid')

            # Fill area
            ax.fill(angles, values, 'b', alpha=0.1)

            columns = ('Position', '# Targets')
            #rows = ('General', 'Exons', 'Introns', 'Promoters', 'DNAse', 'CTCF')
            rows = ['General']
            for el in count_files:
                rows.append(el[0])
            rows = list(rows)
            table_distances= []
            for pos, el in enumerate(array_distances):
                table_distances.append(round(np.argmin(
                el[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0), 2))
            # exons_distance = round(np.argmin(
            #     arrayExonsFdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0), 2)
            # introns_distance = round(np.argmin(
            #     arrayIntronsdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0), 2)
            general_distance = round(np.argmin(
                arrayProfiledistance[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0), 2)
            # promoters_distance = round(np.argmin(
            #     arrayPromotersdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0), 2)
            # dnase_distance = round(np.argmin(
            #     arrayDNAsedistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0), 2)
            # ctcf_distance = round(np.argmin(
            #     arrayCTCFdistance[:, uppermm])/np.size(arrayCTCFMissmatchGlobal, 0), 2)

            temp_vstack = np.vstack([f[0][uppermm] for f in arrayexonsMissmatch])
            offtarget_data = np.vstack((arrayprofileMissmatch[uppermm], temp_vstack))
            
            #offtarget_data = np.vstack((arrayprofileMissmatch[uppermm], arrayexonsMissmatch[uppermm], arrayintronsMissmatch[uppermm],
            #                            arraypromotersMissmatch[uppermm], arrayDNAseMissmatch[uppermm], arrayCTCFMissmatch[uppermm]))
            # distance_data = np.vstack((general_distance, exons_distance,
            #                            introns_distance, promoters_distance, dnase_distance, ctcf_distance))
            temp_vstack = np.vstack([f for f in table_distances])
            distance_data = np.vstack((general_distance, temp_vstack))
            table_data = np.concatenate((distance_data, offtarget_data), axis=1)

            plt.subplot(2, 2, 2)
            table = plt.table(cellText=table_data, rowLabels=rows,
                            colLabels=columns, loc='center', colWidths=[0.35 for x in columns])
            table.auto_set_font_size(False)
            table.set_fontsize(18)
            table.scale(1, 3)
            plt.axis('off')
            plt.savefig("test.summary_single_guide_" + str(guide) +
                        "_"+str(uppermm) + "mm" + ".png", format="png")
            quit()

            datacount = arrayguidesExtendedProfile[missmatch*7] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            data = np.array(datacount, dtype=float)
            data = np.around(data, decimals=1)
            data.shape = (1, len(datacount))

            string = guide[0:20]
            strArray = np.array([list(string)])

            A = arrayguidesExtendedProfile[missmatch*7+1] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            C = arrayguidesExtendedProfile[missmatch*7+2] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            G = arrayguidesExtendedProfile[missmatch*7+3] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            T = arrayguidesExtendedProfile[missmatch*7+4] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            DNA = arrayguidesExtendedProfile[missmatch*7+5] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))
            RNA = arrayguidesExtendedProfile[missmatch*7+6] / \
                (max(arrayguidesExtendedProfile[missmatch*7]))

            # the x locations for the groups
            ind = np.arange(0, len(string), 1) + 0.15
            width = 0.7  # the width of the bars: can also be len(x) sequence

            motif = plt.subplot(2, 1, 2, frameon=False)
            # motif=plt.subplot(1,1,1)
            p1 = plt.bar(ind, A, width, color='red', align='edge')
            p2 = plt.bar(ind, C, width, color='blue', bottom=A, align='edge')
            p3 = plt.bar(ind, G, width, color='green', bottom=A+C, align='edge')
            p4 = plt.bar(ind, T, width, color='purple', bottom=C+G+A, align='edge')
            p5 = plt.bar(ind, DNA, width, color='magenta', bottom=C+G+A+T, align='edge')
            p6 = plt.bar(ind, RNA, width, color='gold', bottom=C+G+A+T+DNA, align='edge')
            plt.xlim(0, len(string))
            plt.xticks([])

            plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]),
                    ('A', 'C', 'G', 'T', 'D', 'R'), fontsize=18)

            table = plt.table(cellText=strArray, loc='bottom',
                            cellLoc='center', rowLoc='bottom')
            table.auto_set_font_size(False)
            table.set_fontsize(18)
            table.scale(1, 1.6)
            table.xticks = ([])
            table.yticks = ([])

            plt.suptitle(str(missmatch)+" Mismatches",
                        horizontalalignment='center', color='black', size=25)

            plt.tight_layout()
            plt.subplots_adjust(top=0.90, bottom=0.07, left=0.05,
                                right=0.99, wspace=0.12)

            plt.savefig("summary_single_guide_" + str(guide) +
                        "_"+str(uppermm) + "mm" + ".pdf", format="pdf")
            plt.savefig("summary_single_guide_" + str(guide) +
                        "_"+str(uppermm) + "mm" + ".png", format="png")
            # plt.show()

        else:       #mms value is range (Eg0-5)

            def make_spider(row, title, count):
                data_for_df = {'group': ['A'], 'General': [np.argmin(arrayProfiledistance[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0)]}
                for pos, el in enumerate(array_distances):
                    data_for_df[count_files[pos][0]] = [np.argmin(el[:, uppermm])/np.size(arrayexonsMissmatchGlobal[0], 0) ]
                # df = pd.DataFrame({
                #     'group': ['A'],
                #     'Exons': [np.argmin(arrayExonsFdistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)],
                #     'General': [np.argmin(arrayProfiledistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)],
                #     'Introns': [np.argmin(arrayIntronsdistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)],
                #     'Promoters': [np.argmin(arrayPromotersdistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)],
                #     'DNAse': [np.argmin(arrayDNAsedistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)],
                #     'CTCF': [np.argmin(arrayCTCFdistance[:, row])/np.size(arrayCTCFMissmatchGlobal, 0)]
                # })
                df = pd.DataFrame(data_for_df)

                # number of variable
                categories = list(df)[1:]
                N = len(categories)

                # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
                angles = [n / float(N) * 2 * pi for n in range(N)]
                angles += angles[:1]

                # Initialise the spider plot
                ax = plt.subplot(2, (uppermm-lowermm)+1, count, polar=True)
                # table=plt.table(cellText=data,colLabels=string,loc='bottom')

                # If you want the first axis to be on top:
                ax.set_theta_offset(pi / 2)
                ax.set_theta_direction(-1)

                # Draw one axe per variable + add labels labels yet
                multi_xticks = ['General (' + str(arrayprofileMissmatch[row])+')']

                for pos, f in enumerate(arrayexonsMissmatch):
                    multi_xticks.append(count_files[pos][0] + ' (' + str(f[0][row])  + ')')
                #multi_xticks.extend([count_files[pos][0] + ' (' + str(max(f[0][row]))  + ')' for pos, f in enumerate(arrayexonsMissmatch)])
                # plt.xticks(angles[:-1], ['CTCF'+' ('+str(arrayCTCFMissmatch[row])+')', 'DNAse'+' ('+str(arrayDNAseMissmatch[row])+')', 'Exons'+' ('+str(arrayexonsMissmatch[row])+')', 'General' +
                #                          ' ('+str(arrayprofileMissmatch[row])+')', 'Introns'+' ('+str(arrayintronsMissmatch[row])+')', 'Promoters'+' ('+str(arraypromotersMissmatch[row])+')'], color='black')
                plt.xticks(angles[:-1], multi_xticks, color='black')
                for label, rot in zip(ax.get_xticklabels(), angles):
                    if (rot == 0):
                        label.set_horizontalalignment("center")
                    if (rot > 0):
                        label.set_horizontalalignment("left")
                    if (rot > 3):
                        label.set_horizontalalignment("center")
                    if (rot > 4):
                        label.set_horizontalalignment("right")

                # Draw ylabels
                ax.set_rlabel_position(0)
                plt.yticks([0, 0.25, 0.50, 0.75], [
                    "0", "0.25", "0.50", "0.75"], color="black")
                plt.ylim(0, 1)

                # Ind1
                values = df.loc[0].drop('group').values.flatten().tolist()
                values += values[:1]
                ax.plot(angles, values, linewidth=2, linestyle='solid')
                ax.fill(angles, values, alpha=0.4)

                
                texts = ['General: ' + str(max(arrayprofileMissmatchGlobal[:, row]))]
                texts.extend([count_files[pos][0] + ': ' + str(max(f[:,row])) for pos, f in enumerate(arrayexonsMissmatchGlobal)])
                # texts = ["CTCF:"+' '+str(max(arrayCTCFMissmatchGlobal[:, row])), "DNAse:"+' '+str(max(arrayDNAseMissmatchGlobal[:, row])), "Exons:"+' '+str(max(arrayexonsMissmatchGlobal[:, row])), "General:" +
                #          ' '+str(max(arrayprofileMissmatchGlobal[:, row])), "Introns:"+' '+str(max(arrayintronsMissmatchGlobal[:, row])), "Promoters:"+' '+str(max(arraypromotersMissmatchGlobal[:, row]))]
                
                colors = []
                for i in texts:
                    colors.append('white')
                #colors = ["white", "white", "white", "white", "white", "white"]
                patches = [mpatches.Patch(color=colors[i], label="{:s}".format(
                    texts[i])) for i in range(len(texts))]
                plt.legend(handles=patches, loc=(-0.65, 0.86), labelspacing=0.1,
                        ncol=1, handlelength=0, handletextpad=0, title="Max Value")

                # Add a title
                plt.title(title, y=1.2)

            def make_motif(row, count):

                datacount = arrayguidesExtendedProfile[row*7] / \
                    (max(arrayguidesExtendedProfile[row*7]))
                data = np.array(datacount, dtype=float)
                data = np.around(data, decimals=1)
                data.shape = (1, len(datacount))

                string = guide[0:20]

                A = arrayguidesExtendedProfile[row*7+1] / \
                    (max(arrayguidesExtendedProfile[row*7]))
                C = arrayguidesExtendedProfile[row*7+2] / \
                    (max(arrayguidesExtendedProfile[row*7]))
                G = arrayguidesExtendedProfile[row*7+3] / \
                    (max(arrayguidesExtendedProfile[row*7]))
                T = arrayguidesExtendedProfile[row*7+4] / \
                    (max(arrayguidesExtendedProfile[row*7]))

                # the x locations for the groups
                ind = np.arange(0, len(string), 1) + 0.15
                width = 0.7  # the width of the bars: can also be len(x) sequence

                plt.subplot(2, (uppermm-lowermm)+1, count+(uppermm-lowermm)+1)

                p1 = plt.bar(ind, A, width, color='#d62728', align='edge')
                p2 = plt.bar(ind, C, width, bottom=A, align='edge')
                p3 = plt.bar(ind, G, width, bottom=A+C, align='edge')
                p4 = plt.bar(ind, T, width, bottom=C+G+A, align='edge')

                plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A', 'C', 'G', 'T'))

                plt.xlim(0, len(string))
                plt.xticks([])
                table = plt.table(cellText=data, colLabels=string,
                                loc='bottom', cellLoc='center')
                # plt.xticks(ind)
                # table.set_fontsize(20))
                table.auto_set_font_size(False)
                table.set_fontsize(14)
                table.scale(1, 1.4)

            # ------- PART 2: Apply to all individuals
            # initialize the figure
            my_dpi = 96
            plt.figure(figsize=(5000/my_dpi, 5000/my_dpi), dpi=my_dpi)

            # Create a color palette:
            # my_palette = plt.cm.get_cmap("Set2",(uppermm-lowermm)+1)

            count = 1
            # Loop to plot
            for row in range(lowermm, uppermm+1):
                make_spider(row=row, title=str(row) + ' Mismatches', count=count)
                make_motif(row=row, count=count)
                # plt.subplot(3,4,count*2)
                count = count+1

            # plt.tight_layout()
            #plt.subplots_adjust(top=0.85, bottom=0.05, left=0.05,right=0.99, wspace=0.1)

            plt.savefig("summary_multiple_guides_" + str(guide) + "_" +
                        str(lowermm) + "-" + str(uppermm) + "mm" + ".pdf", format="pdf")
            plt.savefig("summary_multiple_guides_" + str(guide) + "_" +
                        str(lowermm) + "-" + str(uppermm) + "mm" + ".png", format="png")
            # plt.show()
