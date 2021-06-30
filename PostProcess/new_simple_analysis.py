#!/usr/bin/env python
from math import trunc
import sys
import json
import os
import pickle
import numpy as np
import pandas as pd
import time

inFasta = open(sys.argv[1], 'r')  # lettura fasta del chr
current_chr = inFasta.readline().strip().replace('>', '')  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = ''.join(genomeStr).upper()
genomeStr = genomeStr.replace('\n', '')
inTarget = open(sys.argv[3], 'r')
haplotype_check = False

try:
    inDict = open(sys.argv[2], 'r')
    mydict = json.load(inDict)
    for entry in mydict:
        if '|' in mydict[entry]:
            haplotype_check = True
            break
        elif '/' in mydict[entry]:
            break
    print('haplo check', haplotype_check)
except:
    print("No dict found for", current_chr)

iupac_code_set = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "r": {"A", "G"},
    "y": {"C", "T"},
    "s": {"G", "C"},
    "w": {"A", "T"},
    "k": {"G", "T"},
    "m": {"A", "C"},
    "b": {"C", "G", "T"},
    "d": {"A", "G", "T"},
    "h": {"A", "C", "T"},
    "v": {"A", "C", "G"},
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
    "a": {"a"},
    "t": {"t"},
    "c": {"c"},
    "g": {"g"},
    'N': {'A', 'T', 'G', 'C'}
}


iupac_code = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "r": ("A", "G"),
    "y": ("C", "T"),
    "s": ("G", "C"),
    "w": ("A", "T"),
    "k": ("G", "T"),
    "m": ("A", "C"),
    "b": ("C", "G", "T"),
    "d": ("A", "G", "T"),
    "h": ("A", "C", "T"),
    "v": ("A", "C", "G"),
    'N': ('A', 'T', 'C', 'G')
}

# For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv",
                    "TGACYRSWKMDHVBtgacyrswkmdhvb")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


class reversor:
    '''
    Nel caso debba ordinare più campi però con reverse diversi, eg uno True e l'altro False, posso usare questa classe nella chiave per 
    simulare il contrario del reverse applicato
    '''

    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return other.obj == self.obj

    def __lt__(self, other):
        return other.obj < self.obj


def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores, do_scores):
    if do_scores == False:
        # print("NON CALCOLO")
        score = -1
        return score
    score = 1
    sg = sg.replace('T', 'U')
    guide_seq = guide_seq.replace('T', 'U')
    s_list = list(sg)
    guide_seq_list = list(guide_seq)

    for i, sl in enumerate(s_list):
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            try:  # Catch exception if IUPAC character
                key = 'r' + guide_seq_list[i] + \
                    ':d' + revcom(sl) + ',' + str(i + 1)
            except Exception as e:
                score = 0
                break
            try:
                score *= mm_scores[key]
            except Exception as e:  # If '-' is in first position, i do not have the score for that position
                pass

    score *= pam_scores[pam]
    return score


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', '-': '-'}
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None  # If some IUPAC were not translated
    return ''.join(letters)


def get_mm_pam_scores():
    # print(os.path.dirname(os.path.realpath(__file__)))
    try:
        mm_scores = pickle.load(open(os.path.dirname(
            os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
        pam_scores = pickle.load(open(os.path.dirname(
            os.path.realpath(__file__)) + '/PAM_scores.pkl', 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception(
            "Could not find file with mismatch scores or PAM scores")


def retrieveFromDict(chr_pos):
    try:
        entry = mydict[current_chr+','+str(chr_pos+1)]
    except:
        snp_list = []
        sample_list = []
        AF_list = []
        rsID_list = []
        snp_info_list = []
        sample_list.append([])  # no samples
        snp_list.append('C')  # fake snp
        rsID_list.append('.')  # no rsid
        AF_list.append('0')  # fake AF
        snp_info_list.append(
            current_chr+'_'+str(chr_pos+1)+'_'+'C'+'_'+'G')  # fake snp info list
        return snp_list, sample_list, rsID_list, AF_list, snp_info_list
    multi_entry = entry.split('$')
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for entry in multi_entry:
        split_entry = entry.split(';')
        samples = split_entry[0].strip().split(',')
        if samples[0] == '':
            samples = []
        sample_list.append(samples)
        snp_list.append(split_entry[1].strip().split(',')[1])
        rsID_list.append(split_entry[2].strip())
        AF_list.append(split_entry[3].strip())
        snp_info_list.append(
            current_chr+'_'+str(chr_pos+1)+'_'+split_entry[1].split(',')[0]+'_'+split_entry[1].split(',')[1])
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


pam_at_beginning = False
with open(sys.argv[4], 'r') as pam:
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

global_start = time.time()

cluster_to_save = []
target_no_bulges_in_cluster = []
bulge_pos = 8

header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference'
outputFile = sys.argv[5]

cfd_best = open(outputFile + '.bestCFD.txt', 'w+')
cfd_best.write(header + '\tCFD\n')  # Write header

mmblg_best = open(outputFile + '.bestmmblg.txt', 'w+')
mmblg_best.write(header + '\tCFD\n')  # Write header

mm_scores, pam_scores = get_mm_pam_scores()

do_scores = True
if len_pam != 3 or guide_len != 20 or pam_at_beginning:
    # sys.stderr.write('CFD SCORE IS NOT CALCULATED WITH GUIDES LENGTH != 20 OR PAM LENGTH !=3 OR UPSTREAM PAM')
    do_scores = False

allowed_mms = int(sys.argv[6])
cfd_for_graph = {'ref': [0]*101, 'var': [0]*101}
current_guide_chr_pos_direction = ''
inTarget.readline()
for line in inTarget:
    split = line.strip().split('\t')
    guide = split[1]
    guide_no_bulge = split[1].replace('-', '')
    guide_no_pam = guide[pos_beg:pos_end]
    # SAVE TO FILE
    if (guide_no_bulge + split[3] + split[5] + split[6]) != current_guide_chr_pos_direction:
        current_guide_chr_pos_direction = guide_no_bulge + \
            split[3] + split[5] + split[6]  # update next cluster key

        for t in cluster_to_save:

            if t[0] == 'DNA':
                cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(
                    t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
                # t.append(str(cfd_score))
                t.append("{:.3f}".format(cfd_score))
            else:
                cfd_score = calc_cfd(t[1], t[2].upper()[
                                     :-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
                # t.append(str(cfd_score))
                t.append("{:.3f}".format(cfd_score))

        cluster_to_save.sort(key=lambda x: (
            float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse=True)

        cluster_to_save_mmbl = cluster_to_save.copy()
        cluster_to_save_mmbl.sort(key=lambda x: (int(x[8]), int(x[7])))

        keys_seen = []
        saved = False

        if cluster_to_save:
            c = cluster_to_save[0]
            c.pop(-2)
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]

            if c[12] == 'n':
                cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
            else:
                cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1

            cfd_best.write('\t'.join(c))
            cluster_to_save.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True

        list_for_alt = []
        for c in cluster_to_save:
            c.pop(-2)
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                if c[12] == 'n':
                    cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                else:
                    cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                list_for_alt.append('\t'.join(c))
        if saved:
            cfd_best.write('\t' + str(len(list_for_alt)) + '\n')
        for ele in list_for_alt:
            cfd_best.write(ele+'\t'+str(len(list_for_alt)) + '\n')

        keys_seen = []
        saved = False
        if cluster_to_save_mmbl:
            c = cluster_to_save_mmbl[0]
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]

            mmblg_best.write('\t'.join(c))
            cluster_to_save_mmbl.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True

        list_for_alt = []
        for c in cluster_to_save_mmbl:
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                list_for_alt.append('\t'.join(c))

        if saved:
            mmblg_best.write('\t' + str(len(list_for_alt)) + '\n')
        for ele in list_for_alt:
            mmblg_best.write(ele+'\t'+str(len(list_for_alt)) + '\n')

        cluster_to_save = []

    realTarget = split[2]
    replaceTarget = split[2].replace('-', '')
    refSeq = genomeStr[int(split[4]):int(split[4])+len(replaceTarget)].upper()
    replaceTargetsDict = dict()

    revert = False
    if split[6] == '-':
        revert = True
        replaceTarget = reverse_complement_table(replaceTarget)

    # print('ref', refSeq)
    # print('var original', replaceTarget)

    totalDict = dict()
    totalDict[0] = dict()
    totalDict[0][0] = dict()
    if haplotype_check:
        totalDict[1] = dict()
        totalDict[1][0] = dict()
    countIUPAC = 0
    for pos_c, c in enumerate(replaceTarget):
        if c in iupac_code:
            # print(c)
            countIUPAC += 1
            snpToReplace, sampleSet, rsID, AF_var, snpInfo = retrieveFromDict(
                pos_c+int(split[4]))
            for i, elem in enumerate(snpToReplace):
                listReplaceTarget = list(refSeq)
                listReplaceTarget[pos_c] = elem
                listInfo = [[rsID[i], AF_var[i], snpInfo[i]]]
                if haplotype_check:
                    haploSamples = {0: [], 1: []}
                    for count, sample in enumerate(sampleSet[i]):
                        sampleInfo = sample.split(':')
                        for haplo in totalDict:
                            if sampleInfo[1].split('|')[haplo] != '0':
                                haploSamples[haplo].append(sampleInfo[0])
                    totalDict[0][0][(pos_c, elem)] = [
                        listReplaceTarget, set(haploSamples[0]), listInfo]
                    totalDict[1][0][(pos_c, elem)] = [
                        listReplaceTarget, set(haploSamples[1]), listInfo]
                else:
                    sampleList = list()
                    for count, sample in enumerate(sampleSet[i]):
                        sampleList.append(sample.split(':')[0])
                    totalDict[0][0][(pos_c, elem)] = [
                        listReplaceTarget, set(sampleList), listInfo]

    if countIUPAC > 0:
        if revert:
            refSeq = reverse_complement_table(refSeq)
        for count in totalDict:
            createdNewLayer = True
            for size in range(countIUPAC):  # the time of the universe
                if createdNewLayer:
                    createdNewLayer = False
                else:
                    break
                totalDict[count][size+1] = dict()
                # for each snp in target (fixpoint)
                for key in totalDict[count][size]:
                    # for each other snp in target (> fixpoint)
                    for newkey in totalDict[count][0]:
                        if newkey[-2] > key[-2]:
                            resultSet = totalDict[count][size][key][1].intersection(
                                totalDict[count][0][newkey][1])  # extract intersection of sample to generate possible multisnp target
                            if len(resultSet) > 0:  # if set is not null
                                createdNewLayer = True
                                # add new snp to preceding target seq with snp
                                replaceTarget1 = totalDict[count][0][newkey][0].copy(
                                )
                                replaceTarget2 = totalDict[count][size][key][0].copy(
                                )
                                replaceTarget2[newkey[0]
                                               ] = replaceTarget1[newkey[0]]
                                listInfo2 = totalDict[count][size][key][2].copy(
                                )
                                listInfo2.extend(
                                    totalDict[count][0][newkey][2])
                                # add to next level the modified seq and set of samples and info of snp
                                combinedKey = key + newkey
                                totalDict[count][size+1][combinedKey] = [
                                    replaceTarget2, resultSet, listInfo2]
                                # remove the new generated sample set from all lower levels
                                totalDict[count][size][key][1] = totalDict[count][size][key][1] - \
                                    totalDict[count][size +
                                                     1][combinedKey][1]
                                totalDict[count][0][newkey][1] = totalDict[count][0][newkey][1] - \
                                    totalDict[count][size +
                                                     1][combinedKey][1]

            refSeq_with_bulges = list(refSeq)
            for pos, char in enumerate(realTarget):
                if char == '-':
                    refSeq_with_bulges.insert(pos, '-')

            for position_t, char_t in enumerate(refSeq_with_bulges[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != '-':
                        refSeq_with_bulges[pos_beg +
                                           position_t] = char_t.lower()

            refSeq_with_bulges = ''.join(refSeq_with_bulges)
            if split[0] == 'DNA':
                cfd_score = calc_cfd(split[1][int(split[bulge_pos]):], refSeq_with_bulges.upper()[int(
                    split[bulge_pos]):-3], refSeq_with_bulges.upper()[-2:], mm_scores, pam_scores, do_scores)
                cfd_ref_seq = "{:.3f}".format(cfd_score)  # str(cfd_score)
            else:
                cfd_score = calc_cfd(split[1], refSeq_with_bulges.upper()[
                    :-3], refSeq_with_bulges.upper()[-2:], mm_scores, pam_scores, do_scores)
                cfd_ref_seq = "{:.3f}".format(cfd_score)  # str(cfd_score)

            for level in totalDict[count]:
                for key in totalDict[count][level]:
                    if len(totalDict[count][level][key][1]) > 0:
                        if revert:
                            totalDict[count][level][key][0] = reverse_complement_table(
                                ''.join(totalDict[count][level][key][0]))
                        else:
                            totalDict[count][level][key][0] = ''.join(
                                totalDict[count][level][key][0])

                        final_line = split.copy()

                        target_to_list = list(totalDict[count][level][key][0])
                        for pos, char in enumerate(realTarget):
                            if char == '-':
                                target_to_list.insert(pos, '-')

                        mm_new_t = 0
                        tmp_pos_mms = 0
                        for position_t, char_t in enumerate(target_to_list[pos_beg:pos_end]):
                            if char_t.upper() != guide_no_pam[position_t]:
                                mm_new_t += 1
                                tmp_pos_mms = position_t
                                if guide_no_pam[position_t] != '-':
                                    target_to_list[pos_beg +
                                                   position_t] = char_t.lower()

                        pam_ok = True
                        for pam_chr_pos, pam_chr in enumerate(target_to_list[pam_begin:pam_end]):
                            if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                                pam_ok = False

                        target_pam_ref = refSeq_with_bulges[pam_begin:pam_end]
                        found_creation = False
                        for pos_pam, pam_char in enumerate(target_pam_ref):
                            # ref char not in set of general pam char
                            if not iupac_code_set[pam[pos_pam]] & iupac_code_set[pam_char]:
                                found_creation = True

                        if mm_new_t - int(split[8]) > allowed_mms:
                            continue
                        elif pam_ok:
                            final_line[2] = ''.join(target_to_list)
                            final_line[7] = str(mm_new_t - int(final_line[8]))
                            # total differences between targets and guide (mismatches + bulges)
                            final_line[9] = str(mm_new_t)
                            if found_creation:
                                final_line[10] = ''.join(
                                    target_to_list[pam_begin:pam_end])
                            final_line[12] = ','.join(
                                totalDict[count][level][key][1])
                            tmp_matrix = np.array(
                                totalDict[count][level][key][2])
                            if tmp_matrix.shape[0] > 1:
                                final_line[15] = ','.join(tmp_matrix[:, 0])
                                final_line[16] = ','.join(tmp_matrix[:, 1])
                                final_line[17] = ','.join(tmp_matrix[:, 2])
                            else:
                                final_line[15] = str(tmp_matrix[0][0])
                                final_line[16] = str(tmp_matrix[0][1])
                                final_line[17] = str(tmp_matrix[0][2])

                            final_line.append(refSeq_with_bulges)
                            final_line.append(cfd_ref_seq)
                            final_line.append(tmp_pos_mms)
                            cluster_to_save.append(final_line)
    else:
        final_line = split.copy()

        tmp_pos_mms = 0
        for position_t, char_t in enumerate(final_line[2][pos_beg:pos_end]):
            if char_t.upper() != guide_no_pam[position_t]:
                tmp_pos_mms = position_t
        if split[0] == 'DNA':
            cfd_score = calc_cfd(split[1][int(split[bulge_pos]):], split[2].upper()[int(
                split[bulge_pos]):-3], split[2].upper()[-2:], mm_scores, pam_scores, do_scores)
            cfd_ref_seq = "{:.3f}".format(cfd_score)  # str(cfd_score)
        else:
            cfd_score = calc_cfd(split[1], split[2].upper()[
                :-3], split[2].upper()[-2:], mm_scores, pam_scores, do_scores)
            cfd_ref_seq = "{:.3f}".format(cfd_score)  # str(cfd_score)

        final_line.append("n")
        final_line.append(cfd_ref_seq)
        final_line.append(tmp_pos_mms)
        cluster_to_save.append(final_line)


for t in cluster_to_save:

    if t[0] == 'DNA':
        cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(
            t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        # t.append(str(cfd_score))
        t.append("{:.3f}".format(cfd_score))
    else:
        cfd_score = calc_cfd(t[1], t[2].upper()[
            :-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        # t.append(str(cfd_score))
        t.append("{:.3f}".format(cfd_score))

cluster_to_save.sort(key=lambda x: (
    float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse=True)

cluster_to_save_mmbl = cluster_to_save.copy()
cluster_to_save_mmbl.sort(key=lambda x: (int(x[8]), int(x[7])))

keys_seen = []
saved = False

if cluster_to_save:
    c = cluster_to_save[0]
    c.pop(-2)
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]

    if c[12] == 'n':
        cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
    else:
        cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1

    cfd_best.write('\t'.join(c))
    cluster_to_save.pop(0)
    keys_seen.append(cfd_clus_key)
    saved = True

list_for_alt = []
for c in cluster_to_save:
    c.pop(-2)
    new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]
    if new_cfd_clus_key not in keys_seen:
        keys_seen.append(new_cfd_clus_key)
        if c[12] == 'n':
            cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
        else:
            cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
        list_for_alt.append('\t'.join(c))
if saved:
    cfd_best.write('\t' + str(len(list_for_alt)) + '\n')
for ele in list_for_alt:
    cfd_best.write(ele+'\t'+str(len(list_for_alt)) + '\n')

keys_seen = []
saved = False
if cluster_to_save_mmbl:
    c = cluster_to_save_mmbl[0]
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]

    mmblg_best.write('\t'.join(c))
    cluster_to_save_mmbl.pop(0)
    keys_seen.append(cfd_clus_key)
    saved = True

list_for_alt = []
for c in cluster_to_save_mmbl:
    new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[17]
    if new_cfd_clus_key not in keys_seen:
        keys_seen.append(new_cfd_clus_key)
        list_for_alt.append('\t'.join(c))

if saved:
    mmblg_best.write('\t' + str(len(list_for_alt)) + '\n')
for ele in list_for_alt:
    mmblg_best.write(ele+'\t'+str(len(list_for_alt)) + '\n')

cluster_to_save = []


cfd_best.close()
mmblg_best.close()

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCFD.txt')

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestmmblg.txt')


cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
cfd_dataframe.to_csv(outputFile + '.CFDGraph.txt', sep='\t', index=False)

print('ANALYSIS COMPLETE IN', time.time() - global_start)
