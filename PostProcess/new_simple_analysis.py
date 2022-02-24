#!/usr/bin/env python
from math import trunc
from operator import index
import sys
import json
import os
import pickle
import numpy as np
import pandas as pd
import time
from CRISTA_score import CRISTA_predict_list


# For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv",
                    "TGACYRSWKMDHVBtgacyrswkmdhvb")

iupac_nucleotides = set('RYSWKMBDHVryswkmbdhv')
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
    # if do_scores == False:
    #     # print("NON CALCOLO")
    #     score = -1
    #     return score
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


def iupac_decomposition(split, guide_no_bulge,
                        guide_no_pam, cluster_to_save):
    realTarget = split[2]
    replaceTarget = split[2].replace('-', '')
    refSeq = genomeStr[int(split[4]): int(split[4])+len(replaceTarget)].upper()

    revert = False
    if split[6] == '-':
        revert = True
        replaceTarget = reverse_complement_table(replaceTarget)

    # dict with IUPAC scompositions
    totalDict = dict()
    totalDict[0] = dict()
    totalDict[0][0] = dict()
    if haplotype_check:
        totalDict[1] = dict()
        totalDict[1][0] = dict()
    countIUPAC = 0
    for pos_c, c in enumerate(replaceTarget):
        if c in iupac_code:
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

    if countIUPAC > 0:  # if found valid alternative targets
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
                                # (this should be done only with phased VCF since unphased cannot be verified)
                                if haplotype_check:
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

            for position_t, char_t in enumerate(refSeq_with_bulges[pos_beg: pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != '-':
                        refSeq_with_bulges[pos_beg +
                                           position_t] = char_t.lower()
            # ref sequence with bulges
            refSeq_with_bulges = ''.join(refSeq_with_bulges)

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
                        for position_t, char_t in enumerate(target_to_list[pos_beg: pos_end]):
                            if char_t.upper() != guide_no_pam[position_t]:
                                mm_new_t += 1
                                tmp_pos_mms = position_t
                                if guide_no_pam[position_t] != '-':
                                    target_to_list[pos_beg +
                                                   position_t] = char_t.lower()

                        # pam respect input PAM after IUPAC resolution
                        pam_ok = True
                        for pam_chr_pos, pam_chr in enumerate(target_to_list[pam_begin: pam_end]):
                            if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                                pam_ok = False

                        target_pam_ref = refSeq_with_bulges[pam_begin: pam_end]
                        found_creation = False
                        for pos_pam, pam_char in enumerate(target_pam_ref):
                            # ref char not in set of general pam char
                            if not iupac_code_set[pam[pos_pam]] & iupac_code_set[pam_char]:
                                found_creation = True
                        # value of mm and bulges is over allowed threshold, discard target
                        if mm_new_t - int(split[8]) > allowed_mms:
                            continue
                        elif pam_ok:
                            final_line[2] = ''.join(target_to_list)
                            final_line[7] = str(mm_new_t - int(final_line[8]))
                            # total differences between targets and guide (mismatches + bulges)
                            final_line[9] = str(mm_new_t)
                            if found_creation:
                                final_line[10] = ''.join(
                                    target_to_list[pam_begin: pam_end])
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
                            # report ref DNA seq of target
                            final_line.append(refSeq_with_bulges)
                            # number to activate ref score calculation (active if target is alternative)
                            final_line.append(33)
                            # position of tmp_mms (removed later after processing)
                            final_line.append(tmp_pos_mms)
                            # append processed target to cluster to save
                            cluster_to_save.append(final_line)


def preprocess_CFD_score(target):
    # preprocess target then calculate CFD score
    if do_scores:
        if target[0] == 'DNA':
            cfd_score = calc_cfd(target[1][int(target[bulge_pos]):], target[2].upper()[int(
                target[bulge_pos]):-3], target[2].upper()[-2:], mm_scores, pam_scores, do_scores)
            # append to target the CFD score of the aligned sequence (alt or ref)
            target.append("{:.3f}".format(cfd_score))
            # -3 position is a placeholder for ref score
            if target[-3] == 55:  # if 55 sequence is ref so no score have to be calculated
                target[-3] = "{:.3f}".format(cfd_score)
            if target[-3] == 33:  # if 33 sequence is alt so ref score must be calculated
                cfd_ref_score = calc_cfd(target[1][int(target[bulge_pos]):], target[-4].upper()[int(
                    target[bulge_pos]):-3], target[-4].upper()[-2:], mm_scores, pam_scores, do_scores)
                target[-3] = "{:.3f}".format(cfd_ref_score)
        else:
            cfd_score = calc_cfd(target[1], target[2].upper()[
                :-3], target[2].upper()[-2:], mm_scores, pam_scores, do_scores)
            target.append("{:.3f}".format(cfd_score))
            if target[-3] == 55:
                target[-3] = "{:.3f}".format(cfd_score)
            if target[-3] == 33:
                cfd_ref_score = calc_cfd(target[1], target[-4].upper()[
                    :-3], target[-4].upper()[-2:], mm_scores, pam_scores, do_scores)
                target[-3] = "{:.3f}".format(cfd_ref_score)
    else:
        # no score calculated, append -1 in CFD score and in position -3 insert -1 value (-1 means no score calculated)
        cfd_score = -1
        target.append("{:.3f}".format(cfd_score))
        target[-3] = "{:.3f}".format(cfd_score)
    # print(target)
    # print('cfd', cfd_score)
    return target


def preprocess_CRISTA_score(cluster_targets):
    # list with scored targets
    cluster_scored = list()
    index_to_null = list()

    # skip scoring for CRISTA, remove to activate scoring
    # do_scores = False

    if do_scores:
        pass
    else:
        for target in cluster_targets:
            target_CRISTA = target.copy()
            crista_score = -1  # null score
            target_CRISTA[-2] = "{:.3f}".format(crista_score)
            target_CRISTA.append("{:.3f}".format(crista_score))
            cluster_scored.append(target_CRISTA)
        return cluster_scored

    # preprocess target then calculate CRISTA score
    sgRNA_non_aligned_list = list()
    DNA_aligned_list = list()
    DNAseq_from_genome_list = list()
    # process all found targets
    for index, target in enumerate(cluster_targets):
        # list with non-aligned sgRNA
        sgRNA_non_aligned_list.append(
            str(target[1])[:len(str(target[1]))-3]+'NGG')
        # list with aligned DNA
        DNA_aligned_list.append(str(target[2]))
        # first 5 nucleotide to add to protospacer
        pre_protospacer_DNA = genomeStr[int(
            target[4])-5:int(target[4])].upper()
        # protospacer taken directly from the aligned target
        protospacerDNA = str(target[2]).replace('-', '')
        if target[6] == '-':
            protospacerDNA = reverse_complement_table(protospacerDNA)
        # last 5 nucleotides to add to protospacer
        post_protospacer_DNA = genomeStr[int(
            target[4])+len(target[1]):int(target[4])+len(target[1])+5].upper()

        # DNA seq extracted from genome and append to aligned DNA seq from CRISPRme
        complete_DNA_seq = str(pre_protospacer_DNA) + \
            protospacerDNA+str(post_protospacer_DNA)

        for elem in iupac_nucleotides:
            if elem in complete_DNA_seq:
                complete_DNA_seq = complete_DNA_seq.replace(elem, '')

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq/2)-14:int(len_DNA_seq/2)]
        second_half = complete_DNA_seq[int(
            len_DNA_seq/2):int(len_DNA_seq/2)+15]
        complete_DNA_seq = first_half+second_half
        if target[6] == '-':
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

        # if 'N' is present in the reference DNA seq, we must use a fake DNA seq to complete the aligned
        # that will be discarded after
        if 'N' in complete_DNA_seq or 'n' in complete_DNA_seq or 'N' in DNA_aligned_list[-1] or 'n' in DNA_aligned_list[-1]:
            complete_DNA_seq = 'A'*29
            DNA_aligned_list[-1] = 'A'*len(str(target[2]))
            index_to_null.append(index)

        # append sequence to DNA list
        DNAseq_from_genome_list.append(complete_DNA_seq)

    # calculate scores for alt sequence
    crista_score_list_alt = list()
    if do_scores:
        crista_score_list_alt = CRISTA_predict_list(
            sgRNA_non_aligned_list, DNA_aligned_list, DNAseq_from_genome_list)

    # preprocess target then calculate CRISTA score
    sgRNA_non_aligned_list = list()
    DNA_aligned_list = list()
    DNAseq_from_genome_list = list()
    # process all ref sequences in targets
    for index, target in enumerate(cluster_targets):
        # list with non-aligned sgRNA
        sgRNA_non_aligned_list.append(
            str(target[1])[:len(str(target[1]))-3]+'NGG')
        # list with aligned DNA
        if 'n' not in target[-3]:
            DNA_aligned_list.append(str(target[-3]))
        else:
            DNA_aligned_list.append(str(target[2]))
        # first 5 nucleotide to add to protospacer
        pre_protospacer_DNA = genomeStr[int(
            target[4])-5:int(target[4])]
        # protospacer taken directly from the ref genome
        protospacerDNA = genomeStr[int(target[4]):int(
            target[4])+len(target[1])]
        # last 5 nucleotides to add to protospacer
        post_protospacer_DNA = genomeStr[int(
            target[4])+len(target[1]):int(target[4])+len(target[1])+5]

        # DNA seq extracted from genome and append to aligned DNA seq from CRISPRme
        complete_DNA_seq = str(pre_protospacer_DNA) + \
            protospacerDNA+str(post_protospacer_DNA)

        for elem in iupac_nucleotides:
            if elem in complete_DNA_seq:
                complete_DNA_seq = complete_DNA_seq.replace(elem, '')

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq/2)-14:int(len_DNA_seq/2)]
        second_half = complete_DNA_seq[int(
            len_DNA_seq/2):int(len_DNA_seq/2)+15]
        complete_DNA_seq = first_half+second_half
        if target[6] == '-':
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

        # if 'N' is present in the reference DNA seq, we must use a fake DNA seq to complete the aligned
        # that will be discarded after
        if 'N' in complete_DNA_seq or 'n' in complete_DNA_seq or 'N' in DNA_aligned_list[-1] or 'n' in DNA_aligned_list[-1]:
            complete_DNA_seq = 'A'*29
            DNA_aligned_list[-1] = 'A'*len(str(target[2]))
            index_to_null.append(index)

        # append sequence to DNA list
        DNAseq_from_genome_list.append(complete_DNA_seq)

    # calculate score
    crista_score_list_ref = list()
    if do_scores:
        crista_score_list_ref = CRISTA_predict_list(
            sgRNA_non_aligned_list, DNA_aligned_list, DNAseq_from_genome_list)

    for index, target in enumerate(cluster_targets):
        target_CRISTA = target.copy()
        # if any of the scored target is not valid, due to Ns in the sequence, return a -1 score
        if index in index_to_null:
            crista_score = -1  # null score
            target_CRISTA[-2] = "{:.3f}".format(crista_score)
            target_CRISTA.append("{:.3f}".format(crista_score))
        else:
            # else report the correct score
            if target_CRISTA[-2] == 55:  # reference target have duplicate score
                target_CRISTA[-2] = "{:.3f}".format(
                    crista_score_list_alt[index])
                target_CRISTA.append("{:.3f}".format(
                    crista_score_list_alt[index]))
            if target_CRISTA[-2] == 33:  # alternative target scoring
                target_CRISTA[-2] = "{:.3f}".format(
                    crista_score_list_ref[index])
                target_CRISTA.append("{:.3f}".format(
                    crista_score_list_alt[index]))
        # append to final score cluster
        cluster_scored.append(target_CRISTA)

    return cluster_scored


def calculate_scores(cluster_to_save):
    # function to calculate score for each input target
    # input is target line splitted in list format
    # list of functions to calculate specific score (to add a score, simply add your function to this call and update the clusters list in return)
    cluster_with_CFD_score = list()
    cluster_with_CRISTA_score = list()

    for target in cluster_to_save:  # calculate CFD score for each target
        target_CFD = target.copy()
        cluster_with_CFD_score.append(preprocess_CFD_score(target_CFD))

    # process score for each target in cluster, at the same time to improve execution time
    cluster_with_CRISTA_score = preprocess_CRISTA_score(cluster_to_save)

    # analyze CFD scored targets, returning for each guide,chr,cluster_pos the highest scoring target (or multiple in case of equal)
    df_CFD = pd.DataFrame(cluster_with_CFD_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
                                                           'Position', 'Cluster_Position', 'Direction', 'Mismatches',
                                                           'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
                                                           'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
                                                           'Seq_in_cluster', 'CFD_ref'])
    # group by over real_guide,chr,cluster_pos to avoid mixing targets
    # select lowest count of mm+bul
    idx_fewest_mm_bul = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
        'Total'].transform(min) == df_CFD['Total']
    df_CFD_fewest = df_CFD[idx_fewest_mm_bul]
    df_CFD_fewest.drop_duplicates(
        ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # select highest score
    idx_max_score = df_CFD.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
        'CFD'].transform(max) == df_CFD['CFD']
    df_CFD_best_score = df_CFD[idx_max_score]
    # remove duplicate rows (possible due to haplotypes and variants)
    df_CFD_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
                                       'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    frames = [df_CFD_fewest, df_CFD_best_score]
    df_CFD = pd.concat(frames)
    cluster_with_CFD_score = df_CFD.values.tolist()

    df_CRISTA = pd.DataFrame(cluster_with_CRISTA_score, columns=['Bulge_type', 'crRNA', 'DNA', 'Chromosome',
                                                                 'Position', 'Cluster_Position', 'Direction', 'Mismatches',
                                                                 'Bulge_Size', 'Total', 'PAM_gen', 'Var_uniq', 'Samples', 'Annotation_Type',
                                                                 'Real_Guide', 'rsID', 'AF', 'SNP', 'Reference_target', 'CFD',
                                                                 'Seq_in_cluster', 'CFD_ref'])
    # select lowest count of mm+bul
    idx_fewest_mm_bul = df_CRISTA.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
        'Total'].transform(min) == df_CRISTA['Total']
    df_CRISTA_fewest = df_CRISTA[idx_fewest_mm_bul]
    df_CRISTA_fewest.drop_duplicates(
        ['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples', 'Mismatches', 'Bulge_Size'], inplace=True)
    # select highest score
    idx_max_score = df_CRISTA.groupby(['Real_Guide', 'Chromosome', 'Cluster_Position', 'SNP', 'Samples'])[
        'CFD'].transform(max) == df_CRISTA['CFD']
    df_CRISTA_best_score = df_CRISTA[idx_max_score]
    # remove duplicate rows (possible due to haplotypes and variants)
    df_CRISTA_best_score.drop_duplicates(['Real_Guide', 'Chromosome',
                                          'Cluster_Position', 'SNP', 'Samples', 'CFD'], inplace=True)
    frames = [df_CRISTA_fewest, df_CRISTA_best_score]
    df_CRISTA = pd.concat(frames)
    cluster_with_CRISTA_score = df_CRISTA.values.tolist()

    return [cluster_with_CFD_score, cluster_with_CRISTA_score]


# INPUT AND SETTINGS
# fasta of the reference chromosome
inFasta = open(sys.argv[1], 'r')
current_chr = inFasta.readline().strip().replace('>', '')  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = ''.join(genomeStr).upper()
# string of the whole chromosome on single line
genomeStr = genomeStr.replace('\n', '')
# targets clusterized by chr and ordered by position
inTarget = open(sys.argv[3], 'r')
# text file with PAM sequence and length
inPAMfile = open(sys.argv[4], 'r')
# outfile path
outputFile = sys.argv[5]
# max allowed mismatches in search (to validate ref targets in alternative case)
allowed_mms = int(sys.argv[6])
# column of bulges count
bulge_pos = 8
# header to insert into final file
header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference'
# cfd graphs pre-processing (deprecated)
cfd_for_graph = {'ref': [0]*101, 'var': [0]*101}

# OUT BEST FILES FOR EACH SCORING SYSTEM

# file with best CFD targets
cfd_best = open(outputFile + '.bestCFD.txt', 'w+')
cfd_best.write(header + '\tCFD\n')  # Write header

# file with best mm+bul targets
mmblg_best = open(outputFile + '.bestmmblg.txt', 'w+')
mmblg_best.write(header + '\tCFD\n')  # Write header

# file with best CRISTA targets
crista_best = open(outputFile + '.bestCRISTA.txt', 'w+')
crista_best.write(header + '\tCFD\n')  # Write header

# check if dictionaries has haplotypes
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
    print('Haplotype processing', haplotype_check)
except:
    print("No dict found for", current_chr)

# check PAM position and relative coordinates on targets
pam_at_beginning = False
line = inPAMfile.read().strip()
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

# start time counter
global_start = time.time()

# open mm and pam scores matrices for CFD
mm_scores, pam_scores = get_mm_pam_scores()

# if conditions, execute score (guidelen==20,pamlen==3,pam_at_beginning==FALSE)
do_scores = True
if len_pam != 3 or guide_len != 20 or pam_at_beginning:
    # sys.stderr.write('CFD SCORE IS NOT CALCULATED WITH GUIDES LENGTH != 20 OR PAM LENGTH !=3 OR UPSTREAM PAM')
    do_scores = False

# START TARGET PROCESSING

# keep track of current analyzed cluster (necessary to check if the cluster is terminated)
# current_guide_chr_pos_direction = ''
# count_cluster_dimension = 0

# skip header
inTarget.readline()

# list with clusterized targets in list format (contains ref seq and all other alternative targets)
cluster_to_save = list()
# read lines from target file
for line in inTarget:
    # split target into list
    split = line.strip().split('\t')
    # sgRNA sequence (with bulges and PAM)
    guide = split[1]
    # found target on DNA (with bulges, mismatches and PAM)
    target = split[2]
    guide_no_bulge = split[1].replace('-', '')
    guide_no_pam = guide[pos_beg:pos_end]

    # check if targets cointains IUPAC nucleotide
    if any((c in iupac_nucleotides) for c in target):
        iupac_decomposition(split, guide_no_bulge,
                            guide_no_pam, cluster_to_save)
    else:
        # process_iupac = False
        # append to respect file format for post analysis
        # null ref sequence
        split.append('n')
        # specific value to represent a ref target to avoid recount score
        split.append(55)
        # count of mm_bul for ref sequence in case of alternative target
        split.append(0)
        cluster_to_save.append(split)

    if len(cluster_to_save) >= 100000:
        # after reading 100k lines from file and creating the cluster, start processing it
        clusters_with_scores = calculate_scores(cluster_to_save)

        for count, cluster in enumerate(clusters_with_scores):
            for target in cluster:
                if count == 0:  # CFD target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CFD targets
                    cfd_best.write(
                        '\t'.join(target)+'\t'+str(0)+'\n')
                    # save mm-bul targets
                    mmblg_best.write(
                        '\t'.join(target)+'\t'+str(0)+'\n')
                if count == 1:  # CRISTA target
                    # remove count of tmp_mms
                    target.pop(-2)
                    # save CRISTA targets
                    crista_best.write('\t'.join(target)+'\t' +
                                      str(0)+'\n')
        cluster_to_save = list()


if len(cluster_to_save):
    pass
else:
    # if cluster to save is empty, skip processing
    # close all files
    cfd_best.close()
    mmblg_best.close()
    crista_best.close()
    # rewrite header file
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCFD.txt')
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestmmblg.txt')
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCRISTA.txt')
    # cfd dataframe write
    cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
    cfd_dataframe.to_csv(outputFile + '.CFDGraph.txt', sep='\t', index=False)
    # print complete and exit with no error
    print('ANALYSIS COMPLETE IN', time.time() - global_start)
    exit(0)

# process cluster of targets if less then 1mln rows total
clusters_with_scores = calculate_scores(cluster_to_save)

for count, cluster in enumerate(clusters_with_scores):
    for target in cluster:
        if count == 0:  # CFD target
            # remove count of tmp_mms
            target.pop(-2)
            # save CFD targets
            cfd_best.write(
                '\t'.join(target)+'\t'+str(0)+'\n')
            # save mm-bul targets
            mmblg_best.write(
                '\t'.join(target)+'\t'+str(0)+'\n')
        if count == 1:  # CRISTA target
            # remove count of tmp_mms
            target.pop(-2)
            # save CRISTA targets
            crista_best.write('\t'.join(target)+'\t' +
                              str(0)+'\n')

cfd_best.close()
mmblg_best.close()
crista_best.close()

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCFD.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestmmblg.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCRISTA.txt')


cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
cfd_dataframe.to_csv(outputFile + '.CFDGraph.txt', sep='\t', index=False)

print('ANALYSIS COMPLETE IN', time.time() - global_start)
