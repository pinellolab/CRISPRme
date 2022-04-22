#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 19:55:54 2020

@author: francesco
"""

'''
Calcolo semi bruteforce dei sample.
Per ogni target, se è un nuovo cluster, salva il cluster precedente, poi calcola le possibili scomposizioni esistenti, e passa al prossimo target.
Se stesso cluster: se il target è di una categoria già analizzata (X0, DNA1, DNA2 ... RNA1,RNA2 ... dove il numero indica i bulges presenti), usa
i dati della scomposizione già effettuata per evitare il ricalcolo della scomposizione.
Se invece la categoria non è stata già analizzata, fa una scomposizione completa (come il top1), e salva i dati per i prossimi target

#TODO parallelizzare l'analisi
#NOTE da verificare se è compatibile com pam all'inizio
Added compatibility with dictionary chr_pos -> s1,s2;A,C/sNew;A,T
'''

# argv1 è il file .bed con le annotazioni
# argv2 è il file .cluster.txt, che è ordinato per cromosoma. Era (03/03) il file top1 ordinato per chr
# argv3 è nome del file in output
# argv4 è directory dei dizionari
#argv5 is pamfile
# argv 6 is max allowed mms
# argv 7 is genome reference directory (Eg ../../Genomes/hg38_ref)
# argv8 is guide file
# argv9 is max allowed DNA bulges
# argv10 is max allowed RNA bulges
# argv11 is absolute path of sample file to load dictionary sample -> pop
# NOTE 06/03  -> removed PAM Disruption calculation
# NOTE 29/03 -> le colonne min max sono rimosse, dal file total.cluster sono già presenti colonne sample, annotation, real guide
# 29/03 colonne in input #Bulge_type     crRNA   DNA     Chromosome      Position        Cluster Position        Direction       Mismatches      Bulge_Size      Total   PAM_gen Var_uniq        Samples Annotation Type Real Guide
# NOTE1 can happend that a iupac falls under the N char of NGG, meaning that a target can have the same number of mms both in his REF and VAR part:
# CACTGCAACCTCTGTCTCCCKGG
# CACTGCAACCTCTGTCTCCCGGG        REF
# CACTGCAACCTCTGTCTCCCTGG        VAR
# So check if decrease_ref_count is not empty to avoid this (in this case +1 will be added to samples for VAR part of target and +1 for all the
# other samples for REF part)

import pickle
from typing import final  # to read CFD matrices
from supportFunctions.loadSample import associateSample
import re
import multiprocessing
import string
import numpy as np
import pandas as pd
import subprocess
import concurrent.futures
from intervaltree import Interval, IntervalTree
import os
import itertools
import time
import json
import sys
import warnings
from CRISTA_score import CRISTA_predict_list
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
# import azimuth.model_comparison
SIZE_DOENCH = 10000
N_THR = 3

# Return max doench value among list of extended targets
# def doenchParallel(targets, model, result):
#     doench_score =  azimuth.model_comparison.predict(targets,None, None, model= model, pam_audit=False)
#     doench_score = [np.around(i * 100) for i in doench_score]
#     max_doench = int(max(doench_score))
#     result.append(max_doench)


def doenchForIupac(sequence_doench, guide_seq, genome_type):
    pos_iupac = []
    var = []
    for pos, c in enumerate(sequence_doench):
        if c in iupac_code:
            pos_iupac.append(pos)
            var.append(iupac_code[c])

    if var:
        for i in itertools.product(*var):
            t = list(sequence_doench)
            for p, el in enumerate(pos_iupac):
                t[el] = i[p]
            targets_for_doench[guide_seq][genome_type].append(''.join(t))
    else:
        targets_for_doench[guide_seq][genome_type].append(sequence_doench)


def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open(os.path.dirname(
            os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
        pam_scores = pickle.load(open(os.path.dirname(
            os.path.realpath(__file__)) + '/PAM_scores.pkl', 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception(
            "Could not find file with mismatch scores or PAM scores")


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', '-': '-'}
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None  # If some IUPAC were not translated
    return ''.join(letters)

# Calculates CFD score


def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores, do_scores):
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
                # print(key)
            except Exception as e:
                score = 0
                break
            try:
                if 'N' == sl:
                    score *= 1
                else:
                    score *= mm_scores[key]
            except Exception as e:  # If '-' is in first position, i do not have the score for that position
                pass
    # print(pam)
    if 'N' in pam:
        score *= 1
    else:
        score *= pam_scores[pam]
    return score


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


def convert_fasta_to_dict(f):
    fasta = {}
    with open(f) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = ''
                continue
            sequence = line
            fasta[active_sequence_name] = sequence
    return fasta


def alignRefFromVar(line):  # chr_fake, start_pos, len_guide, bulge):
    t = line.copy()
    #chr_fake = t[10].split('_')
    len_guide = len(t[2])
    start_pos = int(t[4])
    true_chr = t[3]
    good_chr_fake = true_chr+':'+str(start_pos)+'-'+str(start_pos+len_guide)
    # dict_ref_seq[good_chr_fake].upper() #file_fasta.readline().strip().upper()
    sequence = genomeStr[start_pos:start_pos+len_guide].upper()
    # print('ref sequence', sequence, len(sequence), len_guide)

    if t[6] == '-':  # if negative strand reverse ref sequence
        # target = t[2][::-1]
        sequence = reverse_complement_table(sequence)
        target = t[2]
    else:
        target = t[2]

    if t[0] == "RNA":  # if RNA bulges in sequences find gaps and update ref sequence
        tmp_gap_position = [g.start() for g in re.finditer('-', target)]
        sequence = list(sequence)
        # cut sequence to extract only the part interested by counting the RNA bulges
        sequence = sequence[len(tmp_gap_position):]
        for tmp_g_p in tmp_gap_position:
            # sequence[tmp_g_p] = '-'
            sequence.insert(tmp_g_p, '-')
        # if t[6] == '-':
        #     sequence = reverse_complement_table(''.join(sequence))
        # else:

        # sequence = sequence[len(tmp_gap_position):]
        sequence = ''.join(sequence)
        # if t[6] == '+':
        #     sequence = sequence[0:len_guide]
        # else:
        #     sequence = reverse_complement_table(''.join(sequence))
        #     sequence = sequence[len(tmp_gap_position):]
    # elif t[6] == '-':
    #     sequence = reverse_complement_table(sequence)
    guide_no_pam = t[1][pos_beg:pos_end]
    list_t = list(sequence)

    # print('ref seq with bulges', sequence)
    # align ref sequence with var
    for position_t, char_t in enumerate(sequence[pos_beg:pos_end]):
        if char_t.upper() != guide_no_pam[position_t]:
            if guide_no_pam[position_t] != '-':
                list_t[sum_for_mms + position_t] = char_t.lower()
    sequence = ''.join(list_t)

    # print('ref seq realigned', sequence)
    # t[2] = sequence
    # if t[0] == 'DNA':
    #     if do_scores:
    #         cfd_score = calc_cfd(t[1][int(line[bulge_pos]):], t[2].upper()[
    #                              int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
    #     else:
    #         cfd_score = -1
    #     cfd = "{:.3f}".format(cfd_score)  # str(cfd_score)
    # else:
    #     if do_scores:
    #         cfd_score = calc_cfd(t[1], t[2].upper()[:-3],
    #                              t[2].upper()[-2:], mm_scores, pam_scores)
    #     else:
    #         cfd_score = -1
    #     cfd = "{:.3f}".format(cfd_score)  # str(cfd_score)
    # return [sequence, cfd]
    return sequence

# print('ESECUZIONE DI ANNOTATION E CALC SAMPLE INSIEME')
# print('TEST PER ANNOTAZIONE COMPLETA: I TARGET SENZA ANNOTAZIONE SONO SALVATI COME \"n\"')
# print('SE UN  TARGET HA 1+ ANNOTAZIONI, LE SALVA IN SINGOLA UNICA RIGA')
# print('RIMOZIONE DEI TARGET CHE NON HANNO SAMPLES')
# print('CALCOLO SCORES')
# print('SOSTITUZIONE IUPAC DI TUTTI I TARGET CON CHAR DEL TOP1SCOMPOSTO')
# print("READING INPUT FILES")
# Dictionaries for annotating samples

# Dict for populations
# pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
# all_samples = pop_file.Sample.to_list()
# all_pop = pop_file.Population.to_list()
# dict_sample_to_pop = dict()
# for  pos, i in enumerate(all_samples):
#     try:
#         dict_sample_to_pop[i] = all_pop[pos]        #{'S1':'POP1', 'S2':'POP1', ...}
#     except:
#         dict_sample_to_pop[i] = all_pop[pos]

# #Dict for superpopulation
# dict_pop_to_sup = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
#                     'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
#                     'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
#                     'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
#                     'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
# }
# superpopulation = ['EAS', 'EUR', 'AFR', 'AMR','SAS']

#dict_sample_to_pop, dict_pop_to_sup, dict_superpop_to_pop, dict_pop_to_sample, all_samples, all_pop, superpopulation, gender_sample = associateSample.loadSampleAssociation(sys.argv[11])


# READ INPUT FILES
annotationFile = sys.argv[1]  # file with annotation
resultsFile = sys.argv[2]  # file with results from search
outputFile = sys.argv[3]  # file with annotated results

# Get pam and guide length for new count mismatch samples
pam_at_beginning = False
with open(sys.argv[5]) as pam:
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

do_scores = True
if guide_len != 20 or len_pam != 3 or pam_at_beginning:
    # sys.stderr.write('CFD SCORE IS NOT CALCULATED WITH GUIDES LENGTH != 20 OR PAM LENGTH !=3 OR UPSTREAM PAM')
    do_scores = False
    with open(outputFile + '.acfd.txt', 'w+') as result:
        result.write('NO SCORES')

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


# OPEN INPUT FILES AND PREPARE OUTPUT FILE
inResult = open(resultsFile, "r")  # resultfile open
# inAnnotationFile = open(annotationFile, "r")  # file with annotations open
# outFileSampleAll = open(outputFile + '.samples.all.annotation.txt', 'w')  # outfile open (file with IUPAC targets and associated samples and annotation)

count_removed_target = 0
process = subprocess.Popen(['wc', '-l', resultsFile],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
total_line = int(out.decode('UTF-8').split(' ')[0])
if total_line < 1:
    print('WARNING! Input file has no targets')
    sys.exit()
if total_line < 10:
    mod_tot_line = 1
else:
    mod_tot_line = int(total_line/10)
# VARIABLE INIT
guideDict = {}
totalDict = {}

start_time = time.time()

#print("EXECUTING PRELIMINARY OPERATIONS")

# annotationsTree = IntervalTree()
# annotationsSet = set()
# #guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

# for line in inAnnotationFile:
#     x = line.split('\t')
#     x[3] = str(x[3]).rstrip("\n")
#     annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+str(x[3])
#     annotationsSet.add(str(x[3]))

# totalDict['targets'] = [0]*10
# for item in annotationsSet:
#     totalDict[item] = [0]*10

#print("PRELIMINARY OPERATIONS COMPLETED IN: %s seconds" % (time.time() - start_time))

start_time = time.time()

#print("EXECUTING ANNOTATION")

with open(resultsFile, 'r') as resFile:
    header_len = len(resFile.readline().strip().split('\t'))

# if header_len == 15:    #'Both' case : comparison variant/ref is active
# header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide'
header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference'


mm_pos = 7  # position of mismatch column
bulge_pos = 8
max_dna_bulges = int(sys.argv[9])
max_rna_bulges = int(sys.argv[10])
max_bulges = max_dna_bulges
if max_rna_bulges > max_bulges:
    max_bulges = max_rna_bulges

blank_add_begin = ' '  # Needed for replacing IUPAC in cluster targets
blank_add_end = ''
pam_multiplier = 1
pam_multiplier_negative = 0
start_sample_for_cluster = 0
cluster_step = 1  # If PAM end, go left to right
# when updatig lowercase for nem_mm, this value represents the offset for the pam position (mainly needed only if pam at beginning)
sum_for_mms = 0
# Values to check new iupac when working on cluster targets
end_sample_for_cluster = max_dna_bulges + max_rna_bulges
if pam_at_beginning:
    blank_add_begin = ''
    blank_add_end = ' '
    pam_multiplier = 0  # Since ' ' are at end, and '-' to reinsert are before the ' ', need to put max_dna_bulges and rna_bulges of target to 0
    pam_multiplier_negative = 1
    # For PAM at beginning, start from last nucleotide and go to left
    end_sample_for_cluster = len_pam + guide_len - max_rna_bulges
    start_sample_for_cluster = len_pam + guide_len + max_dna_bulges
    cluster_step = -1  # If PAM beginning, go right to left
    sum_for_mms = len_pam
# outFileSampleAll.write(header + '\n')
summary_samples = True

header_list = header.strip().split('\t')
# Variables for summary samples code
'''
{
    GUIDE1 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
    GUIDE2 -> {
        SAMPLE/POP/SUPERPOP1 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        },
        SAMPLE/POP/SUPERPOP2 ->{
            targets -> [0 0 0 0 0 0 0 0 0],
            ann1 -> [0 0 0 0 0 0 0 0 0],
            ann2 -> [0 0 0 0 0 0 0 0 0],
        }
    }
}

Per pop e superpop, se ho due sample stessa famiglia stesso target, conto solo una volta (visited_pop and visited_superpop array)
'''
count_sample = dict()  # NOTE cout_sample -> GUIDE -> SAMPLE -> has targets + ann1 + ann2 ... + refposition
# refposition is a key unrelated to the other keys (targets ann1 ann2 ...) and it's used to classify the sample (0 0+ 1 1+).
# it's put in here just to avoid to duplicate the entire guide -> sample ->     structure
# refposition -> [class , number of specific VAR on target to add/remove]  #Save class (0 at start) and number of ontarget var specific
# for that sample.
# Count number of REF target and REF part of semicommon
ontarget_reference_count = dict()
count_pop = dict()
# NOTE added key 'distributions' for population distribution images
count_superpop = dict()
# count_superpop-> GUIDE -> SUPERPOP -> targets ann1 ann2 ... distributions
# distributions is an array of len mms+bulge, each position contains an array [0,0,0] of len bulge+1 (indicating no bulge, 1 bulge, 2bulge ...)

# Create -Summary_total for a file ref.Annotation.summary.txt from the y and n values of Var_uniq column
summary_barplot_from_total = False
if 'Var_uniq' in header:
    summary_barplot_from_total = True
    vu_pos = header_list.index('Var_uniq')
# count_unique = dict()
# count_unique['targets'] = [0]*10
# count_unique_for_guide = dict()
# for item in annotationsSet:
#     count_unique[item] = [0]*10

# Variables for samples calculation
total_error = 0


current_chr = 'none'
chr_name = 'none'


def rev_comp(a):
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'


def preprocess_CFD_score(target):
    # preprocess target then calculate CFD score
    if do_scores:
        if target[0] == 'DNA':
            cfd_score = calc_cfd(target[1][int(target[bulge_pos]):], target[2].upper()[int(
                target[bulge_pos]):-3], target[2].upper()[-2:], mm_scores, pam_scores, do_scores)
            target.append("{:.3f}".format(cfd_score))
            if target[-3] == 55:
                target[-3] = "{:.3f}".format(cfd_score)
            if target[-3] == 33:
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

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq/2)-14:int(len_DNA_seq/2)]
        second_half = complete_DNA_seq[int(
            len_DNA_seq/2):int(len_DNA_seq/2)+15]
        complete_DNA_seq = first_half+second_half
        if target[6] == '-':
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

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

        # trim the 3' and 5' end to avoid sequences longer than 29
        len_DNA_seq = len(complete_DNA_seq)
        first_half = complete_DNA_seq[int(len_DNA_seq/2)-14:int(len_DNA_seq/2)]
        second_half = complete_DNA_seq[int(
            len_DNA_seq/2):int(len_DNA_seq/2)+15]
        complete_DNA_seq = first_half+second_half
        if target[6] == '-':
            complete_DNA_seq = reverse_complement_table(complete_DNA_seq)

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

    # time_start = time.time()
    # print('CALCULATING CFD SCORE')

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


mm_scores, pam_scores = get_mm_pam_scores()
guides_dict = dict()  # For CFD score
guides_dict_doench = dict()
targets_for_doench = dict()

N_THR = multiprocessing.cpu_count() // 2
refgenomedir = sys.argv[7]

# with open( os.path.dirname(os.path.realpath(__file__)) + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
#     model = pickle.load(f)
max_doench = 0
sum_cfd = 0
cfd_scores = []

# open fasta to read ref chromosome
inFasta = open(sys.argv[7], 'r')
current_chr = inFasta.readline().strip().replace('>', '')  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = ''.join(genomeStr).upper()
# string of the whole chromosome on single line
genomeStr = genomeStr.replace('\n', '')


start_time_total = time.time()
lines_processed = 0
allowed_mms = int(sys.argv[6])
current_guide_chr_pos_direction = 'no'

# save best files
cfd_best = open(outputFile + '.bestCFD_INDEL.txt', 'w')
cfd_best.write(header + '\tCFD\n')  # Write header

crista_best = open(outputFile+'.bestCRISTA_INDEL.txt', 'w')
crista_best.write(header+'\tCFD\n')

mmblg_best = open(outputFile + '.bestmmblg_INDEL.txt', 'w')
mmblg_best.write(header + '\tCFD\n')  # Write header


INDELS_tree = IntervalTree()
with open(os.path.realpath(sys.argv[4]) + '/log' + current_chr + '.txt', 'r') as log:
    log.readline()
    for entry in log:
        splitted = entry.strip().split('\t')
        fake_pos = splitted[5].strip().split(
            ',')  # inizio,fine del fakechr
        # start       end fakechr          #pos_vera,samples,rsID,MAF,INDELinfo,start_fake,ref_seq
        INDELS_tree[int(fake_pos[0]):int(fake_pos[1])] = [
            splitted[0], splitted[1], splitted[2], splitted[3], splitted[4], int(fake_pos[0]), splitted[6]]
#datastore = datastore.to_dict(orient='index')
print('Analysis of ' + current_chr)

save_cluster_targets = True
remove_iupac = False
save_total_general_table = False
# chiave ['add'] = target semicommon da aggiungere alla tab generale delle guide, metto i valori di total presi dal primo target REF nel cluster di un semicommon che esiste
add_to_general_table = dict()
# chiave ['distribution'] = array len total, ogni cella divisa per mm,1B,2B..., per populationDistribution
last_annotation = ''  # in un cluster, l'annotazione è la stessa, quindi la calcolo solo una volta e poi la riscrivo per gli altri target del cluster
last_samples = []  # se due target hanno la stessa scomposizione, hanno gli stessi samples, quindi non li ricalcolo
next(inResult)  # Skip header
# contains all targets for each cluster, it's then sorted by CFD. Only the target with highest CFD is saved
# cluster_to_save = []
# sum number times a cfd value appears, TODO aggiungere anche il conteggio per ogni sample
cfd_for_graph = {'ref': [0]*101, 'var': [0]*101}

# Categorie per la scomposizione, chiavi: 'BulgeType' + 'BulgeSize'
dictionary_entries = ['X0']
for i in range(int(max_dna_bulges)):
    dictionary_entries.append('DNA' + str(i+1))
for i in range(int(max_rna_bulges)):
    dictionary_entries.append('RNA' + str(i+1))


#dict_ref_seq = convert_fasta_to_dict(sys.argv[12])

cluster_class = None
datastore = None


global_start = time.time()

#list saving all the targets reported in chr##
cluster_to_save = list()

for line in inResult:
    # print('line in INDEL analysis', line)

    line = line.strip().split('\t')
    # print(line)
    guide_no_bulge = line[1].replace('-', '')
    # copy line to avoid problem with memory managament and intrisic pointers
    final_result = line.copy()
    # extract indel data from the INDEL tree
    try:
        indel_data = sorted(INDELS_tree[int(line[4])])[0].data
    except:
        # target found has no INDEL, so it's REF and must be skipped
        continue
    # assign extracted data to final target
    final_result[3] = current_chr  # current analyzed chr
    final_result[12] = indel_data[1]  # samples
    final_result[15] = indel_data[2]  # rsID
    final_result[16] = indel_data[3]  # MAF
    final_result[17] = indel_data[4]  # INDELinfo

    # correct position if PAM at beginning and strand
    if pam_at_beginning:
        if line[0] == 'RNA' and line[6] == '-':
            # cluster_position
            fake_start_target = int(line[5]) - int(indel_data[5])
        else:
            # real_position
            fake_start_target = int(line[4]) - int(indel_data[5])
    else:
        if line[0] == 'RNA' and line[6] == '+':
            # cluster_position
            fake_start_target = int(line[5]) - int(indel_data[5])
        else:
            # real_position
            fake_start_target = int(line[4]) - int(indel_data[5])

    # correct start with real chr spatial info and length of the target
    true_start_target = int(indel_data[0].split(
        '_')[1].split('-')[0]) + fake_start_target
    # correct cluster position
    diff_pos_clus = int(line[4]) - int(line[5])
    # real_start
    final_result[4] = str(true_start_target)
    # real_start_cluster
    final_result[5] = str(true_start_target - diff_pos_clus)

    # real_target
    # t = final_result[2]
    # mm_new_t = 0
    # tmp_pos_mms = None  # lista posizione dei mms
    # guide_no_pam = line[1][pos_beg:pos_end]
    # list_t = list(t)
    # for position_t, char_t in enumerate(t[pos_beg:pos_end]):
    #     if char_t.upper() != guide_no_pam[position_t]:
    #         mm_new_t += 1
    #         tmp_pos_mms = position_t
    #         if guide_no_pam[position_t] != '-':
    #             list_t[sum_for_mms + position_t] = char_t.lower()
    # final_result[2] = ''.join(list_t)  # t

    # print('target with no process', final_result[2])
    # if tmp_pos_mms:
    #     final_result.append(int(tmp_pos_mms))
    # else:  # NO mms found
    #     final_result.append(-1)

    # print('results before calc and ref alignment', final_result)

    final_result.append(alignRefFromVar(final_result))
    # number to activate ref score calculation (active if target is alternative)
    final_result.append(33)
    # position of tmp_mms (removed later after processing)
    final_result.append(0)
    cluster_to_save.append(final_result)

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
    # skip processing and return empty file
    # close open files
    cfd_best.close()
    mmblg_best.close()
    crista_best.close()
    # update header
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCFD_INDEL.txt')
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCRISTA_INDEL.txt')
    os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestmmblg_INDEL.txt')
    # print exit messagge and exit 0
    print('Done', current_chr, time.time() - start_time)
    print('ANALYSIS COMPLETE IN', time.time() - global_start)
    exit(0)

# start processing scores after the whole target file is read
clusters_with_scores = calculate_scores(cluster_to_save)

for count, cluster in enumerate(clusters_with_scores):
    for target in cluster:
        if count == 0:  # CFD target
            # remove count of tmp_mms
            target.pop(-2)
            # save CFD targets
            # print('CFD', target)
            cfd_best.write(
                '\t'.join(target)+'\t'+str(0)+'\n')
            # save mm-bul targets
            mmblg_best.write(
                '\t'.join(target)+'\t'+str(0)+'\n')
        if count == 1:  # CRISTA target
            # remove count of tmp_mms
            target.pop(-2)
            # print('CRISTA', target)
            # save CRISTA targets
            crista_best.write('\t'.join(target)+'\t' +
                              str(0)+'\n')

cfd_best.close()
mmblg_best.close()
crista_best.close()

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCFD_INDEL.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestCRISTA_INDEL.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster/' "+outputFile + '.bestmmblg_INDEL.txt')


print('Done', current_chr, time.time() - start_time)

print('ANALYSIS COMPLETE IN', time.time() - global_start)
