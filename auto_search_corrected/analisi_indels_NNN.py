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

#argv1 è il file .bed con le annotazioni
#argv2 è il file .cluster.txt, che è ordinato per cromosoma. Era (03/03) il file top1 ordinato per chr
#argv3 è nome del file in output
#argv4 è directory dei dizionari
#argv5 is pamfile
#argv 6 is max allowed mms
#argv 7 is genome reference directory (Eg ../../Genomes/hg38_ref)
#argv8 is guide file
#argv9 is max allowed DNA bulges
#argv10 is max allowed RNA bulges
#argv11 is absolute path of sample file to load dictionary sample -> pop
# NOTE 06/03  -> removed PAM Disruption calculation
#NOTE 29/03 -> le colonne min max sono rimosse, dal file total.cluster sono già presenti colonne sample, annotation, real guide
# 29/03 colonne in input #Bulge_type     crRNA   DNA     Chromosome      Position        Cluster Position        Direction       Mismatches      Bulge_Size      Total   PAM_gen Var_uniq        Samples Annotation Type Real Guide
#NOTE1 can happend that a iupac falls under the N char of NGG, meaning that a target can have the same number of mms both in his REF and VAR part:
#CACTGCAACCTCTGTCTCCCKGG
#CACTGCAACCTCTGTCTCCCGGG        REF
#CACTGCAACCTCTGTCTCCCTGG        VAR
#So check if decrease_ref_count is not empty to avoid this (in this case +1 will be added to samples for VAR part of target and +1 for all the
# other samples for REF part)

import warnings 
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore',category=UserWarning)
import sys
import json
import time
import itertools
import os
from intervaltree import Interval, IntervalTree
import concurrent.futures
import subprocess
import pandas as pd
import pickle       #to read CFD matrices
import numpy as np
import azimuth.model_comparison
import string
import multiprocessing
import re
from supportFunctions.loadSample import associateSample
SIZE_DOENCH = 10000
N_THR = 3

#Return max doench value among list of extended targets
def doenchParallel(targets, model, result):
    doench_score =  azimuth.model_comparison.predict(targets,None, None, model= model, pam_audit=False)
    doench_score = [np.around(i * 100) for i in doench_score]
    max_doench = int(max(doench_score))
    result.append(max_doench)

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
        mm_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
        pam_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) +'/PAM_scores.pkl', 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', '-':'-'} 
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None     #If some IUPAC were not translated
    return ''.join(letters)

# Calculates CFD score
def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores):
    score = 1
    sg = sg.replace('T', 'U')
    guide_seq = guide_seq.replace('T', 'U')
    s_list = list(sg)
    guide_seq_list = list(guide_seq)
    
    for i, sl in enumerate(s_list):  
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            try:    #Catch exception if IUPAC character
                key = 'r' + guide_seq_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
                #print(key)
            except Exception as e:
                score = 0
                break
            try:
                if 'N' == sl:
                    score *= 1
                else:
                    score *= mm_scores[key]
            except Exception as e : #If '-' is in first position, i do not have the score for that position
                pass
    #print(pam)
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

def alignRefFromVar(line, dict_ref_seq):#chr_fake, start_pos, len_guide, bulge):
    t = line.copy()
    #chr_fake = t[10].split('_')
    len_guide = len(t[2])
    start_pos = int(t[4])
    true_chr = t[3]
    good_chr_fake = true_chr+':'+str(start_pos)+'-'+str(start_pos+len_guide)
    sequence = dict_ref_seq[good_chr_fake].upper() #file_fasta.readline().strip().upper()
    if t[6] == '-':
        target = t[2][::-1]
    else:
        target = t[2]
    if t[0] == "RNA":
        tmp_gap_position = [g.start() for g in re.finditer('-', target)]
        sequence = list(sequence)
        for tmp_g_p in tmp_gap_position:
            sequence.insert(tmp_g_p, '-')
        if t[6] == '+':
            sequence = sequence[0:len_guide]
        else:
            sequence = reverse_complement_table(''.join(sequence))
            sequence = sequence[len(tmp_gap_position):]
    elif t[6] == '-':
        sequence = reverse_complement_table(sequence)
    guide_no_pam = t[1][pos_beg:pos_end]  
    list_t = list(sequence)  
    for position_t, char_t in enumerate(sequence[pos_beg:pos_end]): 
        if char_t.upper() != guide_no_pam[position_t]:
            if guide_no_pam[position_t] != '-':
                list_t[sum_for_mms + position_t] = char_t.lower()
    sequence = ''.join(list_t)
    t[2] = sequence
    if t[0] == 'DNA':
        cfd_score = calc_cfd(t[1][int(line[bulge_pos]):], t[2].upper()[int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
        cfd = str(cfd_score)
    else:
        cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
        cfd = str(cfd_score)
    return [sequence, cfd]


'''
def alignRefFromVar(line):#chr_fake, start_pos, len_guide, bulge):
    t = line.copy()
    chr_fake = t[10]
    len_guide = len(t[1])
    start_pos = int(t[5])
    true_chr = t[3]
    file_bed = open(true_chr+"_tmp_fake.bed", 'w')
    #file_bed.write("CHROM\tSTART\tEND\n")
    file_bed.write(true_chr+"\t"+str(start_pos)+"\t"+str(start_pos+len_guide)+"\n")
    file_bed.close()
    os.system("bedtools getfasta -fi "+sys.argv[7]+"/"+true_chr+".fa -bed "+true_chr+"_tmp_fake.bed -fo "+true_chr+".true.fa")
    file_fasta = open(true_chr+".true.fa", 'r')
    header = file_fasta.readline()
    sequence = file_fasta.readline().strip().upper()
    if t[6] == '-':
        sequence = reverse_complement_table(sequence)
    file_fasta.close()
    if t[0] == "RNA":
        tmp_gap_position = [g.start() for g in re.finditer('-', t[2])]
        sequence = list(sequence)
        for tmp_g_p in tmp_gap_position:
            sequence.insert(tmp_g_p, '-')
        sequence = ''.join(sequence[0:len_guide])
    guide_no_pam = t[1][pos_beg:pos_end]  
    list_t = list(sequence)  
    for position_t, char_t in enumerate(sequence[pos_beg:pos_end]): 
        if char_t.upper() != guide_no_pam[position_t]:
            if guide_no_pam[position_t] != '-':
                list_t[sum_for_mms + position_t] = char_t.lower()
    sequence = ''.join(list_t)
    t[2] = sequence
    if t[0] == 'DNA':
        cfd_score = calc_cfd(t[1][int(line[bulge_pos]):], t[2].upper()[int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
        cfd = str(cfd_score)
    else:
        cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
        cfd = str(cfd_score)
    os.system("rm "+true_chr+"_tmp_fake.bed")
    os.system("rm "+true_chr+".true.fa")
    return [sequence, cfd]
'''

print('ESECUZIONE DI ANNOTATION E CALC SAMPLE INSIEME')
print('TEST PER ANNOTAZIONE COMPLETA: I TARGET SENZA ANNOTAZIONE SONO SALVATI COME \"n\"')
print('SE UN  TARGET HA 1+ ANNOTAZIONI, LE SALVA IN SINGOLA UNICA RIGA')
print('RIMOZIONE DEI TARGET CHE NON HANNO SAMPLES')
print('CALCOLO SCORES')
print('SOSTITUZIONE IUPAC DI TUTTI I TARGET CON CHAR DEL TOP1SCOMPOSTO')
print("READING INPUT FILES")
#Dictionaries for annotating samples

#Dict for populations
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

#READ INPUT FILES
annotationFile = sys.argv[1] #file with annotation
resultsFile = sys.argv[2] #file with results from search
outputFile = sys.argv[3] #file with annotated results

#Get pam and guide length for new count mismatch samples
pam_at_beginning = False
with open (sys.argv[5]) as pam:
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
if guide_len != 20:
    with open(outputFile + '.acfd.txt', 'w+') as result:
        result.write('NO SCORES')
        do_scores = False

iupac_code_set = {
          "R":{"A", "G"},
          "Y":{"C", "T"},
          "S":{"G", "C"},
          "W":{"A", "T"},
          "K":{"G", "T"},
          "M":{"A", "C"},
          "B":{"C", "G", "T"},
          "D":{"A", "G", "T"},
          "H":{"A", "C", "T"},
          "V":{"A", "C", "G"},
          "r":{"A", "G"},
          "y":{"C", "T"},
          "s":{"G", "C"},
          "w":{"A", "T"},
          "k":{"G", "T"},
          "m":{"A", "C"},
          "b":{"C", "G", "T"},
          "d":{"A", "G", "T"},
          "h":{"A", "C", "T"},
          "v":{"A", "C", "G"},
          "A":{"A"},
          "T":{"T"},
          "C":{"C"},
          "G":{"G"},
          "a":{"a"},
          "t":{"t"},
          "c":{"c"},
          "g":{"g"},
          'N':{'A','T','G','C'}
        }


#OPEN INPUT FILES AND PREPARE OUTPUT FILE
inResult = open(resultsFile, "r")  # resultfile open
inAnnotationFile = open(annotationFile, "r")  # file with annotations open
# outFileSampleAll = open(outputFile + '.samples.all.annotation.txt', 'w')  # outfile open (file with IUPAC targets and associated samples and annotation)

count_removed_target = 0
process = subprocess.Popen(['wc', '-l', resultsFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
total_line = int(out.decode('UTF-8').split(' ')[0])
if total_line < 1:
    print('WARNING! Input file has no targets')
    sys.exit()
if total_line < 10:
    mod_tot_line = 1
else:
    mod_tot_line = int(total_line/10)
#VARIABLE INIT
guideDict = {}
totalDict = {}

start_time = time.time()

print("EXECUTING PRELIMINARY OPERATIONS")

annotationsTree = IntervalTree()
annotationsSet = set()
#guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

for line in inAnnotationFile:
    x = line.split('\t')
    x[3] = str(x[3]).rstrip("\n")
    annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+str(x[3])
    annotationsSet.add(str(x[3]))

totalDict['targets'] = [0]*10
for item in annotationsSet:
    totalDict[item] = [0]*10

print("PRELIMINARY OPERATIONS COMPLETED IN: %s seconds" % (time.time() - start_time))

start_time = time.time()

print("EXECUTING ANNOTATION")

with open(resultsFile, 'r') as resFile:
    header_len = len(resFile.readline().strip().split('\t'))

# if header_len == 15:    #'Both' case : comparison variant/ref is active
#header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide'
header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tReference'


mm_pos = 7      #position of mismatch column
bulge_pos = 8
max_dna_bulges = int(sys.argv[9])
max_rna_bulges = int(sys.argv[10])
max_bulges = max_dna_bulges
if max_rna_bulges > max_bulges:
    max_bulges = max_rna_bulges

blank_add_begin = ' '               #Needed for replacing IUPAC in cluster targets
blank_add_end = ''
pam_multiplier = 1
pam_multiplier_negative = 0
start_sample_for_cluster = 0
cluster_step = 1    #If PAM end, go left to right
sum_for_mms = 0 #when updatig lowercase for nem_mm, this value represents the offset for the pam position (mainly needed only if pam at beginning)
end_sample_for_cluster = max_dna_bulges + max_rna_bulges  #Values to check new iupac when working on cluster targets
if pam_at_beginning:
    blank_add_begin = ''
    blank_add_end = ' '
    pam_multiplier = 0              #Since ' ' are at end, and '-' to reinsert are before the ' ', need to put max_dna_bulges and rna_bulges of target to 0
    pam_multiplier_negative = 1
    end_sample_for_cluster = len_pam + guide_len - max_rna_bulges       #For PAM at beginning, start from last nucleotide and go to left
    start_sample_for_cluster = len_pam + guide_len + max_dna_bulges
    cluster_step = -1   #If PAM beginning, go right to left
    sum_for_mms = len_pam
# outFileSampleAll.write(header + '\n')
summary_samples = True

header_list = header.strip().split('\t')
#Variables for summary samples code
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
count_sample = dict()       #NOTE cout_sample -> GUIDE -> SAMPLE -> has targets + ann1 + ann2 ... + refposition
    #refposition is a key unrelated to the other keys (targets ann1 ann2 ...) and it's used to classify the sample (0 0+ 1 1+).
    #it's put in here just to avoid to duplicate the entire guide -> sample ->     structure
    # refposition -> [class , number of specific VAR on target to add/remove]  #Save class (0 at start) and number of ontarget var specific
    #for that sample.
ontarget_reference_count = dict() #Count number of REF target and REF part of semicommon
count_pop = dict()
count_superpop = dict()     #NOTE added key 'distributions' for population distribution images
    #count_superpop-> GUIDE -> SUPERPOP -> targets ann1 ann2 ... distributions
    # distributions is an array of len mms+bulge, each position contains an array [0,0,0] of len bulge+1 (indicating no bulge, 1 bulge, 2bulge ...)

#Create -Summary_total for a file ref.Annotation.summary.txt from the y and n values of Var_uniq column
summary_barplot_from_total = False
if 'Var_uniq' in header:
    summary_barplot_from_total = True
    vu_pos = header_list.index('Var_uniq')
count_unique = dict()
count_unique['targets'] = [0]*10
count_unique_for_guide = dict()
for item in annotationsSet:
    count_unique[item] = [0]*10

#Variables for samples calculation
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

iupac_code = {
            "R":("A", "G"),
            "Y":("C", "T"),
            "S":("G", "C"),
            "W":("A", "T"),
            "K":("G", "T"),
            "M":("A", "C"),
            "B":("C", "G", "T"),
            "D":("A", "G", "T"),
            "H":("A", "C", "T"),
            "V":("A", "C", "G"),
            "r":("A", "G"),
            "y":("C", "T"),
            "s":("G", "C"),
            "w":("A", "T"),
            "k":("G", "T"),
            "m":("A", "C"),
            "b":("C", "G", "T"),
            "d":("A", "G", "T"),
            "h":("A", "C", "T"),
            "v":("A", "C", "G"),
            'N':('A', 'T', 'C', 'G')
            }

#For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv", "TGACYRSWKMDHVBtgacyrswkmdhvb") 

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

mm_scores, pam_scores = get_mm_pam_scores()
guides_dict = dict()    #For CFD score
guides_dict_doench = dict()
targets_for_doench = dict()

N_THR = multiprocessing.cpu_count() // 2
refgenomedir = sys.argv[7]

with open( os.path.dirname(os.path.realpath(__file__)) + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
    model = pickle.load(f)
max_doench = 0
sum_cfd = 0
cfd_scores = []


start_time_total = time.time()
lines_processed = 0
allowed_mms = int(sys.argv[6])
current_guide_chr_pos_direction = 'no'
cfd_best = open(outputFile + '.bestCFD.txt', 'w+')
cfd_best.write(header + '\tCFD\n')        #Write header

cfd_alt = open(outputFile + '.altCFD.txt', 'w+')
cfd_alt.write(header + '\tCFD\n')        #Write header

mmblg_best = open(outputFile + '.bestmmblg.txt', 'w+')
mmblg_best.write(header + '\tCFD\n')        #Write header

mmblg_alt = open(outputFile + '.altmmblg.txt', 'w+')
mmblg_alt.write(header + '\tCFD\n')        #Write header

save_cluster_targets = True
remove_iupac = False
save_total_general_table = False
add_to_general_table = dict()   #chiave ['add'] = target semicommon da aggiungere alla tab generale delle guide, metto i valori di total presi dal primo target REF nel cluster di un semicommon che esiste
                                #chiave ['distribution'] = array len total, ogni cella divisa per mm,1B,2B..., per populationDistribution
last_annotation = ''    #in un cluster, l'annotazione è la stessa, quindi la calcolo solo una volta e poi la riscrivo per gli altri target del cluster
last_samples = []       #se due target hanno la stessa scomposizione, hanno gli stessi samples, quindi non li ricalcolo
next(inResult)      #Skip header
cluster_to_save = []        #contains all targets for each cluster, it's then sorted by CFD. Only the target with highest CFD is saved
cfd_for_graph = {'ref':[0]*101, 'var':[0]*101}  #sum number times a cfd value appears, TODO aggiungere anche il conteggio per ogni sample

dictionary_entries = ['X0']                     #Categorie per la scomposizione, chiavi: 'BulgeType' + 'BulgeSize'
for i in range(int(max_dna_bulges)):
    dictionary_entries.append('DNA' + str(i+1))
for i in range(int(max_rna_bulges)):
    dictionary_entries.append('RNA' + str(i+1))


dict_ref_seq = convert_fasta_to_dict(sys.argv[12])  

cluster_class = None
datastore = None


global_start = time.time()

for line in inResult:
    line = line.strip().split('\t')
    # print(line)
    guide_no_bulge = line[1].replace('-','')
    activate_cluster_scomposition = False   #se ho iupac in posizioni che non c'erano nel top1, fai la scomposizione completa
    if (guide_no_bulge + line[3] + line[5] + line[6]) == current_guide_chr_pos_direction:      #In current cluster
        if cluster_class[line[0] + line[8]] is None:
            activate_cluster_scomposition = True
        if activate_cluster_scomposition:
            final_result = line.copy()
            var = False
            if line[10] != "n":
                a = datastore[line[10]]
                var = True
                final_result[13] = a['SAMPLES']
                final_result[16] = a['rsID']
                final_result[17] = str(a['AF'])
                final_result[18] = a['indel']
                
            t = final_result[2]
            mm_new_t = 0
            tmp_pos_mms =  None       #lista posizione dei mms
            guide_no_pam = line[1][pos_beg:pos_end]  
            list_t = list(t)  
            for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != '-':
                        list_t[sum_for_mms + position_t] = char_t.lower()
            final_result[2] = ''.join(list_t)#t
            if tmp_pos_mms:
                final_result.append(int(tmp_pos_mms))
            else: #NO mms found
                final_result.append(-1)
            
            final_result[14] = last_annotation
            
            if not var:
                final_result[13] = "n" #last_samples
                final_result[16] = "n"#last_IDs
                final_result[17] = "n"#last_AF
                final_result[18] = "n"#last_INDpos
            
            last_samples = final_result[13]
            last_IDs = final_result[16]
            last_AF = final_result[17]
            last_INDpos = final_result[18]
            
            cluster_to_save.append(final_result)
            
            clusterkey = line[0] + line[8]
            cluster_class[clusterkey] = dict()
            #cluster_class[clusterkey]['poscharforcluster'] = pos_char_for_cluster
            cluster_class[clusterkey]['lastsamples'] = last_samples
            cluster_class[clusterkey]['lastID'] = last_IDs #ID
            cluster_class[clusterkey]['lastAF'] = last_AF #AF (come last_samples)
            cluster_class[clusterkey]['lastINDpos'] = last_INDpos
            #cluster_class[clusterkey]['ref_pos_alignement'] = ref_pos_alignement
            #cluster_class[clusterkey]['lastHaplotype'] = last_haplotype
          
        else:
            #Use previous scompositions to create possible scomposition in cluster
            clusterkey = line[0] + line[8]
            #pos_char_for_cluster = cluster_class[clusterkey]['poscharforcluster']
            last_samples = cluster_class[clusterkey]['lastsamples']
            last_IDs = cluster_class[clusterkey]['lastID'] 
            last_AF = cluster_class[clusterkey]['lastAF']
            last_INDpos = cluster_class[clusterkey]['lastINDpos']
            #last_haplotype = cluster_class[clusterkey]['lastHaplotype']
            t = line[2]
            #if line[6] == '-':
            #    t = t[::-1]
            mm_new_t = 0
            tmp_pos_mms =  None       #lista posizione dei mms
            guide_no_pam = line[1][pos_beg:pos_end]  
            list_t = list(t)  
            for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != '-':
                        list_t[sum_for_mms + position_t] = char_t.lower()
            line[2] = ''.join(list_t)#t
            if tmp_pos_mms:
                line.append(int(tmp_pos_mms))
            else: #NO mms found
                line.append(-1)
            
            line[13] = last_samples
            line[14] = last_annotation
            line[16] = last_IDs
            line[17] = last_AF
            line[18] = last_INDpos
            
            cluster_to_save.append(line)
    else:       #New cluster
        #New cluster, order previous by CFD, and save the highest CFD
        
        # check_position = 0

        #Calcolo lo score per ogni target nel cluster
        #print(len(cluster_to_save))
        for t in cluster_to_save:
           
            if t[0] == 'DNA':
                cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
                t.append(str(cfd_score))
            else:
                cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
                t.append(str(cfd_score))

        cluster_to_save.sort(key = lambda x : (float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse = True)
        cluster_to_save_mmbl = cluster_to_save.copy()
        cluster_to_save_mmbl.sort(key = lambda x : (int(x[8]), int(x[7])))
        
        keys_seen = []
        saved = False
        if cluster_to_save:
            c = cluster_to_save[0]

            c.pop(-2)
            #best_seq = c[12]
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            #print("KEY ORIGINALE", cfd_clus_key)
            if c[13] == 'n':
                cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                cfd_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
            else:
                cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                ref_generated = alignRefFromVar(c, dict_ref_seq)#, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                ref_this_clus = ref_generated[0]
                cfd_ref = ref_generated[1]
                reference_pam = ref_this_clus[pam_begin:pam_end]
                found_creation = False
                for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
                    if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:     #ref char not in set of general pam char
                        found_creation = True
                if found_creation:
                    c[11] = c[2][pam_begin:pam_end]
                #print(c, ref_this_clus, cfd_ref)
                cfd_best.write('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
            cluster_to_save.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True
            
        list_for_alt = []
        for c in cluster_to_save:
            c.pop(-2)
            #print(c[3], c[5], c[6], c[-2])
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            #print("POSSIBLE NEW KEY", new_cfd_clus_key)
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                if c[13] == 'n':
                    cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                    #cfd_alt.write('\t'.join(c)+'\tn\n')
                    list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
                else:
                    cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                    #cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
                    ref_generated = alignRefFromVar(c, dict_ref_seq)#, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                    ref_this_clus = ref_generated[0]
                    reference_pam = ref_this_clus[pam_begin:pam_end]
                    cfd_ref = ref_generated[1]
                    found_creation = False
                    for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
                        if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:     #ref char not in set of general pam char
                            found_creation = True
                    if found_creation:
                        c[11] = c[2][pam_begin:pam_end]
                    list_for_alt.append('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
        if saved:
            cfd_best.write('\t' + str(len(list_for_alt))+ '\n')
        for ele in list_for_alt:
            cfd_alt.write(ele+'\t'+str(len(list_for_alt))+ '\n')

        keys_seen = []
        saved = False
        if cluster_to_save_mmbl:
            c = cluster_to_save_mmbl[0]

            #c.pop(-2)
            best_seq = c[12]
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            
            if c[12] == 'n':
                #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                mmblg_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
            else:
                #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                ref_this_clus = ref_generated[2]
                cfd_ref = ref_generated[-1]
                mmblg_best.write('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
            cluster_to_save_mmbl.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True
            
        list_for_alt = []
        for c in cluster_to_save_mmbl:
            #c.pop(-2)
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                if c[12] == 'n':
                    #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                    #cfd_alt.write('\t'.join(c)+'\tn\n')
                    list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
                else:
                    #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                    #cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
                    ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                    ref_this_clus = ref_generated[2]
                    cfd_ref = ref_generated[-1]
                    list_for_alt.append('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
        if saved:
            mmblg_best.write('\t' + str(len(list_for_alt))+ '\n')
        for ele in list_for_alt:
            mmblg_alt.write(ele+'\t'+str(len(list_for_alt))+ '\n')

        cluster_to_save = []
        cluster_class = dict.fromkeys(dictionary_entries)       #Dizionario dove salvare le variabile per ogni tipo di scomposizione (X,DNA1,DNA2,RNA1 etc...)
                                                                #Viene re-inizializzato per ogni cluster
        current_guide_chr_pos_direction = guide_no_bulge + line[3] + line[5] + line[6]
        if line[3].split('_')[0] != current_chr:      #New chr section, load a new dictionary
            if not os.path.exists(os.path.realpath(sys.argv[4]) + '/log' + line[3].split('_')[0] + '.txt'):
                raise Exception("ERROR, NO LOG PASSED")
            else:
                #if current_chr != 'none':
                #    print('Done', current_chr, time.time() - start_time)
                current_chr = line[3].split('_')[0]
                #chr_name = line[3]
                datastore = pd.read_csv(os.path.realpath(sys.argv[4]) + '/log' + current_chr + '.txt', sep='\t', index_col=0).fillna('n')
                datastore = datastore.loc[~datastore.index.duplicated(keep='first')]
                datastore = datastore.to_dict(orient='index')
                print ('Analysis of ' + current_chr)
                start_time = time.time()

        #Calculate samples and all possible true scompositions
        pos_snp = []        # lista delle posizioni degli iupac nel target
        tuple_var_ref = []  #array che contiene le tuple composte da (var1,var2..., ref) (al momento abbiamo solo var1 nei dizionari) per ogni iupac nel target
        # top1_ref = False
        #target_combination = []     #lista di tutte le possibili combinazioni di scomposizione di IUPAC
        pos_snp_chr = []        # lista delle posizioni degli iupac nel target riferite alla posizione nel cromosoma
        #set_list_top1 = []      # lista di set di samples associati ad ogni posizione IUPAC
        last_samples = []
        last_IDs = []
        last_AF = []
        last_SNPpos = []
        #substitute_ref = []     #List of ref_char used to create the aligned ref target from a var
        #last_haplotype = []
        #pos_char_for_cluster = [] # lista di liste. Ogni lista corrisponde ad una scomposizione. Per ogni scomposizione mi salvo la posizione nel target e il carattere che ho trovato nella scomposizione
                                # andrò poi a sostituire questi caratteri negli altri target del cluster che sono nella stessa categoria
        final_result = line.copy()
        last_samples = []
        last_IDs = []
        last_AF = []
        last_INDpos = []
        if line[10] != "n":
            a = datastore[line[10]]
            final_result[13] = a['SAMPLES']
            final_result[16] = a['rsID']
            final_result[17] = str(a['AF'])
            final_result[18] = a['indel']
            
        mm_new_t = 0
        t = line[2]
        #if line[6] == '-':
        #    t = t[::-1]
        guide_no_pam = final_result[1][pos_beg:pos_end]  
        list_t = list(t)  
        for position_t, char_t in enumerate(final_result[2][pos_beg:pos_end]):
            if char_t.upper() != guide_no_pam[position_t]:
                mm_new_t += 1
                tmp_pos_mms = position_t
                if guide_no_pam[position_t] != '-':
                    list_t[sum_for_mms + position_t] = char_t.lower()
        final_result[2] = ''.join(list_t)#t
        if tmp_pos_mms:
            final_result.append(int(tmp_pos_mms))
        else: #NO mms found
            final_result.append(-1)
            
        #Calcolo annotazioni
        foundAnnotations = sorted(annotationsTree[int(line[4]):(int(line[4])+int(len(guide_no_bulge))+1)])
        string_annotation = []
        found_bool = False
        for found in range(0, len(foundAnnotations)):
            guide = foundAnnotations[found].data
            guideSplit = guide.split('\t')
            if(str(guideSplit[0]) == str(line[3])):
                found_bool = True
                string_annotation.append(str(guideSplit[1]))              
        if not found_bool:
            last_annotation = 'n'
        else:
            last_annotation = ','.join(string_annotation)

        final_result[14] = last_annotation
        last_samples= final_result[13]
        last_IDs = final_result[16]
        last_AF = final_result[17]
        last_INDpos = final_result[18]
        
        cluster_to_save.append(final_result)
        
        clusterkey = line[0] + line[8]
        cluster_class[clusterkey] = dict()
        #cluster_class[clusterkey]['poscharforcluster'] = pos_char_for_cluster       #aggiorno la categoria del target per i prossimi target
        cluster_class[clusterkey]['lastsamples'] = last_samples
        cluster_class[clusterkey]['lastID'] = last_IDs #ID
        cluster_class[clusterkey]['lastAF'] = last_AF #AF (come last_samples)
        cluster_class[clusterkey]['lastINDpos'] = last_INDpos
        #cluster_class[clusterkey]['ref_pos_alignement'] = ref_pos_alignement
        
#LAST CLUSTER
for t in cluster_to_save:
    if t[0] == 'DNA':
        cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores)
        t.append(str(cfd_score))
    else:
        cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
        t.append(str(cfd_score))


cluster_to_save.sort(key = lambda x : (float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse = True)
cluster_to_save_mmbl = cluster_to_save.copy()
cluster_to_save_mmbl.sort(key = lambda x : (int(x[8]), int(x[7])))

keys_seen = []
saved = False

if cluster_to_save:
    c = cluster_to_save[0]

    c.pop(-2)
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    if c[13] == 'n':
        cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
        cfd_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
    else:
        cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
        ref_generated = alignRefFromVar(c, dict_ref_seq)#, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
        ref_this_clus = ref_generated[0]
        reference_pam = ref_this_clus[pam_begin:pam_end]
        cfd_ref = ref_generated[1]
        found_creation = False
        for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
            if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:     #ref char not in set of general pam char
                found_creation = True
        if found_creation:
            c[11] = c[2][pam_begin:pam_end]        
        cfd_best.write('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
    cluster_to_save.pop(0)
    keys_seen.append(cfd_clus_key)
    saved = True

list_for_alt = []
for c in cluster_to_save:
    c.pop(-2)
    new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    if new_cfd_clus_key not in keys_seen:
        keys_seen.append(new_cfd_clus_key)
        if c[13] == 'n':
            cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
            #cfd_alt.write('\t'.join(c)+'\tn\n')
            list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
        else:
            cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
            ref_generated = alignRefFromVar(c, dict_ref_seq)#, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
            ref_this_clus = ref_generated[0]
            reference_pam = ref_this_clus[pam_begin:pam_end]
            cfd_ref = ref_generated[1]
            found_creation = False
            for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
                if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:     #ref char not in set of general pam char
                    found_creation = True
            if found_creation:
                c[11] = c[2][pam_begin:pam_end]
            list_for_alt.append('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
if saved:
     cfd_best.write('\t' + str(len(list_for_alt))+ '\n')

for ele in list_for_alt:
    cfd_alt.write(ele+'\t'+str(len(list_for_alt))+'\n')

keys_seen = []
saved = False
if cluster_to_save_mmbl:
    c = cluster_to_save_mmbl[0]

    #c.pop(-2)
    best_seq = c[12]
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    
    if c[12] == 'n':
        #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
        mmblg_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
    else:
        #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
        ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
        ref_this_clus = ref_generated[2]
        cfd_ref = ref_generated[-1]
        mmblg_best.write('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
    cluster_to_save_mmbl.pop(0)
    keys_seen.append(cfd_clus_key)
    saved = True
    
list_for_alt = []
for c in cluster_to_save_mmbl:
    #c.pop(-2)
    new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    if new_cfd_clus_key not in keys_seen:
        keys_seen.append(new_cfd_clus_key)
        if c[12] == 'n':
            #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
            #cfd_alt.write('\t'.join(c)+'\tn\n')
            list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
        else:
            #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
            #cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
            ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
            ref_this_clus = ref_generated[2]
            cfd_ref = ref_generated[-1]
            list_for_alt.append('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
if saved:
    mmblg_best.write('\t' + str(len(list_for_alt))+ '\n')
for ele in list_for_alt:
    mmblg_alt.write(ele+'\t'+str(len(list_for_alt))+ '\n')

cfd_best.close()
cfd_alt.close()
mmblg_best.close()
mmblg_alt.close()

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.bestCFD.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.altCFD.txt')

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.bestmmblg.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.altmmblg.txt')

#Save dataframe for graph
cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
cfd_dataframe.to_csv(outputFile + '.CFDGraph.txt', sep = '\t', index = False)
print('Done', current_chr, time.time() - start_time)

print('ANALYSIS COMPLETE IN', time.time() - global_start)
        
            
