#!/usr/bin/env python


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
import pickle  # to read CFD matrices
from supportFunctions.loadSample import associateSample
import re
import multiprocessing
import string
# import azimuth.model_comparison
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
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

SIZE_DOENCH = 10000
N_THR = 3

# Return max doench value among list of extended targets


# def doenchParallel(targets, model, result):
#     doench_score = azimuth.model_comparison.predict(
#         targets, None, None, model=model, pam_audit=False)
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
    if do_scores == False:
        score = 0
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


def alignRefFromVar(t, refpos):
    '''
    Dato in input una lista del tipo
    t = ['BulgeType','Guide','Target','chr','pos','clusterpos'....]
    refpos = [(3,'T'),(5,'C')]
    Sostituisce i caratteri di refpos nelle posizioni indicate, calcola il CFD del nuovo reftarget e modifica i valori ricalcolando il mms value
    ritorna line con il target, il mms count aggiornati, e ID AF SNPPos messi a 'n' #TODO controllare se il reset di questi valori va bene.
    Return: ['BulgeType','Guide','TargetVersioneRef','chr','pos','clusterpos','mms aggiornati'...]
    '''
    line = t.copy()
    tmp_gap_position = []
    if line[0] == 'X':
        # Eg if max bulge dna is 1: pam end -> [' ', 'A','C','G',...,'A','G','G']; pam begin -> ['T','G','T','C','A','G',...,'A','T',' ']
        target_to_modify_scomposition = list(
            blank_add_begin * max_dna_bulges + line[2] + blank_add_end * max_dna_bulges)
    elif line[0] == 'DNA':
        target_to_modify_scomposition = list(blank_add_begin * (max_dna_bulges - int(
            line[bulge_pos])) + line[2] + blank_add_end * (max_dna_bulges - int(line[bulge_pos])))
    else:
        tmp_gap_position = [g.start() for g in re.finditer('-', line[2])]
        target_to_modify_scomposition = list(blank_add_begin * (max_dna_bulges + int(
            line[bulge_pos])) + line[2].replace('-', '') + blank_add_end * (max_dna_bulges + int(line[bulge_pos])))

    # for scom_position, pcfc in enumerate(pos_char_for_cluster):
    # target_to_modify_scomposition = target_to_modify.copy()     #Copy the datastructure of the target to substitute the iupac for each valid scomposition

    for elem in refpos:  # elem is tuple (position in target, ref char)
        # if target to modify is ' '
        if target_to_modify_scomposition[elem[0]] == ' ':
            continue
        # if target_to_modify_scomposition[elem[0]] not in iupac_code:                 #Non IUPAC char should not be targeted
        #     print('Error: Substituting into a NON IUPAC character:', line, str(elem[0]) , pos_char_for_cluster)
        target_to_modify_scomposition[elem[0]] = elem[1]

    # Fix the removed '-' in RNA targets
    for tmp_g_p in tmp_gap_position:
        target_to_modify_scomposition.insert(
            tmp_g_p + (max_dna_bulges * pam_multiplier) + int(line[bulge_pos]) * pam_multiplier, '-')

    line[2] = ''.join(target_to_modify_scomposition).strip()

    mm_new_t = 0
    tmp_pos_mms = None  # lista posizione dei mms

    guide_no_pam = line[1][pos_beg:pos_end]
    list_t = list(line[2])
    for position_t, char_t in enumerate(line[2][pos_beg:pos_end]):
        if char_t.upper() != guide_no_pam[position_t]:
            mm_new_t += 1
            tmp_pos_mms = position_t
            if guide_no_pam[position_t] != '-':
                list_t[sum_for_mms + position_t] = char_t.lower()

    # if allowed_mms < (mm_new_t - int(line[bulge_pos])):
    #     continue                #Remove target since id does not respect mms constrains
    line[2] = ''.join(list_t)

    line[mm_pos] = str(mm_new_t - int(line[bulge_pos]))
    # total differences between targets and guide (mismatches + bulges)
    line[bulge_pos + 1] = str(mm_new_t)
    # if do_print:
    #     print('last samples', last_samples[scom_position])
    #     print('sample remove', samples_to_remove)

    line[12] = 'n'
    # line[13] = last_annotation
    # last_IDs[scom_position]      #TODO controllo se il reset di queste informazioni va bene
    line[15] = 'n'
    # last_AF[scom_position]          #TODO controllo se il reset di queste informazioni va bene
    line[16] = 'n'
    # last_SNPpos[scom_position]   #TODO controllo se il reset di queste informazioni va bene
    line[17] = 'n'
    #line[18] = last_haplotype[scom_position]
    # line = '\t'.join(line)
    # cluster_update.write('\t'.join(line)  + '\n')
    # final_result = line.copy()
    # if tmp_pos_mms:
    #     final_result.append(int(tmp_pos_mms))
    # else: #NO mms found
    #     final_result.append(-1)

    # cluster_to_save.append(final_result)
    if line[0] == 'DNA':
        cfd_score = calc_cfd(line[1][int(line[bulge_pos]):], line[2].upper()[
                             int(line[bulge_pos]):-3], line[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        line[-1] = str(cfd_score)
    else:
        cfd_score = calc_cfd(line[1], line[2].upper()[
                             :-3], line[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        line[-1] = str(cfd_score)
    return line


print('ESECUZIONE DI ANNOTATION E CALC SAMPLE INSIEME')
print('TEST PER ANNOTAZIONE COMPLETA: I TARGET SENZA ANNOTAZIONE SONO SALVATI COME \"n\"')
print('SE UN  TARGET HA 1+ ANNOTAZIONI, LE SALVA IN SINGOLA UNICA RIGA')
print('RIMOZIONE DEI TARGET CHE NON HANNO SAMPLES')
print('CALCOLO SCORES')
print('SOSTITUZIONE IUPAC DI TUTTI I TARGET CON CHAR DEL TOP1SCOMPOSTO')
print("READING INPUT FILES")
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
if guide_len != 20:
    with open(outputFile + '.acfd.txt', 'w+') as result:
        result.write('NO SCORES')
        do_scores = False
if len_pam != 3:
    do_scores = False

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
inAnnotationFile = open(annotationFile, "r")  # file with annotations open
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
# totalDict = {}

start_time = time.time()

print("EXECUTING PRELIMINARY OPERATIONS")

# annotationsTree = IntervalTree()
# annotationsSet = set()
# # guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

# for line in inAnnotationFile:
#     x = line.split('\t')
#     x[3] = str(x[3]).rstrip("\n")
#     annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+str(x[3])
#     annotationsSet.add(str(x[3]))

# totalDict['targets'] = [0]*10
# for item in annotationsSet:
#     totalDict[item] = [0]*10

print("PRELIMINARY OPERATIONS COMPLETED IN: %s seconds" %
      (time.time() - start_time))

start_time = time.time()

print("EXECUTING ANNOTATION")

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

# with open(os.path.dirname(os.path.realpath(__file__)) + "/azimuth/saved_models/V3_model_nopos.pickle", 'rb') as f:
#     model = pickle.load(f)
max_doench = 0
sum_cfd = 0
cfd_scores = []


start_time_total = time.time()
lines_processed = 0
allowed_mms = int(sys.argv[6])
current_guide_chr_pos_direction = 'no'
cfd_best = open(outputFile + '.bestCFD.txt', 'w+')
cfd_best.write(header + '\tCFD\n')  # Write header

cfd_alt = open(outputFile + '.altCFD.txt', 'w+')
cfd_alt.write(header + '\tCFD\n')  # Write header

mmblg_best = open(outputFile + '.bestmmblg.txt', 'w+')
mmblg_best.write(header + '\tCFD\n')  # Write header

mmblg_alt = open(outputFile + '.altmmblg.txt', 'w+')
mmblg_alt.write(header + '\tCFD\n')  # Write header

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
cluster_to_save = []
# sum number times a cfd value appears, TODO aggiungere anche il conteggio per ogni sample
cfd_for_graph = {'ref': [0]*101, 'var': [0]*101}

# Categorie per la scomposizione, chiavi: 'BulgeType' + 'BulgeSize'
dictionary_entries = ['X0']
for i in range(int(max_dna_bulges)):
    dictionary_entries.append('DNA' + str(i+1))
for i in range(int(max_rna_bulges)):
    dictionary_entries.append('RNA' + str(i+1))


global_start = time.time()

for line in inResult:
    line = line.strip().split('\t')
    # print(line)
    guide_no_bulge = line[1].replace('-', '')
    # se ho iupac in posizioni che non c'erano nel top1, fai la scomposizione completa
    activate_cluster_scomposition = False
    # In current cluster
    if (guide_no_bulge + line[3] + line[5] + line[6]) == current_guide_chr_pos_direction:
        for pos_c, c in enumerate(line[2]):
            if c in iupac_code:
                if cluster_class[line[0] + line[8]] is None:
                    # Non ho le info per questa catergoria di target nel cluster corrente
                    activate_cluster_scomposition = True
                    break
                break
        else:  # no break triggered -> Nessun IUPAC
            #line[13] = last_annotation

            mm_new_t = 0
            guide_no_pam = line[1][pos_beg:pos_end]
            list_t = list(line[2])
            tmp_pos_mms = None  # lista posizione dei mms
            for position_t, char_t in enumerate(line[2][pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
                    tmp_pos_mms = position_t
                    if guide_no_pam[position_t] != '-':
                        list_t[sum_for_mms + position_t] = char_t.lower()
            if tmp_pos_mms:
                line.append(int(tmp_pos_mms))
            else:  # NO mms found
                line.append(-1)

            cluster_to_save.append(line)  # CHECK MMS BEFORE
            continue
        if activate_cluster_scomposition:
            # Calculate samples and all possible true scompositions
            pos_snp = []
            tuple_var_ref = []
            target_combination = []
            pos_snp_chr = []
            set_list = []
            last_samples = []
            last_IDs = []
            last_AF = []
            last_SNPpos = []
            #last_haplotype = []
            pos_char_for_cluster = []
            substitute_ref = []
            target_string = line[2]
            if line[6] == '-':
                target_string = target_string[::-1]
            bulge_found = 0
            for pos, char in enumerate(target_string):
                if char == '-':
                    bulge_found = bulge_found + 1
                if char in iupac_code:
                    iupac_pos = str(int(line[4]) + pos + 1 - bulge_found)
                    try:
                        # NOTE se non ha samples, ritorna ;ref,var
                        a = (datastore[chr_name + ',' + iupac_pos])
                    except Exception as e:  # NOTE this error can occure if i have an IUPAC in a target that has no vcf file
                        print(e)
                        print('Error at ' + ' '.join(line).rstrip() + ', with char ' + char + ', at pos ',
                              iupac_pos, '. No corresponding SNP position was found in the vcf file')
                        samples_this_position = []
                        total_error = total_error + 1
                    else:
                        a = a.split('/')
                        # Set of samples in this position (even if they have different var character)
                        samples_this_position = []
                        character_list = []  # List of var characters, later add the reference character
                        for samples_chars in a:  # samples_char can be 'sample1,sample2;A,T' or ';A,T'
                            samples_chars = samples_chars.split(';')
                            ref_char = samples_chars[1].split(',')[0]
                            var_char = samples_chars[1].split(',')[1]
                            if line[6] == '-':
                                ref_char = rev_comp(ref_char)
                                var_char = rev_comp(var_char)
                            character_list.append(var_char)
                            if samples_chars[0]:
                                samples_this_position.extend(
                                    samples_chars[0].split(','))
                        character_list.append(ref_char)
                        pos_snp.append(pos)
                        pos_snp_chr.append(iupac_pos)
                        tuple_var_ref.append(tuple(character_list))
                        substitute_ref.append(ref_char)

                    if samples_this_position:
                        set_list.append(set(samples_this_position))
                    else:
                        set_list.append(set())

            # Get all combinations to remove targets that have no haplotype
            # Create all combinations
            for i in itertools.product(*tuple_var_ref):
                t = list(target_string)
                for p, el in enumerate(pos_snp):
                    t[el] = i[p]
                target_combination.append(''.join(t))

            target_scomposti_salvare = []
            target_scomposti_non_validi = []

            samples_already_assigned = set()
            false_targets = 0

            if not tuple_var_ref:
                print('Should never activate because i have at least one IUPAC')
                target_scomposti_salvare.append(line)
            # Get PAM of reference
            if line[6] == '-':
                tmp_t = target_combination[-1][::-1]
                reference_pam = tmp_t[pam_begin:pam_end]
            else:
                reference_pam = target_combination[-1][pam_begin:pam_end]
            found_creation = False
            for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
                # ref char not in set of general pam char
                if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:
                    found_creation = True

            for t in target_combination:
                set_list2 = []
                list_rsID = []
                list_af = []
                list_snp = []
                #list_haplotype = []
                final_result = line.copy()
                for ele_pos, p in enumerate(pos_snp_chr):
                    try:
                        a = (datastore[chr_name + ',' + p])
                    except:
                        set_list2.append(set())
                    else:
                        a = a.split('/')
                        for samples_chars in a:  # samples_char can be 'sample1,sample2;A,T' or ';A,T'
                            samples_chars = samples_chars.split(';')
                            ref_char = samples_chars[1].split(',')[0]
                            var_char = samples_chars[1].split(',')[1]
                            if line[6] == '-':
                                var_char = rev_comp(var_char)

                            if t[pos_snp[ele_pos]].upper() == var_char:
                                if samples_chars[0]:
                                    set_list2.append(
                                        set(samples_chars[0].split(',')))
                                    list_rsID.append(samples_chars[2])
                                    list_af.append(samples_chars[3])
                                    list_snp.append(
                                        chr_name + '_' + p + '_' + ref_char + '_' + samples_chars[1].split(',')[1])
                                    # list_haplotype.append(var_char)
                                    # forse append ID e AF
                                else:
                                    set_list2.append(set())
                                break

                if set_list2:
                    common_samples = set.intersection(*set_list2)
                    common_samples = common_samples - samples_already_assigned
                    samples_already_assigned = samples_already_assigned.union(
                        common_samples)
                    if common_samples:
                        final_result[12] = ','.join(common_samples)
                        final_result[15] = ','.join(list_rsID)
                        final_result[16] = ','.join(list_af)
                        final_result[17] = ','.join(list_snp)
                        #final_result[18] = ','.join(list_haplotype)
                        # final_result.append(','.join(list_rsID))
                        # final_result.append(','.join(list_af))
                        # final_result.append(pos_snp_chr)
                        # se ho sample in comune allora mi salvo gli ID delle varianti della stessa posizione
                    else:
                        # final_result.append('No common samples')
                        final_result = []  # DO not save results without samples
                        false_targets += 1
                else:
                    # final_result.append('No samples')         #DO not save results without samples
                    final_result = []
                    if set_list:  # Increase false_targets on targets that have at least 1 IUPAC
                        false_targets += 1
                if line[6] == '-':
                    t = t[::-1]
                mm_new_t = 0
                tmp_pos_mms = None  # lista posizione dei mms
                if final_result:
                    guide_no_pam = final_result[1][pos_beg:pos_end]
                    list_t = list(t)
                    for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                        if char_t.upper() != guide_no_pam[position_t]:
                            mm_new_t += 1
                            tmp_pos_mms = position_t
                            if guide_no_pam[position_t] != '-':
                                list_t[sum_for_mms +
                                       position_t] = char_t.lower()
                    final_result[2] = ''.join(list_t)  # t
                    if tmp_pos_mms:
                        final_result.append(int(tmp_pos_mms))
                    else:  # NO mms found
                        final_result.append(-1)
                    # print('f_r', final_result)
                    # Check for pam status
                    pam_ok = True
                    for pam_chr_pos, pam_chr in enumerate(t[pam_begin:pam_end]):
                        if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                            pam_ok = False

                    if allowed_mms < (mm_new_t - int(final_result[8])) and pam_ok:
                        target_scomposti_non_validi.append(final_result)
                    if not pam_ok or allowed_mms < (mm_new_t - int(final_result[8])):
                        false_targets += 1
                        continue  # Remove target since id does not respect PAM or mms constrains

                    final_result[7] = str(mm_new_t - int(final_result[8]))
                    # total differences between targets and guide (mismatches + bulges)
                    final_result[9] = str(mm_new_t)
                    if found_creation:
                        final_result[10] = t[pam_begin:pam_end]
                    #final_result[13] = last_annotation
                    target_scomposti_salvare.append(final_result)

            # almeno una volta nell'iterazione t in tartarget_scomposti_salvare o  target_scomposti_non_validi,
            adjust_ref_position = True
            # aggiorno la variabile ref_pos_alignement con la posizione corretta e il char ref da sostituire
            ref_pos_alignement = []
            for t in target_scomposti_salvare:
                last_samples.append(set(t[12].split(',')))
                last_IDs.append(t[15])
                last_AF.append(t[16])
                last_SNPpos.append(t[17])
                # last_haplotype.append(t[18])
                pos_char_for_cluster.append([])
                original_targ_len = len(t[2])
                if t[0] == 'RNA':
                    # list of indices where '-' is located
                    gap_position = [g.start() for g in re.finditer('-', t[2])]
                idx_for_ref = 0
                for index_pos_snp, elem in enumerate(pos_snp):
                    add_to_count = 0
                    add_blank = int(t[bulge_pos])
                    if t[0] == 'RNA':
                        add_blank = 0
                        for i in gap_position:
                            # TODO try to remove the pam_at_beginning check to improve performances (it's done only on top1 so maybe it's ok)
                            if not pam_at_beginning:
                                if t[6] == '-':
                                    if (original_targ_len - elem - 1) < i:
                                        add_to_count += 1
                                else:
                                    if elem < i:
                                        add_to_count += 1
                            else:
                                if t[6] == '-':
                                    if (original_targ_len - elem - 1) > i:
                                        add_to_count -= 1
                                else:
                                    if elem > i:
                                        add_to_count -= 1

                    if t[6] == '-':
                        pos_char_for_cluster[-1].append(((len(t[2]) - elem - 1) + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][(len(t[2]) - elem - 1)].upper()))
                        if adjust_ref_position:
                            ref_pos_alignement.append(((len(t[2]) - elem - 1) + add_to_count + (
                                max_dna_bulges - int(add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                    else:
                        # save (position, character) of best scomposition
                        pos_char_for_cluster[-1].append((elem + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][elem].upper()))
                        if adjust_ref_position:
                            ref_pos_alignement.append((elem + add_to_count + (max_dna_bulges - int(
                                add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                    idx_for_ref += 1
                adjust_ref_position = False
                #t[13] = last_annotation

                cluster_to_save.append(t)
            for index_pos_snp, t in enumerate(target_scomposti_non_validi):
                # sample_each_targ_combination.append(set_list2_for_targetscompostinonvalidi[pos])
                last_samples.append(set(t[12].split(',')))
                last_IDs.append(t[15])
                last_AF.append(t[16])
                last_SNPpos.append(t[17])
                # last_haplotype.append(t[18])
                pos_char_for_cluster.append([])
                original_targ_len = len(t[2])
                if t[0] == 'RNA':
                    # list of indices where '-' is located
                    gap_position = [g.start() for g in re.finditer('-', t[2])]
                idx_for_ref = 0
                for elem in pos_snp:
                    add_to_count = 0
                    add_blank = int(t[bulge_pos])
                    if t[0] == 'RNA':
                        add_blank = 0
                        for i in gap_position:
                            # TODO try to remove the pam_at_beginning check to improve performances (it's done only on top1 so maybe it's ok)
                            if not pam_at_beginning:
                                if t[6] == '-':
                                    if (original_targ_len - elem - 1) < i:
                                        add_to_count += 1
                                else:
                                    if elem < i:
                                        add_to_count += 1
                            else:
                                if t[6] == '-':
                                    if (original_targ_len - elem - 1) > i:
                                        add_to_count -= 1
                                else:
                                    if elem > i:
                                        add_to_count -= 1

                    if t[6] == '-':
                        pos_char_for_cluster[-1].append(((len(t[2]) - elem - 1) + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][(len(t[2]) - elem - 1)].upper()))
                        if adjust_ref_position:
                            ref_pos_alignement.append(((len(t[2]) - elem - 1) + add_to_count + (
                                max_dna_bulges - int(add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                    else:
                        # save (position, character) of best scomposition
                        pos_char_for_cluster[-1].append((elem + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][elem].upper()))
                        if adjust_ref_position:
                            ref_pos_alignement.append((elem + add_to_count + (max_dna_bulges - int(
                                add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                    idx_for_ref += 1
                adjust_ref_position = False

            clusterkey = line[0] + line[8]
            cluster_class[clusterkey] = dict()
            cluster_class[clusterkey]['poscharforcluster'] = pos_char_for_cluster
            cluster_class[clusterkey]['lastsamples'] = last_samples
            cluster_class[clusterkey]['lastID'] = last_IDs  # ID
            # AF (come last_samples)
            cluster_class[clusterkey]['lastAF'] = last_AF
            cluster_class[clusterkey]['lastSNPpos'] = last_SNPpos
            cluster_class[clusterkey]['ref_pos_alignement'] = ref_pos_alignement
            #cluster_class[clusterkey]['lastHaplotype'] = last_haplotype

        else:
            # Use previous scompositions to create possible scomposition in cluster
            clusterkey = line[0] + line[8]
            pos_char_for_cluster = cluster_class[clusterkey]['poscharforcluster']
            last_samples = cluster_class[clusterkey]['lastsamples']
            last_IDs = cluster_class[clusterkey]['lastID']
            last_AF = cluster_class[clusterkey]['lastAF']
            last_SNPpos = cluster_class[clusterkey]['lastSNPpos']
            #last_haplotype = cluster_class[clusterkey]['lastHaplotype']
            tmp_gap_position = []
            if line[0] == 'X':
                # Eg if max bulge dna is 1: pam end -> [' ', 'A','C','G',...,'A','G','G']; pam begin -> ['T','G','T','C','A','G',...,'A','T',' ']
                target_to_modify = list(
                    blank_add_begin * max_dna_bulges + line[2] + blank_add_end * max_dna_bulges)
            elif line[0] == 'DNA':
                target_to_modify = list(blank_add_begin * (max_dna_bulges - int(
                    line[bulge_pos])) + line[2] + blank_add_end * (max_dna_bulges - int(line[bulge_pos])))
            else:
                tmp_gap_position = [g.start()
                                    for g in re.finditer('-', line[2])]
                target_to_modify = list(blank_add_begin * (max_dna_bulges + int(line[bulge_pos])) + line[2].replace(
                    '-', '') + blank_add_end * (max_dna_bulges + int(line[bulge_pos])))

            for scom_position, pcfc in enumerate(pos_char_for_cluster):
                # Copy the datastructure of the target to substitute the iupac for each valid scomposition
                target_to_modify_scomposition = target_to_modify.copy()

                # elem is tuple (position in target, character that was substituted in top1scomposed)
                for elem in pcfc:
                    # if target to modify is ' '
                    if target_to_modify_scomposition[elem[0]] == ' ':
                        continue
                    # Non IUPAC char should not be targeted
                    if target_to_modify_scomposition[elem[0]] not in iupac_code:
                        print('Error: Substituting into a NON IUPAC character:', line, str(
                            elem[0]), pos_char_for_cluster)
                    target_to_modify_scomposition[elem[0]] = elem[1]

                # Fix the removed '-' in RNA targets
                for tmp_g_p in tmp_gap_position:
                    target_to_modify_scomposition.insert(
                        tmp_g_p + (max_dna_bulges * pam_multiplier) + int(line[bulge_pos]) * pam_multiplier, '-')

                line[2] = ''.join(target_to_modify_scomposition).strip()

                mm_new_t = 0
                tmp_pos_mms = None  # lista posizione dei mms
                guide_no_pam = line[1][pos_beg:pos_end]
                list_t = list(line[2])
                for position_t, char_t in enumerate(line[2][pos_beg:pos_end]):
                    if char_t.upper() != guide_no_pam[position_t]:
                        mm_new_t += 1
                        tmp_pos_mms = position_t
                        if guide_no_pam[position_t] != '-':
                            list_t[sum_for_mms + position_t] = char_t.lower()

                if allowed_mms < (mm_new_t - int(line[bulge_pos])):
                    continue  # Remove target since id does not respect mms constrains
                # TODO add pam creation check
                line[2] = ''.join(list_t)

                line[mm_pos] = str(mm_new_t - int(line[bulge_pos]))
                # total differences between targets and guide (mismatches + bulges)
                line[bulge_pos + 1] = str(mm_new_t)
                # if do_print:
                #     print('last samples', last_samples[scom_position])
                #     print('sample remove', samples_to_remove)

                line[12] = ','.join(last_samples[scom_position])
                #line[13] = last_annotation
                line[15] = last_IDs[scom_position]
                line[16] = last_AF[scom_position]
                line[17] = last_SNPpos[scom_position]
                #line[18] = last_haplotype[scom_position]
                # line = '\t'.join(line)
                # cluster_update.write('\t'.join(line)  + '\n')
                final_result = line.copy()
                if tmp_pos_mms:
                    final_result.append(int(tmp_pos_mms))
                else:  # NO mms found
                    final_result.append(-1)

                cluster_to_save.append(final_result)

    else:  # New cluster
        # New cluster, order previous by CFD, and save the highest CFD

        # check_position = 0

        # Calcolo lo score per ogni target nel cluster
        for t in cluster_to_save:

            if t[0] == 'DNA':
                cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(
                    t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
                t.append(str(cfd_score))
            else:
                cfd_score = calc_cfd(t[1], t[2].upper()[
                                     :-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
                t.append(str(cfd_score))

        cluster_to_save.sort(key=lambda x: (
            float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse=True)

        cluster_to_save_mmbl = cluster_to_save.copy()
        cluster_to_save_mmbl.sort(key=lambda x: (int(x[8]), int(x[7])))

        keys_seen = []
        saved = False
        """
        ref_this_clus = 'n'
        cfd_ref = None
        for c in cluster_to_save:
            if c[12] != 'n':
                ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                ref_this_clus = ref_generated[2]
                cfd_ref = ref_generated[-1]
                break
        """
        if cluster_to_save:
            c = cluster_to_save[0]

            # if c[12] != 'n':        #Calcola ref nello stesso allineamento del var
            #    ref_best_var = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement']) #TODO salvare questa line sul file

            c.pop(-2)
            best_seq = c[12]
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]

            if c[12] == 'n':
                cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                cfd_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
            else:
                cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                ref_generated = alignRefFromVar(
                    c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                ref_this_clus = ref_generated[2]
                cfd_ref = ref_generated[-1]
                cfd_best.write('\t'.join(c)+'\t' +
                               ref_this_clus+'\t'+str(cfd_ref))
            cluster_to_save.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True

        list_for_alt = []
        for c in cluster_to_save:
            c.pop(-2)
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                if c[12] == 'n':
                    cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                    # cfd_alt.write('\t'.join(c)+'\tn\n')
                    list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
                else:
                    cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                    # cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
                    ref_generated = alignRefFromVar(
                        c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                    ref_this_clus = ref_generated[2]
                    cfd_ref = ref_generated[-1]
                    list_for_alt.append(
                        '\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
        if saved:
            cfd_best.write('\t' + str(len(list_for_alt)) + '\n')
        for ele in list_for_alt:
            cfd_alt.write(ele+'\t'+str(len(list_for_alt)) + '\n')

        keys_seen = []
        saved = False
        if cluster_to_save_mmbl:
            c = cluster_to_save_mmbl[0]

            # if c[12] != 'n':        #Calcola ref nello stesso allineamento del var
            #    ref_best_var = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement']) #TODO salvare questa line sul file

            # c.pop(-2)
            best_seq = c[12]
            cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]

            if c[12] == 'n':
                #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                mmblg_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
            else:
                #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                ref_generated = alignRefFromVar(
                    c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                ref_this_clus = ref_generated[2]
                cfd_ref = ref_generated[-1]
                mmblg_best.write('\t'.join(c)+'\t' +
                                 ref_this_clus+'\t'+str(cfd_ref))
            cluster_to_save_mmbl.pop(0)
            keys_seen.append(cfd_clus_key)
            saved = True

        list_for_alt = []
        for c in cluster_to_save_mmbl:
            # c.pop(-2)
            new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
            if new_cfd_clus_key not in keys_seen:
                keys_seen.append(new_cfd_clus_key)
                if c[12] == 'n':
                    #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
                    # cfd_alt.write('\t'.join(c)+'\tn\n')
                    list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
                else:
                    #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
                    # cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
                    ref_generated = alignRefFromVar(
                        c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
                    ref_this_clus = ref_generated[2]
                    cfd_ref = ref_generated[-1]
                    list_for_alt.append(
                        '\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
        if saved:
            mmblg_best.write('\t' + str(len(list_for_alt)) + '\n')
        for ele in list_for_alt:
            mmblg_alt.write(ele+'\t'+str(len(list_for_alt)) + '\n')

        cluster_to_save = []
        # Dizionario dove salvare le variabile per ogni tipo di scomposizione (X,DNA1,DNA2,RNA1 etc...)
        cluster_class = dict.fromkeys(dictionary_entries)
        # Viene re-inizializzato per ogni cluster
        current_guide_chr_pos_direction = guide_no_bulge + \
            line[3] + line[5] + line[6]
        if line[3] != current_chr:  # New chr section, load a new dictionary
            if not os.path.exists(os.path.realpath(sys.argv[4]) + '/my_dict_' + line[3] + '.json'):
                #raise Exception("Dictionary not found"+line[3])
                print("DICT NOT FOUND")
                current_chr = line[3]
                chr_name = line[3]
                pass
            else:
                if current_chr != 'none':
                    print('Done', current_chr, time.time() - start_time)
                current_chr = line[3]
                chr_name = line[3]
                with open(os.path.realpath(sys.argv[4]) + '/my_dict_' + current_chr + '.json', 'r') as f:
                    datastore = json.load(f)
                    print('Analysis of ' + current_chr)
                start_time = time.time()

        # Calculate samples and all possible true scompositions
        pos_snp = []        # lista delle posizioni degli iupac nel target
        # array che contiene le tuple composte da (var1,var2..., ref) (al momento abbiamo solo var1 nei dizionari) per ogni iupac nel target
        tuple_var_ref = []
        # top1_ref = False
        # lista di tutte le possibili combinazioni di scomposizione di IUPAC
        target_combination = []
        # lista delle posizioni degli iupac nel target riferite alla posizione nel cromosoma
        pos_snp_chr = []
        set_list_top1 = []      # lista di set di samples associati ad ogni posizione IUPAC
        last_samples = []
        last_IDs = []
        last_AF = []
        last_SNPpos = []
        substitute_ref = []  # List of ref_char used to create the aligned ref target from a var
        #last_haplotype = []
        pos_char_for_cluster = []  # lista di liste. Ogni lista corrisponde ad una scomposizione. Per ogni scomposizione mi salvo la posizione nel target e il carattere che ho trovato nella scomposizione
        # andrò poi a sostituire questi caratteri negli altri target del cluster che sono nella stessa categoria
        target_string = line[2]
        if line[6] == '-':
            target_string = target_string[::-1]
        bulge_found = 0

        # Prendi i sample associati ad ogni posizione iupac, insieme ai caratteri scomposti
        for pos, char in enumerate(target_string):
            if char == '-':
                bulge_found = bulge_found + 1
            if char in iupac_code:
                iupac_pos = str(int(line[4]) + pos + 1 - bulge_found)
                try:
                    # NOTE se non ha samples, ritorna ;ref,var
                    a = (datastore[chr_name + ',' + iupac_pos])
                except Exception as e:  # NOTE this error can occure if i have an IUPAC in a target that has no vcf file
                    print(e)
                    print('Error at ' + ' '.join(line).rstrip() + ', with char ' + char + ', at pos ',
                          iupac_pos, '. No corresponding SNP position was found in the vcf file')
                    samples_this_position = []
                    total_error = total_error + 1
                else:
                    a = a.split('/')
                    # Set of samples in this position (even if they have different var character)
                    samples_this_position = []
                    character_list = []  # List of var characters, later add the reference character
                    for samples_chars in a:  # samples_char can be 'sample1,sample2;A,T' or ';A,T'
                        samples_chars = samples_chars.split(';')
                        ref_char = samples_chars[1].split(',')[0]
                        var_char = samples_chars[1].split(',')[1]
                        if line[6] == '-':
                            ref_char = rev_comp(ref_char)
                            var_char = rev_comp(var_char)
                        character_list.append(var_char)
                        if samples_chars[0]:
                            samples_this_position.extend(
                                samples_chars[0].split(','))
                    character_list.append(ref_char)
                    pos_snp.append(pos)
                    pos_snp_chr.append(iupac_pos)
                    tuple_var_ref.append(tuple(character_list))
                    # when i found a iupac, save the position and the reference char
                    substitute_ref.append(ref_char)
                if samples_this_position:
                    # posso usarla per prendere i sample che sono associati ad una determinata posizione iupac, se questo iupac non c'è in un target del
                    set_list_top1.append(set(samples_this_position))
                    # cluster, da last_sample tolgo i sample di quella posizione
                    # var char da sostituire con elenco di var char se sample hanno var diversi
                else:
                    set_list_top1.append(set())
                # print('samples this pos', samples_this_position)

        # Get all combinations to remove targets that have no haplotype
        # Create all combinations
        for i in itertools.product(*tuple_var_ref):
            t = list(target_string)
            for p, el in enumerate(pos_snp):
                t[el] = i[p]
            target_combination.append(''.join(t))

        target_scomposti_salvare = []  # lista di target che esistono realmente
        # target che non esistono per colpa del mms threshold, ma la scomposizione potrebbe essere ok per un altro traget nel cluster
        target_scomposti_non_validi = []
        samples_already_assigned = set()
        false_targets = 0

        # Get PAM of reference per vedere se ho una PAM creation. Nota che per come è stato costruito, il target reference è sempre l'ultimo di target combination
        if line[6] == '-':
            tmp_t = target_combination[-1][::-1]
            reference_pam = tmp_t[pam_begin:pam_end]
        else:
            reference_pam = target_combination[-1][pam_begin:pam_end]
        found_creation = False
        for pos_pam_ref, pam_char_ref in enumerate(reference_pam):
            # ref char not in set of general pam char
            if not iupac_code_set[pam[pos_pam_ref]] & iupac_code_set[pam_char_ref]:
                found_creation = True

        # sample_each_targ_combination = [] #lista di liste, una lista per ogni target combination valida. Dentro ogni lista c'è una lista di set
                # di sample per ogni posizione iupac
                # EG [[{s1,s2}, {s3},{s1,s4}], [...]] la prima combinazione ha 3 iupac, e quelli sono i sample associati
        # set_list2_for_targetscompostinonvalidi = []     #come sample_each_targ_comb, ma per le scomposizioni che non rispettano il mms threshold

        # Per ogni possibile scomposizione, salva solo quelle che esistono (ovvero hanno sample associati e rispettano il mms threshold e la pam)
        for t in target_combination:
            # print(t)
            set_list2 = []
            list_rsID = []
            list_af = []
            list_snp = []
            #list_haplotype = []
            final_result = line.copy()
            for ele_pos, p in enumerate(pos_snp_chr):
                try:
                    a = (datastore[chr_name + ',' + p])
                except:
                    set_list2.append(set())
                else:
                    a = a.split('/')
                    for samples_chars in a:  # samples_char can be 'sample1,sample2;A,T' or ';A,T'
                        samples_chars = samples_chars.split(';')
                        ref_char = samples_chars[1].split(',')[0]
                        var_char = samples_chars[1].split(',')[1]
                        if line[6] == '-':
                            var_char = rev_comp(var_char)

                        if t[pos_snp[ele_pos]].upper() == var_char:
                            if samples_chars[0]:
                                set_list2.append(
                                    set(samples_chars[0].split(',')))
                                list_rsID.append(samples_chars[2])
                                list_af.append(samples_chars[3])
                                list_snp.append(
                                    chr_name + '_' + p + '_' + ref_char + '_' + samples_chars[1].split(',')[1])
                                # list_haplotype.append(var_char)
                            else:
                                set_list2.append(set())
                            break

            if set_list2:
                common_samples = set.intersection(*set_list2)
                common_samples = common_samples - samples_already_assigned
                samples_already_assigned = samples_already_assigned.union(
                    common_samples)
                if common_samples:
                    final_result[12] = ','.join(common_samples)
                    final_result[15] = ','.join(list_rsID)
                    final_result[16] = ','.join(list_af)
                    final_result[17] = ','.join(list_snp)
                    #final_result[18] = ','.join(list_haplotype)
                    # final_result.append(','.join(list_rsID))
                    # final_result.append(','.join(list_af))
                    # final_result.append(pos_snp_chr)
                else:
                    # final_result.append('No common samples')
                    final_result = []  # DO not save results without common samples
                    false_targets += 1
                    # print('no common samples')
            else:
                # final_result.append('No samples')         #DO not save results without samples
                final_result = []
                # print('no samples')
                if set_list_top1:  # Increase false_targets on targets that have at least 1 IUPAC
                    false_targets += 1
            if line[6] == '-':
                t = t[::-1]
            mm_new_t = 0
            tmp_pos_mms = None  # lista posizione dei mms

            if final_result:  # Se ho samples associati
                # Calcola nuovo valori di mms
                guide_no_pam = final_result[1][pos_beg:pos_end]
                list_t = list(t)
                #print("target",t, line[5])
                # print("guide",guide_no_pam)
                for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                    if char_t.upper() != guide_no_pam[position_t]:
                        mm_new_t += 1
                        tmp_pos_mms = position_t
                        if guide_no_pam[position_t] != '-':
                            list_t[sum_for_mms + position_t] = char_t.lower()
                final_result[2] = ''.join(list_t)  # t
                if tmp_pos_mms:
                    final_result.append(int(tmp_pos_mms))
                else:  # NO mms found
                    final_result.append(-1)

                # Check for pam status
                pam_ok = True
                for pam_chr_pos, pam_chr in enumerate(t[pam_begin:pam_end]):
                    if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                        pam_ok = False
                # Se la scomposizione non rispetta il mms threshold, ma la pam è ok, magari può esistere un target nel cluster a cui questa scomposizione può andare bene
                if allowed_mms < (mm_new_t - int(final_result[8])) and pam_ok:
                    # set_list2_for_targetscompostinonvalidi.append(set_list2)
                    target_scomposti_non_validi.append(final_result)
                if not pam_ok or allowed_mms < (mm_new_t - int(final_result[8])):
                    false_targets += 1
                    continue  # Remove target since id does not respect PAM or mms constrains

                final_result[7] = str(mm_new_t - int(final_result[8]))
                # total differences between targets and guide (mismatches + bulges)
                final_result[9] = str(mm_new_t)
                if found_creation:
                    final_result[10] = t[pam_begin:pam_end]
                target_scomposti_salvare.append(final_result)

        # Calcolo annotazioni
        # foundAnnotations = sorted(
        #     annotationsTree[int(line[4]):(int(line[4])+int(len(guide_no_bulge))+1)])
        # string_annotation = []
        # found_bool = False
        # for found in range(0, len(foundAnnotations)):
        #     guide = foundAnnotations[found].data
        #     guideSplit = guide.split('\t')
        #     if(str(guideSplit[0]) == str(line[3])):
        #         found_bool = True
        #         string_annotation.append(str(guideSplit[1]))
        # if not found_bool:
        #     last_annotation = 'n'
        # else:
        #     last_annotation = ','.join(string_annotation)

        # almeno una volta nell'iterazione t in tartarget_scomposti_salvare o  target_scomposti_non_validi,
        adjust_ref_position = True
        # aggiorno la variabile ref_pos_alignement con la posizione corretta e il char ref da sostituire
        ref_pos_alignement = []  # list of (pos, ref) for this cluster class
        # analisi dei target esistenti per evitare poi fare tutte le scomposizioni se un target nel cluster è in questa stessa categoria
        for t in target_scomposti_salvare:
            last_samples.append(set(t[12].split(',')))
            last_IDs.append(t[15])
            last_AF.append(t[16])
            last_SNPpos.append(t[17])
            # last_haplotype.append(t[18])
            pos_char_for_cluster.append([])
            original_targ_len = len(t[2])
            if t[0] == 'RNA':
                # list of indices where '-' is located
                gap_position = [g.start() for g in re.finditer('-', t[2])]
            idx_for_ref = 0
            for index_pos_snp, elem in enumerate(pos_snp):
                add_to_count = 0
                add_blank = int(t[bulge_pos])
                if t[0] == 'RNA':
                    add_blank = 0
                    for i in gap_position:
                        # TODO try to remove the pam_at_beginning check to improve performances (it's done only on top1 so maybe it's ok)
                        if not pam_at_beginning:
                            if t[6] == '-':
                                if (original_targ_len - elem - 1) < i:
                                    add_to_count += 1
                            else:
                                if elem < i:
                                    add_to_count += 1
                        else:
                            if t[6] == '-':
                                if (original_targ_len - elem - 1) > i:
                                    add_to_count -= 1
                            else:
                                if elem > i:
                                    add_to_count -= 1

                if t[6] == '-':
                    pos_char_for_cluster[-1].append(((len(t[2]) - elem - 1) + add_to_count + (
                        max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][(len(t[2]) - elem - 1)].upper()))
                    if adjust_ref_position:
                        ref_pos_alignement.append(((len(t[2]) - elem - 1) + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, substitute_ref[index_pos_snp]))
                else:
                    # save (position, character) of best scomposition
                    pos_char_for_cluster[-1].append((elem + add_to_count + (
                        max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][elem].upper()))
                    if adjust_ref_position:
                        ref_pos_alignement.append((elem + add_to_count + (max_dna_bulges - int(
                            add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                idx_for_ref += 1
            adjust_ref_position = False
            #t[13] = last_annotation
            cluster_to_save.append(t)
        # Stessa cosa ma per i target che non hanno rispettato il mms threshold
        for index_pos_snp, t in enumerate(target_scomposti_non_validi):
            # sample_each_targ_combination.append(set_list2_for_targetscompostinonvalidi[pos])
            last_samples.append(set(t[12].split(',')))
            last_IDs.append(t[15])
            last_AF.append(t[16])
            last_SNPpos.append(t[17])
            # last_haplotype.append(t[18])
            pos_char_for_cluster.append([])
            original_targ_len = len(t[2])
            if t[0] == 'RNA':
                # list of indices where '-' is located
                gap_position = [g.start() for g in re.finditer('-', t[2])]
            idx_for_ref = 0
            for elem in pos_snp:
                add_to_count = 0
                add_blank = int(t[bulge_pos])
                if t[0] == 'RNA':
                    add_blank = 0
                    for i in gap_position:
                        # TODO try to remove the pam_at_beginning check to improve performances (it's done only on top1 so maybe it's ok)
                        if not pam_at_beginning:
                            if t[6] == '-':
                                if (original_targ_len - elem - 1) < i:
                                    add_to_count += 1
                            else:
                                if elem < i:
                                    add_to_count += 1
                        else:
                            if t[6] == '-':
                                if (original_targ_len - elem - 1) > i:
                                    add_to_count -= 1
                            else:
                                if elem > i:
                                    add_to_count -= 1

                if t[6] == '-':
                    pos_char_for_cluster[-1].append(((len(t[2]) - elem - 1) + add_to_count + (
                        max_dna_bulges - int(add_blank)) * pam_multiplier, t[2][(len(t[2]) - elem - 1)].upper()))
                    if adjust_ref_position:
                        ref_pos_alignement.append(((len(t[2]) - elem - 1) + add_to_count + (
                            max_dna_bulges - int(add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                else:
                    pos_char_for_cluster[-1].append((elem + add_to_count + (max_dna_bulges - int(
                        add_blank)) * pam_multiplier, t[2][elem].upper()))  # save (position, character) of scomposition
                    if adjust_ref_position:
                        ref_pos_alignement.append((elem + add_to_count + (max_dna_bulges - int(
                            add_blank)) * pam_multiplier, substitute_ref[idx_for_ref]))
                idx_for_ref += 1
            adjust_ref_position = False
        if not tuple_var_ref:  # Se il target è un reference, tutta l'analisi di prima non si fa
            # target_scomposti_salvare.append(line)
            #line[13] = last_annotation

            mm_new_t = 0
            tmp_pos_mms = None  # lista posizione dei mm
            guide_no_pam = line[1][pos_beg:pos_end]

            for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
                    tmp_pos_mms = position_t
                    # pos mm da salvare

            if tmp_pos_mms:
                line.append(int(tmp_pos_mms))
            else:  # NO mms found
                line.append(-1)
            cluster_to_save.append(line)
            # top1_ref = True
        else:       # at least one iupac
            clusterkey = line[0] + line[8]
            cluster_class[clusterkey] = dict()
            # aggiorno la categoria del target per i prossimi target
            cluster_class[clusterkey]['poscharforcluster'] = pos_char_for_cluster
            cluster_class[clusterkey]['lastsamples'] = last_samples
            cluster_class[clusterkey]['lastID'] = last_IDs  # ID
            # AF (come last_samples)
            cluster_class[clusterkey]['lastAF'] = last_AF
            cluster_class[clusterkey]['lastSNPpos'] = last_SNPpos
            cluster_class[clusterkey]['ref_pos_alignement'] = ref_pos_alignement
            # print(ref_pos_alignement)
            #cluster_class[clusterkey]['lastHaplotype'] = last_haplotype

# LAST CLUSTER
for t in cluster_to_save:
    if t[0] == 'DNA':
        cfd_score = calc_cfd(t[1][int(t[bulge_pos]):], t[2].upper()[int(
            t[bulge_pos]):-3], t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        t.append(str(cfd_score))
    else:
        cfd_score = calc_cfd(t[1], t[2].upper()[:-3],
                             t[2].upper()[-2:], mm_scores, pam_scores, do_scores)
        t.append(str(cfd_score))


cluster_to_save.sort(key=lambda x: (
    float(x[-1]), reversor(int(x[9])), reversor(int(x[-2]))), reverse=True)

cluster_to_save_mmbl = cluster_to_save.copy()
cluster_to_save_mmbl.sort(key=lambda x: (int(x[8]), int(x[7])))

keys_seen = []
saved = False
"""
ref_this_clus = 'n'
for c in cluster_to_save:
    if c[12] != 'n':
        ref_generated = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
        ref_this_clus = ref_generated[2]
        cfd_ref = ref_generated[-1]
        list_for_alt.append('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
        break
"""
if cluster_to_save:
    c = cluster_to_save[0]

    # if c[12] != 'n':        #Calcola ref nello stesso allineamento del var
    #    ref_best_var = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])   #TODO salvare questa line sul file

    c.pop(-2)
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    if c[12] == 'n':
        cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
        cfd_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
    else:
        cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
        ref_generated = alignRefFromVar(
            c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
        ref_this_clus = ref_generated[2]
        cfd_ref = ref_generated[-1]
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
        if c[12] == 'n':
            cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
            # cfd_alt.write('\t'.join(c)+'\tn\n')
            list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
        else:
            cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
            ref_generated = alignRefFromVar(
                c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
            ref_this_clus = ref_generated[2]
            cfd_ref = ref_generated[-1]
            list_for_alt.append('\t'.join(c)+'\t' +
                                ref_this_clus+'\t'+str(cfd_ref))
if saved:
    cfd_best.write('\t' + str(len(list_for_alt)) + '\n')

for ele in list_for_alt:
    cfd_alt.write(ele+'\t'+str(len(list_for_alt))+'\n')

keys_seen = []
saved = False
if cluster_to_save_mmbl:
    c = cluster_to_save_mmbl[0]

    # if c[12] != 'n':        #Calcola ref nello stesso allineamento del var
    #    ref_best_var = alignRefFromVar(c, cluster_class[c[0] + c[8]]['ref_pos_alignement']) #TODO salvare questa line sul file

    # c.pop(-2)
    best_seq = c[12]
    cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]

    if c[12] == 'n':
        #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
        mmblg_best.write('\t'.join(c)+'\tn\t'+str(c[-1]))
    else:
        #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
        ref_generated = alignRefFromVar(
            c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
        ref_this_clus = ref_generated[2]
        cfd_ref = ref_generated[-1]
        mmblg_best.write('\t'.join(c)+'\t'+ref_this_clus+'\t'+str(cfd_ref))
    cluster_to_save_mmbl.pop(0)
    keys_seen.append(cfd_clus_key)
    saved = True

list_for_alt = []
for c in cluster_to_save_mmbl:
    # c.pop(-2)
    new_cfd_clus_key = c[3] + " " + c[5] + " " + c[6] + " " + c[-2]
    if new_cfd_clus_key not in keys_seen:
        keys_seen.append(new_cfd_clus_key)
        if c[12] == 'n':
            #cfd_for_graph['ref'][round(float(c[-1]) * 100)] += 1
            # cfd_alt.write('\t'.join(c)+'\tn\n')
            list_for_alt.append('\t'.join(c)+'\tn\t'+str(c[-1]))
        else:
            #cfd_for_graph['var'][round(float(c[-1]) * 100)] += 1
            # cfd_alt.write('\t'.join(c)+'\t'+ref_this_clus+'\n')
            ref_generated = alignRefFromVar(
                c, cluster_class[c[0] + c[8]]['ref_pos_alignement'])
            ref_this_clus = ref_generated[2]
            cfd_ref = ref_generated[-1]
            list_for_alt.append('\t'.join(c)+'\t' +
                                ref_this_clus+'\t'+str(cfd_ref))
if saved:
    mmblg_best.write('\t' + str(len(list_for_alt)) + '\n')
for ele in list_for_alt:
    mmblg_alt.write(ele+'\t'+str(len(list_for_alt)) + '\n')


cfd_best.close()
cfd_alt.close()
mmblg_best.close()
mmblg_alt.close()

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.bestCFD.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.altCFD.txt')

os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.bestmmblg.txt')
os.system("sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\tCFD\tReference\tCFD_ref\t#Seq_in_cluster/' "+outputFile + '.altmmblg.txt')


# Save dataframe for graph
cfd_dataframe = pd.DataFrame.from_dict(cfd_for_graph)
cfd_dataframe.to_csv(outputFile + '.CFDGraph.txt', sep='\t', index=False)
print('Done', current_chr, time.time() - start_time)

print('ANALYSIS COMPLETE IN', time.time() - global_start)
