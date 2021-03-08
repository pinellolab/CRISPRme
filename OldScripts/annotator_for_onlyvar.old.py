#!/usr/bin/env python

'''
Merge of annotator, calc_samples_faster.py and scores. ONLY FOR VAR SEARCH since no distinction between semicommon etc
Prende in input il file dei top1, ordinati per chr, e estrae i samples corrispondenti. Per ogni target, salva l'insieme dei sample in samples.all.txt, crea le combinazioni tenendo i target reali
in samples.txt, poi calcola l'annotazione corrispondente e crea il file Annotation.targets e  i vari summaries.
'''


#NOTE serve 20130606_sample_info.xlsx nella stessa cartella di questo script 
#argv1 è il file .bed con le annotazioni
#argv2 è il file .cluster.txt, che è ordinato per cromosoma. Era (03/03) il file top1 ordinato per chr
#argv3 è nome del file in output
#argv4 è directory dei dizionari
#argv5 is pamfile
#argv 6 is max allowed mms
#argv 7 is genome reference directory (Eg ../../Genomes/hg38_ref)
#argv8 is guide file
# NOTE 06/03  -> removed PAM Disruption calculation
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

SIZE_DOENCH = 10000
N_THR = 3

#Return max doench value among list of extended targets
def doenchParallel(targets, model, result):
    doench_score =  azimuth.model_comparison.predict(targets,None, None, model= model, pam_audit=False)
    doench_score = [np.around(i * 100) for i in doench_score]
    max_doench = int(max(doench_score))
    result.append(max_doench)

def doenchForIupac(sequence_doench, guide_seq): 
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
            targets_for_doench[guide_seq].append(''.join(t))
    else:
        targets_for_doench[guide_seq].append(sequence_doench)

def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) + '/mismatch_score.pkl', 'rb'))
        pam_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) +'/PAM_scores.pkl', 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

# Calculates CFD score
def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores):
    score = 1
    dna_gp = 0
    sg = sg.replace('T', 'U')
    guide_seq = guide_seq.replace('T', 'U')
    s_list = list(sg)
    guide_seq_list = list(guide_seq)
    for i, sl in enumerate(s_list):  
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            key = 'r' + guide_seq_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
            score *= mm_scores[key]
            if '-' in guide_seq_list[i]:
                dna_gp = dna_gp + 1 
    score *= pam_scores[pam]
    return score

print('ESECUZIONE DI ANNOTATION E CALC SAMPLE INSIEME')
print('TEST PER ANNOTAZIONE COMPLETA: I TARGET SENZA ANNOTAZIONE SONO SALVATI COME \"n\"')
print('SE UN  TARGET HA 1+ ANNOTAZIONI, LE SALVA IN SINGOLA UNICA RIGA')
print('RIMOZIONE DEI TARGET CHE NON HANNO SAMPLES')
print('CALCOLO SCORES')
print("READING INPUT FILES")
#Dictionaries for annotating samples

#Dict for populations
pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
all_samples = pop_file.Sample.to_list()
all_pop = pop_file.Population.to_list()
dict_sample_to_pop = dict()
for  pos, i in enumerate(all_samples):
    try:
        dict_sample_to_pop[i] = all_pop[pos]        #{'S1':'POP1', 'S2':'POP1', ...}
    except:
        dict_sample_to_pop[i] = all_pop[pos]

#Dict for superpopulation
dict_pop_to_sup = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
                    'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
                    'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
                    'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
                    'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
}
superpopulation = ['EAS', 'EUR', 'AFR', 'AMR','SAS']


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
if guide_len != 20 or 'NGG' != pam:
    with open('acfd.txt', 'w+') as result:
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
outFileSample = open(outputFile + '.samples.annotation.txt', 'w') #file with real nucleotides with associated samples and annotation
outFileSummary = open(outputFile + '.Annotation.summary.txt', 'w')  # outfile open (summary file calculated on top1file)

process = subprocess.Popen(['wc', '-l', resultsFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
total_line = int(out.decode('UTF-8').split(' ')[0])
if total_line < 2:
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

if header_len == 14:    #'Both' case : comparison variant/ref is active
    header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tPAM_gen\tVar_uniq\tSamples\tReal Guide\tAnnotation Type'
else:                   #'Var' case: PAM creation and Variant_unique not calculated
    header = '#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tSamples\tReal Guide\tAnnotation Type'

mm_pos = 7      #position of mismatch column
bulge_pos = 8
outFileSample.write(header + '\n')
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
count_sample = dict()
count_pop = dict()
count_superpop = dict()

#Create -Summary_total for a file ref.Annotation.summary.txt from the y and n values of Var_uniq column
summary_barplot_from_total = False
if 'Var_uniq' in header:
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
guides_dict = dict()
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
current_guide_chr_pos = 'no'
cluster_update = open(outputFile + '.cluster.tmp.txt', 'w+')
cluster_update.write(header + '\n')        #Write header
save_cluster_targets = True
remove_iupac = False
save_total_general_table = False
add_to_general_table = dict()   #target semicommon da aggiungere alla tab generale delle guide, metto i valori di total presi dal primo target REF nel cluster di un semicommon che esiste
last_annotation = ''    #needed when counting the ref part of a semicommon in order to not redo the annotation

next(inResult)      #Skip header
for line in inResult:
    x = line.strip().split('\t')
    guide_no_bulge = x[1].replace("-","")
    if (guide_no_bulge + x[3] + x[5]) == current_guide_chr_pos:     #Target is in current cluster, simply save the sample and annotation, discard if status is F
        if save_cluster_targets:
            if remove_iupac:
                for c in x[2]:
                    if c in iupac_code:
                        break
                else:       #no break triggered
                    cluster_update.write(line.strip() + '\t.\t' + guide_no_bulge + '\t.\n')
                continue        
            #Keep the semicommon ref total value to be added to general table
            if save_total_general_table:
                for c in x[2]:
                    if c in iupac_code:
                        break
                else: #Found first REF target in cluster of semicommon
                    add_to_general_table[guide_no_bulge][int(x[mm_pos]) +  int(x[bulge_pos])] += 1
                    #Do annotation to keep numbers consistent between images and general table
                    #conto i target generali per mm threshold
                    totalDict['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    guideDict[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1

                    #Conto per annotazione
                    for ann in last_annotation.split(','):
                        if ann == 'n':
                            break
                        guideDict[guide_no_bulge][ann][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                        totalDict[ann][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    
                    #Calculate scores
                    if do_scores and  x[0] == 'X':       #Calculate scores for reference targets
                        cfd_score = calc_cfd(x[1], x[2].upper()[:-3], x[2].upper()[-2:], mm_scores, pam_scores)
                        sum_cfd = sum_cfd + cfd_score
                        try:
                            guides_dict[x[1]] = guides_dict[x[1]] + cfd_score
                        except:
                            guides_dict[x[1]] = cfd_score

                        if x[mm_pos] == '0':    #DOENCH
                            #estraggo sequenza
                            with open('bedfile_tmp.bed', 'w+') as bedfile:
                                if x[6] == '+':
                                    bedfile.write(x[3] + '\t' + str(int(x[4]) - 4 ) + '\t' + str(int(x[4]) + 23 + 3 ))
                                else:
                                    bedfile.write(x[3] + '\t' + str(int(x[4]) - 3 ) + '\t' + str(int(x[4]) + 23 + 4 ))
                            #Extract sequence from REFERENCE
                            extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + x[3] +'.enriched.fa' ' -bed bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                            extr.wait()
                            out, err = extr.communicate()
                            out = out.decode('UTF-8')
                            if x[6] == '+':
                                sequence_doench = out.strip().split('\n')[-1].upper()
                                # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                            else:
                                sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                                # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                            
                            if x[1] not in targets_for_doench:
                                targets_for_doench[x[1]] = []
                            doenchForIupac(sequence_doench, x[1])  #Get all possible targets with iupac itertools for doench
                    save_total_general_table = False
                    
            cluster_update.write(line.strip() + '\t.\t' + guide_no_bulge + '\t.\n')  #add Sample (.) GuideNoBulge and Annotation(.). Use (.) to save space
        lines_processed +=1
        if lines_processed % (mod_tot_line) == 0:
            print('Annotation: Total progress ' + str(round(lines_processed /total_line *100, 2)) + '%')
        continue
    save_cluster_targets = True
    remove_iupac = False
    current_guide_chr_pos = guide_no_bulge + x[3] + x[5]
    if x[3] != current_chr:
        if not os.path.exists(os.path.realpath(sys.argv[4]) + '/my_dict_' + x[3] + '.json'):
            pass
        else:
            print('Done ', current_chr)
            current_chr = x[3]
            chr_name = x[3]
            with open(os.path.realpath(sys.argv[4]) + '/my_dict_' + current_chr + '.json', 'r') as f:
                start_time = time.time()
                datastore = json.load(f)
                print ('Load ' + current_chr + ' done', time.time() - start_time)
    
    pos_snp = []
    tuple_var_ref = []
    target_combination = []
    pos_snp_chr = []
    set_list = []
    target_string = x[2]
    if x[6] == '-':
        target_string = target_string[::-1]
    bulge_found = 0 
    for pos, char in enumerate(target_string):
        if char == '-':
            bulge_found = bulge_found + 1 
        if char in iupac_code:
            iupac_pos = str(int(x[4]) + pos + 1 - bulge_found)
            try:
                a = (datastore[chr_name + ',' + iupac_pos])   #NOTE se non ha samples, ritorna ;ref,var
                
                ref_char = a.split(';')[-1].split(',')[0]
                var_char = a.split(';')[-1].split(',')[1]

                if x[6] == '-':
                    ref_char = rev_comp(ref_char)
                    var_char = rev_comp(var_char)

                a = a.split(';')[0]
                pos_snp.append(pos)
                pos_snp_chr.append(iupac_pos)
                tuple_var_ref.append((var_char, ref_char))
            except Exception as e:      #NOTE this error can occure if i have an IUPAC in a target that has no vcf file
                print(e)
                print('Error at ' + line.rstrip() + ', with char ' + char + ', at pos ', iupac_pos, '. No corresponding SNP position was found in the vcf file')
                a = []
                total_error = total_error + 1
            if a:
                set_list.append(set(a.split(',')))
            else:
                set_list.append(set())
    #Get Union of all samples
    union_sample = list(set().union(*set_list))
    if union_sample:
        x.append(','.join(union_sample))
    else:
        x.append('n')
    x.append(guide_no_bulge)
    #Get all combinations to remove targets that have no haplotype 
    #Create all combinations
    for i in itertools.product(*tuple_var_ref):
        t = list(target_string)
        for p, el in enumerate(pos_snp):
            t[el] = i[p]
        target_combination.append(''.join(t))
    
    target_scomposti_salvare = []
    samples_already_assigned = set()
    false_targets = 0 
    for t in target_combination:
        set_list2 = []
        final_result = x.copy()
        for ele_pos,p in enumerate(pos_snp_chr):
            a = (datastore[chr_name + ',' + p])        
            samples = a.split(';')[0] #a[:-4] 
            
            ref = a.split(';')[-1].split(',')[0]
            var = a.split(';')[-1].split(',')[1]
            if x[6] == '-':
                ref = rev_comp(ref)
                var = rev_comp(var)
            
            if t[pos_snp[ele_pos]].upper() == var:   
                if samples:
                    set_list2.append(set(samples.split(',')))
                else:
                    set_list2.append(set())
        
        if set_list2:
            common_samples = set.intersection(*set_list2)
            common_samples = common_samples - samples_already_assigned
            samples_already_assigned = samples_already_assigned.union(common_samples)
            if common_samples:
                final_result[-2] = ','.join(common_samples)
            else:
                # final_result.append('No common samples')
                final_result = []                       #DO not save results without samples
                false_targets += 1
        else:
            # final_result.append('No samples')         #DO not save results without samples
            final_result = []
            if set_list:            #Increase false_targets on targets that have at least 1 IUPAC
                false_targets += 1
        if x[6] == '-':
            t = t[::-1]
        mm_new_t = 0
        
        if final_result:
            guide_no_pam = final_result[1][pos_beg:pos_end]    
            for position_t, char_t in enumerate(t[pos_beg:pos_end]):
                if char_t.upper() != guide_no_pam[position_t]:
                    mm_new_t += 1
            final_result[2] = t
        
            #Check for pam status
            pam_ok = True
            for pam_chr_pos, pam_chr in enumerate(t[pam_begin:pam_end]):
                if pam_chr.upper() not in iupac_code_set[pam[pam_chr_pos]]:
                    pam_ok = False

            if not pam_ok or allowed_mms < (mm_new_t - int(final_result[8])):                     
                false_targets += 1
                x[-2] = ','.join(set(x[-2].split(',')) - set(common_samples))
                continue                #Remove target since id does not respect PAM or mms constrains

            final_result[7] = str(mm_new_t - int(final_result[8]))
            final_result[9] = str(mm_new_t) #total differences between targets and guide (mismatches + bulges)
            target_scomposti_salvare.append(final_result)
    if false_targets >= len(target_combination):        #If all the scomposed targets have no sample or do not respect PAM/mm threasold, the iupac target does not really exist
        line = line.strip().split('\t')
                 #Do not do annotation because target does not exists, and do not save his cluster
        save_cluster_targets = False
        continue    #DO NOT save this target because no ref homologous and no sample combination exists
        #target does not exists in enriched, but exists in reference (semi_common), so keep only the reference targets (from his cluster)
        # if line[6] == '-':
        #     reference_semicommon = list(x[2][::-1])
        # else:
        #     reference_semicommon = list(x[2])
        # for tuple_var_ref_pos, tuple_var_ref_chars in enumerate(tuple_var_ref):
        #     reference_semicommon[pos_snp[tuple_var_ref_pos]] = tuple_var_ref_chars[1]
        # new_ref_target = ''.join(reference_semicommon)
        # if line[6] == '-':
        #     new_ref_target = new_ref_target[::-1]
        # line[2] = new_ref_target  #TODO fixare perchè il carattere ref potrebbe non essere corretto (non esiste il ref che ha lo stesso gap che nel top1)
        # x[2] = new_ref_target  
        # guide_no_pam = line[1][pos_beg:pos_end]    
        # for position_t, char_t in enumerate(new_ref_target[pos_beg:pos_end]):
        #     if char_t.upper() != guide_no_pam[position_t]:          #TODO mettere lettere minuscole
        #         mm_new_t += 1     
        x[-2] = 'n'     #Since iupac target has no scomposition, it means it has no sample associated
        # x[7] = str(mm_new_t - int(x[8]))
        # x[9] = str(mm_new_t) #total differences between targets and guide (mismatches + bulges)
        # line[7] = str(mm_new_t - int(line[8]))
        # line[9] = str(mm_new_t) #total differences between targets and guide (mismatches + bulges)
        line = '\t'.join(line)
        save_cluster_targets = True        #TODO salvare il cluster ma togliendo gli iupac
        remove_iupac = True
        x = next(inResult).strip().split('\t')   #get next target of the cluster
        while (x[1].replace('-','') + x[3] + x[5]) == current_guide_chr_pos:   #while still in same cluster
            for c in x[2]:
                if c in iupac_code:
                    break
            else:   #no break triggered in previous for --> x[2] has no iupac char
                break
            x = next(inResult).strip().split('\t')
        line = '\t'.join(x)                 #Fist target in the cluster that is REF
        x.append('n')                       #Fist target in the cluster that is REF
        x.append(x[1].replace('-',''))      #Fist target in the cluster that is REF
        tuple_var_ref = []      #Since this is now a REF target, it has no iupac --> needed to save to sample.annotation file
    if target_scomposti_salvare:        #Keep the target with lowest total and mms as representative of the IUPAC target
        target_scomposti_salvare.sort(key = lambda x : (int(x[mm_pos + 2]), int(x[mm_pos]))) #Order scomposition by total and mms values
        x[2] = target_scomposti_salvare[0][2]       #Adjust Target sequence, from IUPAC to first of scomposition
        x[mm_pos] = target_scomposti_salvare[0][mm_pos]
        x[mm_pos + 2] = target_scomposti_salvare[0][mm_pos + 2]     #Adjust IUPAC with min total and mms of his scomposition
        
    #Annotate target
    visited_pop = []
    visited_superpop = []

    #inserisco la key nel dict se non presente e creo la sua matrice
    if(guide_no_bulge not in guideDict.keys()):
        guideDict[guide_no_bulge] = {}
        guideDict[guide_no_bulge]['targets'] = {}
        guideDict[guide_no_bulge]['targets'] = [0]*10

        add_to_general_table[guide_no_bulge] = [0] * 10    # GUIDE -> [ 0 0 0 0 0 ...] valori per total (mms + bulge)

        count_unique_for_guide[guide_no_bulge] = dict()                 #NOTE count_unique means that the target have at least 1 sample
        count_unique_for_guide[guide_no_bulge]['targets'] = [0]*10

        count_sample[guide_no_bulge] = dict()
        count_pop[guide_no_bulge] = dict()
        count_superpop[guide_no_bulge] = dict()

        for item in annotationsSet:
            guideDict[guide_no_bulge][item]= {}
            guideDict[guide_no_bulge][item] = [0]*10

            count_unique_for_guide[guide_no_bulge][item] = [0]*10
    
    #conto i target generali per mm threshold
    totalDict['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
    guideDict[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1

    if summary_barplot_from_total:
        if x[-2] != 'n': 
            count_unique['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            count_unique_for_guide[guide_no_bulge]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
    
    if summary_samples:
        for sample in x[-2].split(','):
            if sample == 'n':
                continue
            #Initialization if sample, pop or superpop not in dict
            if sample not in count_sample[guide_no_bulge]:
                count_sample[guide_no_bulge][sample] = {'targets': [0]*10}
                for item in annotationsSet:
                    count_sample[guide_no_bulge][sample][item] = [0]*10
                if dict_sample_to_pop[sample] not in count_pop[guide_no_bulge]:
                    count_pop[guide_no_bulge][dict_sample_to_pop[sample]] = {'targets': [0]*10}
                    for item in annotationsSet:
                        count_pop[guide_no_bulge][dict_sample_to_pop[sample]][item] = [0]*10
                if dict_pop_to_sup[dict_sample_to_pop[sample]] not in count_superpop[guide_no_bulge]:
                    count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]] = {'targets': [0]*10}
                    for item in annotationsSet:
                        count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]][item] = [0]*10
            #Add +1 to targets
            count_sample[guide_no_bulge][sample]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            if dict_sample_to_pop[sample] not in visited_pop:
                visited_pop.append(dict_sample_to_pop[sample])
                count_pop[guide_no_bulge][dict_sample_to_pop[sample]]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            if dict_pop_to_sup[dict_sample_to_pop[sample]] not in visited_superpop:
                visited_superpop.append(dict_pop_to_sup[dict_sample_to_pop[sample]])
                count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]]['targets'][int(x[mm_pos]) + int(x[bulge_pos])] += 1
        visited_pop = []
        visited_superpop = []
    
    #faccio match su albero
    foundAnnotations = sorted(annotationsTree[int(x[4]):(int(x[4])+int(len(guide_no_bulge))+1)])
    string_annotation = []
    found_bool = False
    for found in range(0, len(foundAnnotations)):
        guide = foundAnnotations[found].data
        guideSplit = guide.split('\t')
        # print(guide, str(guideSplit[0]), str(x[3]))
        if(str(guideSplit[0]) == str(x[3])):
            found_bool = True
            #outFileTargets.write(line.rstrip() + '\t' + str(guideSplit[1]) + "\n")
            string_annotation.append(str(guideSplit[1]))
            guideDict[guide_no_bulge][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            totalDict[guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            
            if summary_barplot_from_total:
                if x[-2] != 'n':
                    count_unique[guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    count_unique_for_guide[guide_no_bulge][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
            
            if summary_samples:
                for sample in x[-2].split(','):
                    if sample == 'n':
                        continue
                    #Add +1 to annotation
                    count_sample[guide_no_bulge][sample][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    if dict_sample_to_pop[sample] not in visited_pop:
                        visited_pop.append(dict_sample_to_pop[sample])
                        count_pop[guide_no_bulge][dict_sample_to_pop[sample]][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                    if dict_pop_to_sup[dict_sample_to_pop[sample]] not in visited_superpop:
                        visited_superpop.append(dict_pop_to_sup[dict_sample_to_pop[sample]])
                        count_superpop[guide_no_bulge][dict_pop_to_sup[dict_sample_to_pop[sample]]][guideSplit[1]][int(x[mm_pos]) + int(x[bulge_pos])] += 1
                visited_pop = []
                visited_superpop = []
    if not found_bool:
        x.append('n')
        #outFileTargets.write(line.rstrip() + '\tn\n')
    else:
        x.append(','.join(string_annotation))
        #outFileTargets.write(line.rstrip() + '\t' + ','.join(string_annotation) + '\n')
    last_annotation = x[-1]
    #Save union samples + annotation
    # outFileSampleAll.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n')  

    #Save cluster
    # cluster_update.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n')
    if target_scomposti_salvare:
        cluster_update.write('\t'.join(x[:-3]) + '\t' + target_scomposti_salvare[0][-2] + '\t' + '\t'.join(x[-2:]) + '\n')  ##This line does not contain IUPAC, needed for summary by position; Adjust sample list for target scomposed
    else:
        cluster_update.write('\t'.join(x) + '\n')       #This line does not contain IUPAC, needed for summary by position
    cluster_update.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n')   #Write line with iupac (if present)
    #Save scomposed targets
    if do_scores:
        for t in target_scomposti_salvare:
            outFileSample.write('\t'.join(t) + '\t' + x[-1] + '\n')
            
            #Calc scores for scomposed targets
            if t[0] == 'X':
                cfd_score = calc_cfd(t[1], t[2].upper()[:-3], t[2].upper()[-2:], mm_scores, pam_scores)
                sum_cfd = sum_cfd + cfd_score
                try:
                    guides_dict[t[1]] = guides_dict[t[1]] + cfd_score
                except:
                    guides_dict[t[1]] = cfd_score

                if t[mm_pos] == '0':    #DOENCH
                    #estraggo sequenza
                    with open('bedfile_tmp.bed', 'w+') as bedfile:
                        if t[6] == '+':
                            bedfile.write(t[3] + '\t' + str(int(t[4]) - 4 ) + '\t' + str(int(t[4]) + 23 + 3 ))
                        else:
                            bedfile.write(t[3] + '\t' + str(int(t[4]) - 3 ) + '\t' + str(int(t[4]) + 23 + 4 ))
                        
                    extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + t[3] +'.enriched.fa' ' -bed bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                    extr.wait()
                    out, err = extr.communicate()
                    out = out.decode('UTF-8')
                    if t[6] == '+':
                        sequence_doench = out.strip().split('\n')[-1].upper()
                        # sequence_doench = sequence_doench[:4] + t[2] + sequence_doench[-3:]   #Uncomment to use sequence specific for sample
                    else:
                        sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                        # sequence_doench = sequence_doench[:4] + t[2] + sequence_doench[-3:]   #Uncomment to use sequence specific for sample
                    
                    if t[1] not in targets_for_doench:
                        targets_for_doench[t[1]] = []
                    doenchForIupac(sequence_doench, t[1])  #Get all possible targets with iupac itertools for doench
        
        if not tuple_var_ref and x[0] == 'X':       #Calculate scores for reference targets
            cfd_score = calc_cfd(x[1], x[2].upper()[:-3], x[2].upper()[-2:], mm_scores, pam_scores)
            sum_cfd = sum_cfd + cfd_score
            try:
                guides_dict[x[1]] = guides_dict[x[1]] + cfd_score
            except:
                guides_dict[x[1]] = cfd_score

            if x[mm_pos] == '0':    #DOENCH
                #estraggo sequenza
                with open('bedfile_tmp.bed', 'w+') as bedfile:
                    if x[6] == '+':
                        bedfile.write(x[3] + '\t' + str(int(x[4]) - 4 ) + '\t' + str(int(x[4]) + 23 + 3 ))
                    else:
                        bedfile.write(x[3] + '\t' + str(int(x[4]) - 3 ) + '\t' + str(int(x[4]) + 23 + 4 ))
                #Extract sequence from REFERENCE
                extr = subprocess.Popen(['bedtools getfasta -fi ' + refgenomedir + '/' + x[3] +'.enriched.fa' ' -bed bedfile_tmp.bed'], shell = True, stdout=subprocess.PIPE)  #TODO insert option for .fasta
                extr.wait()
                out, err = extr.communicate()
                out = out.decode('UTF-8')
                if x[6] == '+':
                    sequence_doench = out.strip().split('\n')[-1].upper()
                    # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                else:
                    sequence_doench = reverse_complement_table(out.strip().split('\n')[-1].upper())
                    # sequence_doench = sequence_doench[:4] + x[2] + sequence_doench[-3:]
                
                if x[1] not in targets_for_doench:
                    targets_for_doench[x[1]] = []
                doenchForIupac(sequence_doench, x[1])  #Get all possible targets with iupac itertools for doench

    else:
        for t in target_scomposti_salvare:
            outFileSample.write('\t'.join(t) + '\t' + x[-1] + '\n')

    if not tuple_var_ref:
        outFileSample.write(line.rstrip() + '\t' + '\t'.join(x[-3:]) + '\n') #Save REF target in samples.annotation, needed for sum by guide
    lines_processed +=1
    if lines_processed % (mod_tot_line) == 0:
        print('Annotation: Total progress ' + str(round(lines_processed /total_line *100, 2)) + '%')

############ SAVE SUMMARIES ############


#scorro tutto il dict total e scrivo il summary, targets e ogni annotation
outFileSummary.write("-Summary_Total\n")
outFileSummary.write('targets' + '\t'+'\t'.join(str(i) for i in totalDict['targets'])+'\n')
for elem in sorted(totalDict.keys(), key = lambda s : s.lower()):
    if elem == 'targets':
        continue
    outFileSummary.write(str(elem)+'\t'+'\t'.join(str(i) for i in totalDict[elem])+'\n')


for elem in guideDict.keys():
    outFileSummary.write("-Summary_"+str(elem)+'\n')
    outFileSummary.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[elem]['targets'])+'\n')
    for item in sorted(annotationsSet, key = lambda s : s.lower()):
        outFileSummary.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[elem][item])+'\n')

#Write summaries for samples, pop, superpop
if summary_samples:
    for guide in guideDict:
        #Save sample summary
        with open(outputFile + '.sample_annotation.' + guide +'.samples.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write sample specific counting, put [0]*10 if sample was not found
            for sample in all_samples:
                result.write('-Summary_' + sample + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_sample[guide][sample]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_sample[guide][sample][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
        
        #Save population summary
        with open(outputFile + '.sample_annotation.' + guide +'.population.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write population specific counting, put [0]*10 if sample was not found
            for population in set(all_pop):
                result.write('-Summary_' + population + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_pop[guide][population]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_pop[guide][population][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
        
        #Save superpopulation summary
        with open(outputFile + '.sample_annotation.' + guide +'.superpopulation.txt', 'w+') as result:
            result.write('-Summary_Total\n')
            result.write('targets'+'\t'+'\t'.join(str(i) for i in guideDict[guide]['targets'])+'\n')
            for item in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(str(item)+'\t'+'\t'.join(str(i) for i in guideDict[guide][item])+'\n')
            #Write superpopulation specific counting, put [0]*10 if sample was not found
            for superpop in superpopulation:
                result.write('-Summary_' + superpop + '\n')
                try:
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in count_superpop[guide][superpop]['targets']) + '\n')
                except: #Sample not found in targets
                    result.write('targets' + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')
                for item in sorted(annotationsSet, key = lambda s : s.lower()):
                    try:
                        result.write(item + '\t' + '\t'.join(str(i) for i in count_superpop[guide][superpop][item]) + '\n')
                    except:
                        result.write(item + '\t' + '\t'.join(str(i) for i in [0]*10) + '\n')


#Write sumref for barplot for targets in top1 form of var/ref search
if summary_barplot_from_total:
    with open(outputFile + '.sumref.Annotation.summary.txt', 'w+') as result:
        result.write('-Summary_Total\n')
        result.write('targets'+'\t'+'\t'.join(str(i - count_unique['targets'][pos]) for pos,i in enumerate(totalDict['targets'])) + '\n')
        for elem in sorted(annotationsSet, key = lambda s : s.lower()):
            result.write(str(elem)+'\t'+'\t'.join(str(i - count_unique[elem][pos]) for pos, i in enumerate(totalDict[elem]))+'\n')
        for guide in count_unique_for_guide:
            result.write('-Summary_' + guide + '\n')
            result.write('targets' + '\t' + '\t'.join(str(i - count_unique_for_guide[guide]['targets'][pos]) for pos,i in enumerate(guideDict[guide]['targets'])) + '\n')
            for annotation in sorted(annotationsSet, key = lambda s : s.lower()):
                result.write(annotation + '\t' + '\t'.join(str(i - count_unique_for_guide[guide][annotation][pos]) for pos, i in enumerate(guideDict[guide][annotation])) + '\n')

#SAVE SCORES#
with open( 'acfd.txt', 'w+') as res, open(sys.argv[8], 'r') as guides:
    man = multiprocessing.Manager()
    shared_doench = man.list() #list containing max doech for each thread
    guides = guides.read().strip().split('\n')
    for g in guides:
        guides_dict_doench[g] = 0
        if g not in guides_dict:
            guides_dict[g] = 0    
        if g not in targets_for_doench:
            guides_dict_doench[g] = 0
        else:
            if len (targets_for_doench[g]) > SIZE_DOENCH:
                jobs = []
                remaining_splits = (len(targets_for_doench[g])//SIZE_DOENCH) + 1
                for i in range ((len(targets_for_doench[g])//SIZE_DOENCH) + 1):
                    for thr in range (min(N_THR, remaining_splits)):
                        p = multiprocessing.Process(target = doenchParallel, args=(np.asarray(targets_for_doench[g][i*N_THR*SIZE_DOENCH + thr*SIZE_DOENCH : min( i*N_THR*SIZE_DOENCH + (thr+1)*SIZE_DOENCH,len(targets_for_doench[g]))]), model, shared_doench,) )
                        remaining_splits -= 1
                        p.start()
                        jobs.append(p)
                    for i in jobs:
                        i.join()
                
                guides_dict_doench[g] = max(shared_doench)
                shared_doench =  man.list()
            else:
                start_time = time.time()
                doench_score =  azimuth.model_comparison.predict(np.asarray(targets_for_doench[g]), None, None, model= model, pam_audit=False)
                doench_score = [np.around(i * 100) for i in doench_score]
                guides_dict_doench[g] =  int(max(doench_score))
        res.write(g + '\t' + str(guides_dict[g]) + '\t' + str(guides_dict_doench[g]) + '\n')

#Save additional values from semicommon for general guide table
with open(outputFile + '.addToGeneralTable.txt', 'w+') as add_file:
    for guide in add_to_general_table:
        add_file.write(guide + '\t' + '\t'.join([str(x) for x in add_to_general_table[guide]]) + '\n')



if total_error > 0:
    print('Skipped SNP:', total_error)
print('Annotation: Total progress 100%')
print("ANNOTATION COMPLETED IN: %s seconds" % (time.time() - start_time_total))
