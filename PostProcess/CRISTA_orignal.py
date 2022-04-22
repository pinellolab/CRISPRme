#########################################################################
#########################################################################
##                                                                     ##
##                                                                     ##
##                              CRISTA                                 ##
##                                                                     ##
##                                                                     ##
##                 A tool for CRISPR Targets Assessment                ##
##                                                                     ##
##                             v. 1.0                                  ##
##                                                                     ##
##                                                                     ##
#########################################################################
#########################################################################
## This code is provided by Shiran Abadi                               ##
##                                                                     ##
## CRISTA is based on learning a regression model using the Random     ##
## Forest algorithm within the machine learning paradigm. CRISTA can   ##
## be used to determine the propensity of a genomic site to be cleaved ##
## by a given sgRNA. CRISTA was trained on a large dataset assembled   ##
## from published data of genome-wide unbiased methods for CRISPR-Cas9 ##
## cleavage sites profiling [1–5]. It accounts for the possibility of  ##
## bulges and incorporates a wide range of features encompassing those ##
## that are specific to the genomic content, features that define the  ##
## thermodynamics of the sgRNA, and features concerning the pairwise   ##
## similarity between the sgRNA and the genomic target. Altogether,    ##
## these form a complex model that can be used to predict the          ##
## cleavage propensity of a selected genomic site.                     ##
##                                                                     ##
## More functionalities are available at www.crista.tau.ac.il          ##
##                                                                     ##
## For academic use, please cite crista.tau.ac.il.                     ##
## Non-commercial use!                                                 ##
##                                                                     ##
## Please do not change and distribute.                                ##
##                                                                     ##
#########################################################################
#########################################################################
##                                                                     ##
##    usage: command line                                              ##
##    python CRISTA.py -s SGRNA_SEQ -d GENOMIC_SEQ                     ##
##                                                                     ##
##                                                                     ##
##     SGRNA_SEQ: sgRNA sequence of 20 bases (without PAM)             ##
##     GENOMIC_SEQ: DNA target sequence with 3 additional bases at     ##
##                   each end (total of 29 nucleotides)                ##
##                                                                     ##
#########################################################################
#########################################################################
##                                                                     ##
##    Dependencies:                                                    ##
##       python 3                                                      ##
##       numpy, sklearn, pickle, and argparse modules                  ##
##                                                                     ##
#########################################################################
# usage example
# python CRISTA.py -s CTCAGCTGAGGTTGCTGCTG -d GGCCTCAGCTGAGGTTGCTGCTGTGGAAG
#########################################################################

import argparse
import pickle
import random
import re
import numpy as np
import PA_limitedIndel as PA_script

### globals
RF_PICKLE_PATH = "CRISTA_predictors.pkl"
MATCH_SCORE = 1.0
MISMATCH_PENALTY = 0.0
GAP_PENALTY = -1.25
MAX_ALLOWED_GAPS = 3
EXTENSION = 3
DNA_PAIRS_THERMODYNAMICS = {"AA": 9.1, "AT": 8.6, "TA": 6.0, "CA": 5.8, "GT": 6.5, "CT": 7.8, "GA": 5.6, "CG": 11.9,
							"GC": 11.1, "GG": 11.0, "TT": 9.1, "TG": 5.8, "AC": 6.5, "AG": 7.8, "TC": 5.6, "CC": 11.0} #Breslauer et al.
DNASHAPE_DICT_FILE = "dnaShape.pkl"
DNASHAPE_DICT = None
ACGT_REPLACEMENT = {"A": '1', "C": '2', "G": '3', "T": '4', 'N': '0'}

MMS_TYPE_REPLACEMENT = {'0': "match", '-1': "indel", '1': "wobble", '2': "RR transition", '3': "YY transition", '4': "transversion"}


def agct2numerals(st):
	new_st = ""
	for x in st:
		new_st += ACGT_REPLACEMENT[x]
	return new_st


def get_avg(l):
	return sum(l)/float(len(l))


def count_mismatches(aligned_seq1, aligned_seq2):
	"""
	:param seq1, aligned_seq2: aligned sgRNA and genomic target (seq+PAM)
	:return:
	"""
	cnt = 0
	ending = len(aligned_seq1) - 3

	for i in range(ending):
		cnt += int(aligned_seq1[i] != aligned_seq2[i] and aligned_seq1[i] != "-" and aligned_seq2[i] != "-")
	return cnt


def cnt_bulge(aligned_seq):
	return aligned_seq.count("-")


def count_consecutive_inconsistencies(aligned_seq1, aligned_seq2):
	"""
	:param seq1, aligned_seq2: aligned sgRNA and genomic target (seq+PAM)
	:return: number of concatenated-extended mismatches and bulges
	"""
	cnt = 0
	current_cnt = 0

	for i in range(len(aligned_seq2) - 3):
		if aligned_seq2[i] != aligned_seq1[i]:
			current_cnt += 1
		else:
			cnt += current_cnt > 0
			current_cnt = 0
	return cnt


def get_DNAshape_features(dna_seq):
	"""
	:param dna_seq: sequence of nucleotides
	:return: a dictionary with scores of rigidity for Major Groove Width (MGW), ProT (Propeller-Twist), Roll, and HelT (Helical-Twist).
	The values are the scores for each pentamer/hexamer as computed by DNAshape (Zhou et al., doi:10.1093/nar/gkt437)
	across the DNA sequence
	"""

	global DNASHAPE_DICT
	if DNASHAPE_DICT is None:
		DNASHAPE_DICT = pickle.load(open(DNASHAPE_DICT_FILE, "rb"))

	mgw = [None]
	roll = [None]
	prot = [None]
	helt = [None]

	for i in range(2, len(dna_seq)-2):
		current_heptamer = dna_seq[i-2 : i+3]
		current_heptamer = re.sub("N", random.choice(["A", "C", "G", "T"]), current_heptamer)
		current_nucleotide = DNASHAPE_DICT[current_heptamer]
		mgw += current_nucleotide["MGW"]
		roll += current_nucleotide["Roll"]
		prot += current_nucleotide["ProT"]
		helt += current_nucleotide["HelT"]

	helt_modified = [helt[1]]
	for i in range(2, len(helt), 2):
		helt_modified.append(get_avg(helt[i:i+2]))
	roll_modified = [roll[1]]
	for i in range(2, len(roll), 2):
		roll_modified.append(get_avg(roll[i:i+2]))

	return {"MGW":mgw[1:], "ProT": prot[1:], "Roll": roll_modified, "HelT": helt_modified}


def get_features(full_dna_seq, aligned_sgRNA, aligned_offtarget, pa_score):
	"""
	compute CRISTA features
	"""

	# get alignment features
	mms_cnt = count_mismatches(aligned_sgRNA, aligned_offtarget)
	rna_bulges = cnt_bulge(aligned_sgRNA)
	dna_bulges = cnt_bulge(aligned_offtarget)
	gapless_dnaseq = re.sub("-", "", aligned_offtarget)

	# quartets mismatches counts
	rev_rna = (aligned_sgRNA[::-1])[3:]
	rev_dna = (aligned_offtarget[::-1])[3:]
	mismatches_1_4 = count_mismatches(rev_rna[:4], rev_dna[:4])
	mismatches_5_8 = count_mismatches(rev_rna[4:8], rev_dna[4:8])
	mismatches_9_12 = count_mismatches(rev_rna[8:12], rev_dna[8:12])
	mismatches_13_16 = count_mismatches(rev_rna[12:16], rev_dna[12:16])
	mismatches_17_end = count_mismatches(rev_rna[16:], rev_dna[16:])

	# get from alignment mismatches per position
	mismatches = [-2] * 23  # 5' -> 3', without PAM                                             # undefined
	offset = 26 - len(aligned_offtarget)
	for i in range(len(aligned_offtarget) - 3):
		rna_base = aligned_sgRNA[i]
		dna_base = aligned_offtarget[i]
		# Categorization of mismatch type
		if rna_base == dna_base:                                                                # match
			mismatches[offset + i] = 0
		elif rna_base == "-" or dna_base == "-":                                                # indel
			mismatches[offset + i] = -1
		elif (rna_base == "T" and dna_base == "C") or (rna_base == "G" and dna_base == "A"):    # wobble: rG:dA, rT:dC
			mismatches[offset + i] = 1
		elif rna_base in ["G", "A"] and dna_base in ["T", "C"]:                                 # R-R pairing
			mismatches[offset + i] = 2
		elif rna_base in ["C", "T"] and dna_base in ["G", "A"]:                                 # Y-Y pairing
			mismatches[offset + i] = 3
		elif (rna_base == "A" and dna_base == "G") or (rna_base == "C" and dna_base == "T"):    # other transversion
			mismatches[offset + i] = 4

	# total types mismatches
	wobble_total = mismatches.count(1)
	RR_total = mismatches.count(2)
	YY_total = mismatches.count(3)
	Tv_total = mismatches.count(4)

	# pairs of nucleotides in positions 1-5 upstream to PAM (1-2, 2-3, 3-4, 4-5)
	seed_couples = []
	for i in range(4):
		seed_couples.append(agct2numerals(gapless_dnaseq[-8 + i:-6 + i]))

	# PAM and 5'-end nucleotides
	PAM_2_first = agct2numerals(gapless_dnaseq[-2:])
	PAM_N_id = agct2numerals(gapless_dnaseq[-3])
	last_pos_nucleotide = agct2numerals(gapless_dnaseq[0])

	# mismatches and bulges - linked
	consecutive_inconsistencies_cnt = count_consecutive_inconsistencies(aligned_sgRNA, aligned_offtarget)
	avg_inconsistency_length = (mms_cnt + rna_bulges + dna_bulges) / float(
		consecutive_inconsistencies_cnt) if consecutive_inconsistencies_cnt > 0 else 0

	# nucleotides occupancies in DNA target sequence
	nA = gapless_dnaseq.count("A")
	nC = gapless_dnaseq.count("C")
	nG = gapless_dnaseq.count("G")
	nT = gapless_dnaseq.count("T")

	# GC content
	extended_genomic_GC_content = (full_dna_seq.count("C") + full_dna_seq.count("G")) / float(
		len(full_dna_seq))

	# five nucleotides downstream to PAM
	nucleotides_down_pam = [agct2numerals(c) for c in full_dna_seq[-3:]] + [0,0] # the model feature for additional two nucleotides is disregarded (0) but still exists

	# geometry features: dna_enthalpy
	extended_dna_enthalpy = sum([DNA_PAIRS_THERMODYNAMICS[full_dna_seq[i-1:i+1]] for i in range(1, len(full_dna_seq))])
	dna_enthalpy = sum([DNA_PAIRS_THERMODYNAMICS[gapless_dnaseq[i-1:i+1]] for i in range(1, len(gapless_dnaseq))])

	# geometry features: DNA shape per pentamer
	dna_shape_features = get_DNAshape_features(full_dna_seq)

	features = [pa_score, rna_bulges + dna_bulges, rna_bulges, dna_bulges, PAM_2_first,
	        PAM_N_id, last_pos_nucleotide, mms_cnt,
	        consecutive_inconsistencies_cnt, avg_inconsistency_length] \
	       + [mismatches_1_4, mismatches_5_8, mismatches_9_12, mismatches_13_16, mismatches_17_end] + \
	       [wobble_total, YY_total, RR_total, Tv_total] + \
	       seed_couples + [agct2numerals(gapless_dnaseq[i]) for i in range(-8, -3)] + \
	       [extended_genomic_GC_content, #upstream_50_extension_gc, downstream_50_extension_gc,
	        dna_enthalpy, extended_dna_enthalpy,
	        nA, nC, nT, nG] + nucleotides_down_pam + \
	       [min(dna_shape_features["MGW"])] + \
	       [get_avg(dna_shape_features["HelT"])] + \
	       [get_avg(dna_shape_features["Roll"])] + \
	       [get_avg(dna_shape_features["ProT"])] + \
	       [dna_shape_features["MGW"][-3], dna_shape_features["HelT"][-3],
	        dna_shape_features["Roll"][-3], dna_shape_features["ProT"][-3]]
	return np.array(features).reshape((1,len(features)))


def align_sequences(sgRNA, genomic_extended):
	"""
	:param sgRNA: 20-nt long sgRNA
	:param genomic_extended: 23-nt target site + 3-nt at each end (total 29-nt)
	:return: aligned sgRNA, aligned target (only target, no flanking)
	"""

	extended_offtarget_seq = genomic_extended[:-3]
	max_score = float("-inf")

	for i in [0, 6, 1, 5, 2, 4, 3]: #because starting at 3 is the original - we'd prefer that
		current_dna = extended_offtarget_seq[i:-3]
		(alnA, alnB, score) = PA_script.align_pair(seqA=sgRNA[:-3], seqB=current_dna, match_score=MATCH_SCORE, mismatch_score=MISMATCH_PENALTY, gap_score=GAP_PENALTY, gaps_allowed=MAX_ALLOWED_GAPS) #regular pa

		if re.search("^\-", alnA) is None and score >= max_score:
			# the target can begin a '-', it means that the last nt of the sg is not paired
			# the sg cannot begin with a '-' - in that case, a better alignment would be found (shorter DNA).
			#   However, if we first found this target and then another with a different score- we'd prefer the other
			(alignmentA, alignmentB, max_score) = (alnA, alnB, score)

	# add PAM
	aligned_sgRNA = alignmentA + sgRNA[-3:]
	aligned_offtarget = alignmentB + extended_offtarget_seq[-3:]

	print("Aligned sgRNA:  ", aligned_sgRNA)
	print("Aligned target: ", aligned_offtarget)
	return aligned_sgRNA, aligned_offtarget, max_score


def predict_crista_score(features_lst):
	"""
	:param features_df: dataframe: first col: rna, second: dna, the rest are features
	mode: either full, nogenomic or noflanking
	:return: features df + prediction col
	"""
	n_predictors = 5

	path = RF_PICKLE_PATH
	with open(path, "rb") as pklr:
		predictors = pickle.load(pklr)

	predictions = []
	for i in range(n_predictors):
		rf_predictor = predictors[i]
		predictions.append(rf_predictor.predict(features_lst))

	return get_avg(predictions) / 8.22


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='CRISTA, a tool for CRISPR Target Assessment')
	parser.add_argument('--sgseq', '-s', required=True, help="sgRNA seq of 20 bases (without PAM)")
	parser.add_argument('--genomic_seq', '-d', required=True, help="DNA target sequence with 3 additional bases at each end (total of 29 nucleotides)")

	args = parser.parse_args()
	sgRNA_seq = args.sgseq
	full_dna_seq = args.genomic_seq

	### validate input: sgrna, genomic
	try:
		sgRNA_seq_re = re.search("[acgtu]+", sgRNA_seq, re.IGNORECASE)
		assert sgRNA_seq_re is not None and len(sgRNA_seq_re.group())==20, "sgRNA sequence must be 20-nt long sequence of ACGTU nucleotides"
		full_dna_seq_re = re.search("[acgtu]+", full_dna_seq, re.IGNORECASE)
		assert full_dna_seq_re is not None and len(full_dna_seq_re.group())==29, "genomic sequence must be 29-nt long sequence of ACGTU nucleotides"
	except:
		print("Invalid arguments.")
		print(parser.parse_args(['-h']))
		exit()

	sgRNA_seq = sgRNA_seq.upper() + "NGG"
	full_dna_seq = full_dna_seq.upper()

	print("Running CRISTA")
	### align_sequences
	aligned_sgRNA, aligned_offtarget, max_score = align_sequences(sgRNA=sgRNA_seq, genomic_extended=full_dna_seq)

	### get features
	features = get_features(full_dna_seq=full_dna_seq, aligned_sgRNA=aligned_sgRNA, aligned_offtarget=aligned_offtarget,
							pa_score=max_score)
	### predict
	prediction = predict_crista_score(features)
	print("CRISTA predicted score:", prediction[0])