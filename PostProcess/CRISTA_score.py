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
## cleavage sites profiling [1Ã¢â‚¬â€œ5]. It accounts for the possibility of  ##
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

from typing import List

import PA_limitedIndel as PA_script
import numpy as np

import warnings
import pickle
import random
import re
import os

warnings.filterwarnings("ignore")

# globals
RF_PICKLE_PATH = "CRISTA_predictors.pkl"
MATCH_SCORE = 1.0
MISMATCH_PENALTY = 0.0
GAP_PENALTY = -1.25
MAX_ALLOWED_GAPS = 3
EXTENSION = 3
DNA_PAIRS_THERMODYNAMICS = {
    "AA": 9.1,
    "AT": 8.6,
    "TA": 6.0,
    "CA": 5.8,
    "GT": 6.5,
    "CT": 7.8,
    "GA": 5.6,
    "CG": 11.9,
    "GC": 11.1,
    "GG": 11.0,
    "TT": 9.1,
    "TG": 5.8,
    "AC": 6.5,
    "AG": 7.8,
    "TC": 5.6,
    "CC": 11.0,
}  # Breslauer et al.
DNASHAPE_DICT_FILE = "dnaShape.pkl"
DNASHAPE_DICT = None
ACGT_REPLACEMENT = {"A": "1", "C": "2", "G": "3", "T": "4", "N": "0"}

MMS_TYPE_REPLACEMENT = {
    "0": "match",
    "-1": "indel",
    "1": "wobble",
    "2": "RR transition",
    "3": "YY transition",
    "4": "transversion",
}


def agct2numerals(st: str) -> str:
    """
    Convert a nucleotide sequence to its corresponding numeral representation.

    This function replaces each nucleotide in the input string with its mapped numeral value
    according to the ACGT_REPLACEMENT dictionary.

    Args:
        st (str): The nucleotide sequence to convert.

    Returns:
        str: The numeral representation of the nucleotide sequence.
    """

    return "".join(ACGT_REPLACEMENT[x] for x in st)


def get_avg(l: List[float]) -> float:
    """
    Compute the average value of a list of numbers.

    This function returns the arithmetic mean of the input list.

    Args:
        l (list): A list of numeric values.

    Returns:
        float: The average of the list.
    """

    return sum(l) / float(len(l))


def count_mismatches(aligned_seq1: str, aligned_seq2: str) -> int:
    """
    Count the number of mismatches between two aligned sequences, excluding gaps.

    This function compares two aligned sequences and returns the number of positions
    where the characters differ and neither character is a gap, excluding the last three positions.

    Args:
        aligned_seq1 (str): The first aligned sequence.
        aligned_seq2 (str): The second aligned sequence.

    Returns:
        int: The number of mismatches (excluding gaps and the last three positions).
    """

    ending = len(aligned_seq1) - 3
    return sum(
        int(
            aligned_seq1[i] != aligned_seq2[i]
            and aligned_seq1[i] != "-"
            and aligned_seq2[i] != "-"
        )
        for i in range(ending)
    )


def cnt_bulge(aligned_seq: str) -> int:
    """
    Count the number of bulges (gaps) in an aligned sequence.

    This function returns the number of gap characters ('-') present in the input
    sequence.

    Args:
        aligned_seq (str): The aligned sequence to analyze.

    Returns:
        int: The number of bulges (gaps) in the sequence.
    """

    return aligned_seq.count("-")


def count_consecutive_inconsistencies(aligned_seq1: str, aligned_seq2: str) -> int:
    """
    Count the number of consecutive inconsistency runs between two aligned sequences.

    This function returns the number of contiguous runs where the two sequences differ,
    excluding the last three positions.

    Args:
        aligned_seq1 (str): The first aligned sequence.
        aligned_seq2 (str): The second aligned sequence.

    Returns:
        int: The number of consecutive inconsistency runs.
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


def get_DNAshape_features(dna_seq: str):
    """
    Extract DNA shape features from a nucleotide sequence.

    This function computes minor groove width (MGW), propeller twist (ProT), roll,
    and helix twist (HelT) features for the input DNA sequence using a precomputed
    dictionary of DNA shape values.

    Args:
        dna_seq (str): The DNA sequence for which to extract shape features.

    Returns:
        dict: A dictionary containing lists of MGW, ProT, Roll, and HelT values.
    """

    global DNASHAPE_DICT
    if DNASHAPE_DICT is None:
        DNASHAPE_DICT = pickle.load(open(DNASHAPE_DICT_FILE, "rb"))
    mgw = [None]
    roll = [None]
    prot = [None]
    helt = [None]
    for i in range(2, len(dna_seq) - 2):
        current_heptamer = dna_seq[i - 2 : i + 3]
        current_heptamer = re.sub(
            "N", random.choice(["A", "C", "G", "T"]), current_heptamer
        )
        current_nucleotide = DNASHAPE_DICT[current_heptamer]
        mgw += current_nucleotide["MGW"]
        roll += current_nucleotide["Roll"]
        prot += current_nucleotide["ProT"]
        helt += current_nucleotide["HelT"]

    helt_modified = [helt[1]]
    helt_modified.extend(get_avg(helt[i : i + 2]) for i in range(2, len(helt), 2))
    roll_modified = [roll[1]]
    roll_modified.extend(get_avg(roll[i : i + 2]) for i in range(2, len(roll), 2))
    return {
        "MGW": mgw[1:],
        "ProT": prot[1:],
        "Roll": roll_modified,
        "HelT": helt_modified,
    }


def get_features(full_dna_seq, aligned_sgRNA, aligned_offtarget, pa_score):

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
    # 5' -> 3', without PAM                                             # undefined
    mismatches = [-2] * 23
    offset = 26 - len(aligned_offtarget)
    for i in range(len(aligned_offtarget) - 3):
        rna_base = aligned_sgRNA[i]
        dna_base = aligned_offtarget[i]
        # Categorization of mismatch type
        if rna_base == dna_base:  # match
            mismatches[offset + i] = 0
        elif rna_base == "-" or dna_base == "-":  # indel
            mismatches[offset + i] = -1
        elif (rna_base == "T" and dna_base == "C") or (
            rna_base == "G" and dna_base == "A"
        ):  # wobble: rG:dA, rT:dC
            mismatches[offset + i] = 1
        # R-R pairing
        elif rna_base in ["G", "A"] and dna_base in ["T", "C"]:
            mismatches[offset + i] = 2
        # Y-Y pairing
        elif rna_base in ["C", "T"] and dna_base in ["G", "A"]:
            mismatches[offset + i] = 3
        elif (rna_base == "A" and dna_base == "G") or (
            rna_base == "C" and dna_base == "T"
        ):  # other transversion
            mismatches[offset + i] = 4
    # total types mismatches
    wobble_total = mismatches.count(1)
    RR_total = mismatches.count(2)
    YY_total = mismatches.count(3)
    Tv_total = mismatches.count(4)
    # pairs of nucleotides in positions 1-5 upstream to PAM (1-2, 2-3, 3-4, 4-5)
    seed_couples = []
    seed_couples.extend(
        agct2numerals(gapless_dnaseq[-8 + i : -6 + i]) for i in range(4)
    )
    # PAM and 5'-end nucleotides
    PAM_2_first = agct2numerals(gapless_dnaseq[-2:])
    PAM_N_id = agct2numerals(gapless_dnaseq[-3])
    last_pos_nucleotide = agct2numerals(gapless_dnaseq[0])
    # mismatches and bulges - linked
    consecutive_inconsistencies_cnt = count_consecutive_inconsistencies(
        aligned_sgRNA, aligned_offtarget
    )
    avg_inconsistency_length = (
        (mms_cnt + rna_bulges + dna_bulges) / float(consecutive_inconsistencies_cnt)
        if consecutive_inconsistencies_cnt > 0
        else 0
    )
    # nucleotides occupancies in DNA target sequence
    nA = gapless_dnaseq.count("A")
    nC = gapless_dnaseq.count("C")
    nG = gapless_dnaseq.count("G")
    nT = gapless_dnaseq.count("T")
    # GC content
    extended_genomic_GC_content = (
        full_dna_seq.count("C") + full_dna_seq.count("G")
    ) / float(len(full_dna_seq))
    # five nucleotides downstream to PAM
    # the model feature for additional two nucleotides is disregarded (0) but still exists
    nucleotides_down_pam = [agct2numerals(c) for c in full_dna_seq[-3:]] + [0, 0]
    # geometry features: dna_enthalpy
    extended_dna_enthalpy = sum(
        DNA_PAIRS_THERMODYNAMICS[full_dna_seq[i - 1 : i + 1]]
        for i in range(1, len(full_dna_seq))
    )
    dna_enthalpy = sum(
        DNA_PAIRS_THERMODYNAMICS[gapless_dnaseq[i - 1 : i + 1]]
        for i in range(1, len(gapless_dnaseq))
    )
    # geometry features: DNA shape per pentamer
    dna_shape_features = get_DNAshape_features(full_dna_seq)
    features = (
        [
            pa_score,
            rna_bulges + dna_bulges,
            rna_bulges,
            dna_bulges,
            PAM_2_first,
            PAM_N_id,
            last_pos_nucleotide,
            mms_cnt,
            consecutive_inconsistencies_cnt,
            avg_inconsistency_length,
        ]
        + [
            mismatches_1_4,
            mismatches_5_8,
            mismatches_9_12,
            mismatches_13_16,
            mismatches_17_end,
        ]
        + [wobble_total, YY_total, RR_total, Tv_total]
        + seed_couples
        + [agct2numerals(gapless_dnaseq[i]) for i in range(-8, -3)]
        + [
            extended_genomic_GC_content,  # upstream_50_extension_gc, downstream_50_extension_gc,
            dna_enthalpy,
            extended_dna_enthalpy,
            nA,
            nC,
            nT,
            nG,
        ]
        + nucleotides_down_pam
        + [min(dna_shape_features["MGW"])]
        + [get_avg(dna_shape_features["HelT"])]
        + [get_avg(dna_shape_features["Roll"])]
        + [get_avg(dna_shape_features["ProT"])]
        + [
            dna_shape_features["MGW"][-3],
            dna_shape_features["HelT"][-3],
            dna_shape_features["Roll"][-3],
            dna_shape_features["ProT"][-3],
        ]
    )
    return np.array(features).reshape((1, len(features)))


def align_sequences(sgRNA, genomic_extended):
    """
    :param sgRNA: 20-nt long sgRNA
    :param genomic_extended: 23-nt target site + 3-nt at each end (total 29-nt)
    :return: aligned sgRNA, aligned target (only target, no flanking)
    """

    extended_offtarget_seq = genomic_extended[:-3]
    max_score = float("-inf")

    for i in [
        0,
        6,
        1,
        5,
        2,
        4,
        3,
    ]:  # because starting at 3 is the original - we'd prefer that
        current_dna = extended_offtarget_seq[i:-3]
        (alnA, alnB, score) = PA_script.align_pair(
            seqA=sgRNA[:-3],
            seqB=current_dna,
            match_score=MATCH_SCORE,
            mismatch_score=MISMATCH_PENALTY,
            gap_score=GAP_PENALTY,
            gaps_allowed=MAX_ALLOWED_GAPS,
        )  # regular pa

        if re.search("^\-", alnA) is None and score >= max_score:  # type: ignore
            # the target can begin a '-', it means that the last nt of the sg is not paired
            # the sg cannot begin with a '-' - in that case, a better alignment would be found (shorter DNA).
            #   However, if we first found this target and then another with a different score- we'd prefer the other
            (alignmentA, alignmentB, max_score) = (alnA, alnB, score)

    # add PAM
    aligned_sgRNA = alignmentA + sgRNA[-3:]
    aligned_offtarget = alignmentB + extended_offtarget_seq[-3:]

    # print("Aligned sgRNA:  ", aligned_sgRNA)
    # print("Aligned target: ", aligned_offtarget)
    return aligned_sgRNA, aligned_offtarget, max_score


def predict_crista_score(features_lst):
    """
    :param features_df: dataframe: first col: rna, second: dna, the rest are features
    mode: either full, nogenomic or noflanking
    :return: features df + prediction col
    """
    n_predictors = 5

    path = RF_PICKLE_PATH
    if not os.path.isfile(path):  # check for existance of predictors file
        raise FileNotFoundError(
            f"Cannot find predictors file {path}, did you unzip it?"
        )
    with open(path, "rb") as pklr:
        predictors = pickle.load(pklr)

    predictions = []
    for i in range(n_predictors):
        rf_predictor = predictors[i]
        predictions.append(rf_predictor.predict(features_lst))

    return get_avg(predictions) / 8.22


def two_chars_score(
    C1: str, C2: str, match_score: float, mismatch_score: float, gap_score: float
) -> float:
    """
    Compute the score for a pair of nucleotides based on match, mismatch, or gap.

    This function returns a score depending on whether the two characters match,
    represent a gap, or are a mismatch.

    Args:
        C1 (str): The first character (nucleotide or gap).
        C2 (str): The second character (nucleotide or gap).
        match_score (float): The score for a match.
        mismatch_score (float): The score for a mismatch.
        gap_score (float): The score for a gap.

    Returns:
        float: The score for the character pair.
    """
    if C1 == C2:  # match
        return match_score
    elif C1 == "-" or C2 == "-":  # dna/rna bulge
        return gap_score
    return mismatch_score  # mismatch


def get_alignment_score(sgRNA: str, genomic_extended: str) -> float:
    """
    Calculate the alignment score between an sgRNA and a genomic target sequence.

    This function compares each nucleotide of the sgRNA to the corresponding nucleotide
    in the genomic target (excluding the last three bases), using match, mismatch,
    and gap scores.

    Args:
        sgRNA (str): The sgRNA sequence.
        genomic_extended (str): The genomic target sequence, including flanking bases.

    Returns:
        float: The total alignment score.
    """
    offtarget_seq = genomic_extended[:-3]
    score = 0
    for i, nt in enumerate(offtarget_seq):
        score += two_chars_score(
            sgRNA[i], nt, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY
        )
    return score


def CRISTA_predict(sgseq_aligned, offseq_aligned, genomic_seq_29nt):
    sgRNA_seq = sgseq_aligned.upper()
    aligned_off_seq = offseq_aligned.upper()
    dna_seq_29nt = genomic_seq_29nt.upper()
    # print(sgRNA_seq)
    # print(aligned_off_seq)
    # print(len(sgRNA_seq), len(aligned_off_seq), len(dna_seq_29nt))
    max_score = get_alignment_score(sgRNA_seq, aligned_off_seq)
    features = get_features(
        full_dna_seq=dna_seq_29nt,
        aligned_sgRNA=sgRNA_seq,
        aligned_offtarget=aligned_off_seq,
        pa_score=max_score,
    )
    crista_features = [features[0]]
    return predict_crista_score(crista_features)


def CRISTA_predict_list(
    sgseq_aligned_list: List[str],
    offseq_aligned_list: List[str],
    genomic_seq_29nt_list: List[str],
) -> float:
    crista_features = []
    for i in range(len(sgseq_aligned_list)):
        sgRNA_seq = sgseq_aligned_list[i].upper()
        aligned_off_seq = offseq_aligned_list[i].upper()
        dna_seq_29nt = genomic_seq_29nt_list[i].upper()
        max_score = get_alignment_score(sgRNA_seq, aligned_off_seq)
        features = get_features(
            full_dna_seq=dna_seq_29nt,
            aligned_sgRNA=sgRNA_seq,
            aligned_offtarget=aligned_off_seq,
            pa_score=max_score,
        )
        crista_features.append(features[0])
    return predict_crista_score(crista_features)
