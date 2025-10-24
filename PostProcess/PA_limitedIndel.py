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
## Needleman-Wunch global alignemnt implementation with limited gaps.  ##
## This code is a part of CRISTA:  www.crista.tau.ac.il                ##
##                                                                     ##
## For academic use only, please cite crista.tau.ac.il.                ##
## Non-commercial use                                                  ##
##                                                                     ##
## Please do not change and distribute.                                ##
##                                                                     ##
#########################################################################

import numpy as np


def align_pair(
    seqA,
    seqB,
    match_score,
    mismatch_score,
    gap_score,
    gaps_allowed,
    non_gapped_5p_len=0,
):
    """
    returns the global alignment of seqA and seqB
    non_gapped_5p_len sets the length from left that cannot have gaps
    """

    # Initialization
    alignmentA = ""
    alignmentB = ""
    scoring_matrix = init_matrix(
        seqA,
        seqB,
        match_score,
        mismatch_score,
        gap_score,
        gaps_allowed,
        non_gapped_5p_len,
    )
    j = len(seqA)
    i = len(seqB)
    c = gaps_allowed
    score = 0

    # scoring_matrix: lines represent seqB, columns represent seqA
    # meaning: decreasing i = proceeding with seqB
    #          holding j = putting gaps in AlignmentA

    while i > 0 or j > 0:
        if i > 0:
            charB = seqB[i - 1]
        if j > 0:
            charA = seqA[j - 1]

        if i > 0 and j > 0:
            sigma = two_chars_score(charA, charB, match_score, mismatch_score)
            diag_score = scoring_matrix[i - 1, j - 1, c] + sigma
            if c != 0:
                up_score = scoring_matrix[i - 1, j, c - 1] + gap_score
                left_score = scoring_matrix[i, j - 1, c - 1] + gap_score
            else:
                up_score = left_score = float("-inf")
            max_score = max(diag_score, up_score, left_score)

        else:  # have to initiate arguments
            diag_score = up_score = left_score = 0
            max_score = -1

        # check in which direction to head
        if diag_score == max_score:
            # diagonal - both sequences are aligned at position, no gap
            alignmentA = seqA[j - 1] + alignmentA
            alignmentB = seqB[i - 1] + alignmentB
            i -= 1
            j -= 1
            score += sigma

        elif j == 0 or up_score == max_score:
            # up - gap in seqA
            # base case: j==0 , seqA is completed (adding gaps to beginning), seqB not yet
            alignmentA = "-" + alignmentA
            alignmentB = seqB[i - 1] + alignmentB
            i -= 1
            c -= 1
            score += gap_score

        elif i == 0 or left_score == max_score:
            # left - gap in seqB
            # base case: i==0 , seqB id completed (adding gaps to beginning), seqA not yet
            alignmentA = seqA[j - 1] + alignmentA
            alignmentB = "-" + alignmentB
            j -= 1
            c -= 1
            score += gap_score

    return (alignmentA, alignmentB, score)


def init_matrix(
    seqA, seqB, match_score, mismatch_score, gap_score, gaps_allowed, non_gapped_5p_len
):
    """
    initiates a 3D matrix according to the global alignment function
    """

    scoring_matrix = np.full(
        (len(seqB) + 1, len(seqA) + 1, gaps_allowed + 1),
        fill_value=float("-inf"),
        dtype=float,
    )

    for c in range(0, gaps_allowed + 1):
        scoring_matrix[0, 0, c] = 0
        for l in range(max(non_gapped_5p_len, 1), c + 1):
            scoring_matrix[l, 0, c] = scoring_matrix[l - 1, 0, c] + gap_score
            scoring_matrix[0, l, c] = scoring_matrix[0, l - 1, c] + gap_score

    for i in range(1, len(seqB) + 1):
        for j in range(1, len(seqA) + 1):
            for c in range(0, gaps_allowed + 1):
                match = scoring_matrix[i - 1, j - 1, c] + two_chars_score(
                    seqB[i - 1], seqA[j - 1], match_score, mismatch_score
                )
                if c != 0 and i >= non_gapped_5p_len and j >= non_gapped_5p_len:
                    delete = scoring_matrix[i - 1, j, c - 1] + gap_score
                    insert = scoring_matrix[i, j - 1, c - 1] + gap_score
                else:
                    delete = insert = float("-inf")

                scoring_matrix[i, j, c] = max(match, insert, delete)
    return scoring_matrix


def two_chars_score(C1, C2, match_score, mismatch_score):
    """
    returns the score of 2 chars comparison
    """

    if C1 == C2:
        return match_score
    else:
        return mismatch_score
