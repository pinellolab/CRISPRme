#!/usr/bin/env python3
"""
Brute-force off-target ground-truth generator for CRISPRme's validate-test.

This is an exhaustive Needleman-Wunsch aligner: for a given guide (spacer+PAM,
IUPAC allowed) it scans a per-chromosome (optionally variant-enriched) FASTA on
both strands and reports every alignment within the mismatch / DNA-gap / RNA-gap
budget. Because it uses full dynamic-programming alignment (not CRISPRme's TST
index), its output is an independent ground truth for `crisprme.py validate-test`.

Output is tab-separated with columns:
    CHR  RNA  DNA  Strand  Start  END  mismatches  gaps_in_RNA  gaps_in_DNA

Reproducible generation pipeline (see test/benchmark/README.md):
    1. CRISPRme `add-variants` -> per-chromosome variant-enriched FASTA
    2. this script -> brute-force reference TSV
    3. commit the TSV + its MD5 for validate-test

Credit: the dynamic-programming brute-force checker was written by
Benjamin Vyshedskiy (https://github.com/benjaminvyshedskiy/Dynamic_checker).
Vendored into CRISPRme with attribution; only output formatting (tab-delimited,
underscore column names) was adapted to the validate-test reference format.
The alignment algorithm is unchanged.
"""
import csv
import argparse
from Bio import SeqIO

# IUPAC ambiguity codes (variant-enriched genomes encode SNVs as IUPAC bases)
iupac_dict = {
    "A": {"A"}, "C": {"C"}, "G": {"G"},
    "T": {"T", "U"}, "U": {"T", "U"},
    "R": {"A", "G"}, "Y": {"C", "T"},
    "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def iupac_match(a, b):
    return bool(iupac_dict.get(a.upper(), set()) & iupac_dict.get(b.upper(), set()))


def needleman_wunsch_all_alignments(
    RNA, DNA,
    match=0, mismatch=-1, gap=-1,
    max_DNA_gaps=1, max_RNA_gaps=1, max_mismatches=4,
):
    m, n = len(RNA), len(DNA)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[[] for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + gap
        traceback[i][0].append(("up", i - 1, 0))
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap
        traceback[0][j].append(("left", 0, j - 1))

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if iupac_match(RNA[i - 1], DNA[j - 1]):
                diag_score = dp[i - 1][j - 1] + match
            else:
                diag_score = dp[i - 1][j - 1] + mismatch
            up_score = dp[i - 1][j] + gap
            left_score = dp[i][j - 1] if i == m else dp[i][j - 1] + gap
            max_score = max(diag_score, up_score, left_score)
            dp[i][j] = max_score
            if diag_score == max_score:
                traceback[i][j].append(("diag", i - 1, j - 1))
            if up_score == max_score:
                traceback[i][j].append(("up", i - 1, j))
            if left_score == max_score:
                traceback[i][j].append(("left", i, j - 1))

    alignments = []

    def backtrack(i, j, aligned1, aligned2):
        if i == 0 and j == 0:
            RNA_GAPS = aligned1.strip("-").count("-")
            DNA_GAPS = aligned2.count("-")
            Missmatches = 0 - dp[m][n] - RNA_GAPS - DNA_GAPS
            if (aligned1[-1] != "-" and
                    RNA_GAPS <= max_RNA_gaps and
                    DNA_GAPS <= max_DNA_gaps and
                    Missmatches <= max_mismatches):
                alignments.append(
                    (aligned1[::-1], aligned2[::-1], Missmatches, RNA_GAPS, DNA_GAPS)
                )
            return
        for direction, prev_i, prev_j in traceback[i][j]:
            if direction == "diag":
                backtrack(prev_i, prev_j, aligned1 + RNA[i - 1], aligned2 + DNA[j - 1])
            elif direction == "up":
                backtrack(prev_i, prev_j, aligned1 + RNA[i - 1], aligned2 + "-")
            else:  # left
                backtrack(prev_i, prev_j, aligned1 + "-", aligned2 + DNA[j - 1])

    minscore = -max_mismatches - max_DNA_gaps - max_RNA_gaps
    if dp[m][n] >= minscore:
        backtrack(m, n, "", "")
        if alignments:
            return alignments


def reverse_complement(seq):
    complement = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
    }
    return "".join(complement.get(b, b) for b in reversed(seq.upper()))


def scan(RNA, DNA, Strand, max_dna_gaps, max_rna_gaps, max_mismatches, chrom, writer):
    for starting_pos in range(0, len(DNA) - len(RNA) + max_dna_gaps + 1):
        slice_end = min(starting_pos + len(RNA) + max_rna_gaps, len(DNA))
        DNA_Slice = DNA[starting_pos:slice_end]
        if "N" in DNA_Slice:
            continue
        aligns = needleman_wunsch_all_alignments(
            RNA, DNA_Slice,
            max_DNA_gaps=max_dna_gaps, max_RNA_gaps=max_rna_gaps,
            max_mismatches=max_mismatches,
        )
        if not aligns:
            continue
        for r_aln, d_aln, mm, rg, dg in aligns:
            actual_len = len(r_aln.rstrip("-"))
            if Strand == "+":
                start = starting_pos
                end = starting_pos + actual_len
            else:
                start = len(DNA) - starting_pos - actual_len + dg
                end = len(DNA) - starting_pos + dg
            writer.writerow(
                [chrom, r_aln[:actual_len], d_aln[:actual_len], Strand,
                 start, end, mm, rg, dg]
            )


def main():
    p = argparse.ArgumentParser(
        description="Exhaustive brute-force off-target ground truth for validate-test."
    )
    p.add_argument("--fasta", required=True, help="Per-chromosome (enriched) FASTA")
    p.add_argument("--rna", required=True, help="Guide spacer+PAM (IUPAC allowed)")
    p.add_argument("--max-mismatches", type=int, default=4)
    p.add_argument("--max-dna-gaps", type=int, default=1)
    p.add_argument("--max-rna-gaps", type=int, default=1)
    p.add_argument("--chrom", required=True, help="Chromosome name for output (e.g. chr22)")
    p.add_argument("--output", required=True, help="Output TSV path")
    args = p.parse_args()

    DNA = str(SeqIO.read(args.fasta, "fasta").seq)
    with open(args.output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            ["CHR", "RNA", "DNA", "Strand", "Start", "END",
             "mismatches", "gaps_in_RNA", "gaps_in_DNA"]
        )
        scan(args.rna, DNA, "+", args.max_dna_gaps, args.max_rna_gaps,
             args.max_mismatches, args.chrom, writer)
        scan(args.rna, reverse_complement(DNA), "-", args.max_dna_gaps,
             args.max_rna_gaps, args.max_mismatches, args.chrom, writer)


if __name__ == "__main__":
    main()
