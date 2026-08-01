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
The alignment algorithm is vendored from that repository's
`brute_force_works_notcommandline.py` (the script that generated the committed
Cas9 reference), which enumerates *all* alignments within the mismatch/gap
budget (not just the optimal-scoring one). Only the CLI and output formatting
(tab-delimited, underscore column names) were adapted here.
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
    RNA, DNA, max_mismatches=4, max_DNA_gaps=1, max_RNA_gaps=1,
):
    # Budget-bounded enumeration: the DP tracks every (mismatch, dna_gap, rna_gap)
    # count reachable at each cell and backtracks all of them, so non-optimal but
    # in-budget alignments are enumerated too (not just the single best-scoring
    # alignment). This matches the reference generator.
    m, n = len(RNA), len(DNA)
    dp_score = [[float("inf")] * (n + 1) for _ in range(m + 1)]
    dp_cells = [[{} for _ in range(n + 1)] for _ in range(m + 1)]
    dp_score[0][0] = 0
    dp_cells[0][0] = {(0, 0, 0): []}
    for i in range(1, m + 1):
        if i <= max_DNA_gaps:
            dp_score[i][0] = i
            dp_cells[i][0] = {(0, i, 0): [("up", i - 1, 0)]}
    for j in range(1, n + 1):
        if j <= max_RNA_gaps:
            dp_score[0][j] = j
            dp_cells[0][j] = {(0, 0, j): [("left", 0, j - 1)]}
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            candidates = []
            for (mm, dg, rg) in dp_cells[i - 1][j - 1]:
                mis = 0 if iupac_match(RNA[i - 1], DNA[j - 1]) else 1
                new = (mm + mis, dg, rg)
                if new[0] <= max_mismatches:
                    candidates.append((dp_score[i - 1][j - 1] + mis, new, "diag", i - 1, j - 1))
            for (mm, dg, rg) in dp_cells[i - 1][j]:
                new = (mm, dg + 1, rg)
                if new[1] <= max_DNA_gaps:
                    candidates.append((dp_score[i - 1][j] + 1, new, "up", i - 1, j))
            for (mm, dg, rg) in dp_cells[i][j - 1]:
                if i == m:
                    new, cost = (mm, dg, rg), dp_score[i][j - 1]
                else:
                    new, cost = (mm, dg, rg + 1), dp_score[i][j - 1] + 1
                if new[2] <= max_RNA_gaps:
                    candidates.append((cost, new, "left", i, j - 1))
            if not candidates:
                continue
            dp_score[i][j] = min(c[0] for c in candidates)
            cell = {}
            for cost, cnt, move, pi, pj in candidates:
                cell.setdefault(cnt, []).append((move, pi, pj))
            dp_cells[i][j] = cell

    alignments = []
    seen = set()

    def backtrack(i, j, cnt, a1, a2):
        if i == 0 and j == 0:
            aln1, aln2 = a1[::-1], a2[::-1]
            mismatches = sum(0 if iupac_match(x, y) else 1
                             for x, y in zip(aln1, aln2) if x != "-" and y != "-")
            dna_gaps = aln2.count("-")
            rna_gaps = aln1.rstrip("-").count("-")
            key = (aln1, aln2, mismatches, dna_gaps, rna_gaps)
            if key not in seen:
                seen.add(key)
                alignments.append((aln1, aln2, mismatches, rna_gaps, dna_gaps))
            return
        for move, pi, pj in dp_cells[i][j][cnt]:
            if move == "diag":
                mis = 0 if iupac_match(RNA[i - 1], DNA[j - 1]) else 1
                backtrack(pi, pj, (cnt[0] - mis, cnt[1], cnt[2]), a1 + RNA[i - 1], a2 + DNA[j - 1])
            elif move == "up":
                backtrack(pi, pj, (cnt[0], cnt[1] - 1, cnt[2]), a1 + RNA[i - 1], a2 + "-")
            else:  # left
                prev = cnt if i == m else (cnt[0], cnt[1], cnt[2] - 1)
                backtrack(pi, pj, prev, a1 + "-", a2 + DNA[j - 1])

    for final_cnt in dp_cells[m][n]:
        backtrack(m, n, final_cnt, "", "")
    return alignments or None


def reverse_complement(seq):
    complement = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
    }
    return "".join(complement.get(b, b) for b in reversed(seq.upper()))


def keep_alignment(rna_out, dna_out, pam, pam_5prime):
    """Off-target validity filter matching CRISPRme semantics:
      - the concrete PAM bases must match the DNA (no mismatch/gap in the PAM), and
      - the alignment must not start or end with an RNA gap.
    `pam` is the concrete PAM with wildcard N's removed (e.g. "GG" for NGG, or the
    full "TTTV" for Cas12a). Without these filters the raw exhaustive aligner
    over-produces PAM-mismatched and edge-gap alignments the reference excludes.
    """
    if not rna_out or rna_out[0] == "-" or rna_out[-1] == "-":
        return False
    if pam:
        d = dna_out.upper()
        seg = d[:len(pam)] if pam_5prime else d[-len(pam):]
        if len(seg) < len(pam) or "-" in seg:
            return False
        if not all(iupac_match(g, x) for g, x in zip(pam, seg)):
            return False
    return True


def scan(RNA, DNA, Strand, max_dna_gaps, max_rna_gaps, max_mismatches, chrom,
         writer, pam="", pam_5prime=False):
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
            r_out, d_out = r_aln[:actual_len], d_aln[:actual_len]
            if not keep_alignment(r_out, d_out, pam, pam_5prime):
                continue
            if Strand == "+":
                start = starting_pos
                end = starting_pos + actual_len
            else:
                start = len(DNA) - starting_pos - actual_len + dg
                end = len(DNA) - starting_pos + dg
            writer.writerow([chrom, r_out, d_out, Strand, start, end, mm, rg, dg])


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
    p.add_argument("--pam", default="",
                   help="Concrete PAM bases with wildcard N removed (e.g. GG for "
                        "NGG, TTTV for Cas12a); empty disables PAM filtering")
    p.add_argument("--pam-5prime", action="store_true",
                   help="PAM is at the 5' end of the guide (Cas12a); default 3' (Cas9)")
    args = p.parse_args()

    DNA = str(SeqIO.read(args.fasta, "fasta").seq)
    with open(args.output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            ["CHR", "RNA", "DNA", "Strand", "Start", "END",
             "mismatches", "gaps_in_RNA", "gaps_in_DNA"]
        )
        for strand_seq, strand in ((DNA, "+"), (reverse_complement(DNA), "-")):
            scan(args.rna, strand_seq, strand, args.max_dna_gaps, args.max_rna_gaps,
                 args.max_mismatches, args.chrom, writer,
                 pam=args.pam.upper(), pam_5prime=args.pam_5prime)


if __name__ == "__main__":
    main()
