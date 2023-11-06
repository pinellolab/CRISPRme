"""
"""

from utils import IUPACTABLE, reverse_complement
from pam import PAM

from typing import Tuple, Dict, List

import pysam

REPORTCOLS = ["sgRNA", "PAM", "sgRNA_length", "PAM_length", "chr", "pos", "strand"]


def is_valid_sequence(sequence: str) -> bool:
    return "N" not in sequence


def process_seqname(coords: str) -> Tuple[str, int]:
    fields = coords.strip().split(":")
    if len(fields) > 1:  # coordinates given as input
        return fields[0], int(fields[1].split("-")[0])
    return fields[0], 0


def expand_iupac(nt: str, i: int) -> str:
    if nt in IUPACTABLE:
        return IUPACTABLE[nt]
    raise ValueError(f"Forbidden nucleotide ({nt}) at position {i}")


def explode_iupac_sequence(sequence: str) -> List[str]:
    xpanded_sequences = [""]  # expanded IUPAC sequences
    for i, nt in enumerate(sequence):
        xpanded_nt = expand_iupac(nt, i)  # expand input nt in its IUPAC variants
        xpanded_sequences = [s + b for s in xpanded_sequences for b in xpanded_nt]
    assert len(xpanded_sequences) >= 1
    return xpanded_sequences


def pad_guide(guide: str, pam_length: int, pam_at_beginning: bool):
    if pam_at_beginning:
        return "N" * pam_length + guide
    return guide + "N" * pam_length


def extract_guide_pam(
    seqname: str,
    positions: Tuple[List[int], List[int]],
    genome: pysam.FastaFile,
    pam: PAM,
) -> List[List[str]]:
    report = []
    seqname, pos = process_seqname(
        seqname
    )  # recover seqname and relative start position
    for p in positions[0]:  # forward
        p += pos
        start = p + len(pam) if pam.pam_at_beginning else p
        end = (
            p + len(pam) + pam.guide_length
            if pam.pam_at_beginning
            else p + pam.guide_length
        )
        guide_iupac = genome.fetch(seqname, start, end).upper()  # allow IUPAC chars
        if is_valid_sequence(guide_iupac):
            guides = explode_iupac_sequence(guide_iupac.upper())  # replace IUPAC chars
            start = p if pam.pam_at_beginning else p + pam.guide_length
            end = (
                p + pam.guide_length
                if pam.pam_at_beginning
                else p + len(pam) + pam.guide_length
            )
            pam_seq = genome.fetch(seqname, start, end)
            if is_valid_sequence(pam_seq):
                for guide in guides:
                    report.append(
                        list(
                            map(
                                str,
                                [
                                    pad_guide(guide, len(pam), pam.pam_at_beginning),
                                    pam_seq,
                                    pam.guide_length,
                                    len(pam),
                                    seqname,
                                    p,
                                    "+",
                                ],
                            )
                        )
                    )
    for p in positions[1]:  # reverse
        p += pos
        start = p if pam.pam_at_beginning else p + len(pam)
        end = (
            p + pam.guide_length
            if pam.pam_at_beginning
            else p + len(pam) + pam.guide_length
        )
        print(p, start, end)
        guide_iupac = reverse_complement(genome.fetch(seqname, start, end).upper())
        if is_valid_sequence(guide_iupac):
            guides = explode_iupac_sequence(guide_iupac)  # replace IUPAC chars
            start = p + pam.guide_length if pam.pam_at_beginning else p
            end = (
                p + len(pam) + pam.guide_length
                if pam.pam_at_beginning
                else p + len(pam)
            )
            pam_seq = genome.fetch(seqname, start, end)
            if is_valid_sequence(pam_seq):
                for guide in guides:
                    report.append(
                        list(
                            map(
                                str,
                                [
                                    pad_guide(guide, len(pam), pam.pam_at_beginning),
                                    pam_seq,
                                    pam.guide_length,
                                    len(pam),
                                    seqname,
                                    p,
                                    "-",
                                ],
                            )
                        )
                    )
    return report


def write_report(report: List[List[str]], outname: str) -> None:
    try:
        guides_file = f"{outname}.guides.txt"
        report_file = f"{outname}.summary.tsv"
        with open(guides_file, mode="w") as outguide, open(
            report_file, mode="w"
        ) as outreport:
            # write report header
            header = "\t".join(REPORTCOLS)
            outreport.write(f"{header}\n")
            for line in report:
                rline = "\t".join(line)  # report line
                outguide.write(f"{line[0]}\n")
                outreport.write(f"{rline}\n")
    except IOError as e:
        raise IOError("Report writing failed!") from e


def recover_guides(
    pam_positions: Dict[str, Tuple[List[int], List[int]]],
    genome: pysam.FastaFile,
    pam: PAM,
    outname: str,
) -> None:
    report = []
    for seqname in pam_positions:
        report += extract_guide_pam(seqname, pam_positions[seqname], genome, pam)
    write_report(report, outname)


def write_results(
    guide_list: list, guide_file: str, positions: Tuple, chr_name: str
) -> None:
    """_summary_
    Args:
        guide_list (list): _description_
        guide_file (str): _description_
        positions (Tuple): _description_
        chr_name (str): _description_
    """
    guide_file_out = open(guide_file, "w")
    guide_summary_out = open(guide_file.replace(".txt", "") + ".summary.tsv", "w")
    guide_summary_out.write("sgRNA\tPAM\tsgRNA_length\tPAM_length\tchr\tstrand\n")

    for i, guide in enumerate(guide_list):
        guide_file_out.write(guide[0] + "\n")
        guide_summary_out.write(
            guide[0]
            + "\t"
            + guide[1]
            + "\t"
            + str(len(guide[0]) - len(guide[1]))
            + "\t"
            + str(len(guide[1]))
            + "\t"
            + chr_name
            + "\t"
            + guide[2]
            + "\n"
        )
    guide_file_out.close()
    guide_summary_out.close()
