""" """

# from pages_utils import GENOMES_DIR

from typing import Tuple, List
from itertools import product
from Bio.Seq import Seq

import subprocess
import re
import os

GETFASTA = "bedtools getfasta"
PAMDICT = {
    "A": "ARWMDHV",
    "C": "CYSMBHV",
    "G": "GRSKBDV",
    "T": "TYWKBDH",
    "R": "ARWMDHVSKBG",
    "Y": "CYSMBHVWKDT",
    "S": "CYSMBHVKDRG",
    "W": "ARWMDHVYKBT",
    "K": "GRSKBDVYWHT",
    "M": "ARWMDHVYSBC",
    "B": "CYSMBHVRKDGWT",
    "D": "ARWMDHVSKBGYT",
    "H": "ARWMDHVYSBCKT",
    "V": "ARWMDHVYSBCKG",
    "N": "ACGTRYSWKMBDHV",
}

GENOMES_DIR = "Genomes"


def clean_interval_value(value: str) -> int:
    return int(value.replace(",", "").replace(" ", "").replace(".", ""))


def retrieve_coordinates(input_range: str) -> Tuple[str, int, int]:
    chrom, coordinates = input_range.split(":")  # extract chromosome
    start = clean_interval_value(coordinates.split("-")[0])  # extract start position
    stop = clean_interval_value(coordinates.split("-")[1])  # extract stop position
    return chrom, start, stop


def write_mock_bed(seqname: str, chrom: str, start: int, stop: int) -> str:
    # define mock bed file name
    bedfname = os.path.join(os.getcwd(), f"{seqname}.bed")
    try:
        with open(bedfname, mode="w") as outfile:
            outfile.write(f"{chrom}\t{start}\t{stop}\n")
    except OSError as e:
        raise IOError(
            f"An error occurred while writing mock bed file {bedfname}"
        ) from e
    assert os.path.isfile(bedfname) and os.stat(bedfname).st_size > 0
    return bedfname


def retrieve_chromosomes_fasta(genome: str) -> List[str]:
    genomedir = os.path.join(os.getcwd(), GENOMES_DIR, genome)
    if not os.path.isdir(genomedir):
        raise FileNotFoundError(f"Cannot find Genome folder {genomedir}")
    return [
        f
        for f in os.listdir(genomedir)
        if os.path.isfile(os.path.join(genomedir, f)) and not f.endswith(".fai")
    ]


def getfasta(bedfname: str, chromfasta: str) -> str:
    # extract sequence from input coordinates
    sequence = subprocess.check_output(
        [f"{GETFASTA} -fi {chromfasta} -bed {bedfname}"], shell=True
    ).decode("utf-8")
    return sequence.split("\n")[1].strip()  # remove header


def extract_sequence(seqname: str, input_range: str, genome: str) -> str:
    # replace spaces with underscores in sequence name
    seqname = "_".join(seqname.replace(">", "").split())
    # retrieve genomic coordinates
    chrom, start, stop = retrieve_coordinates(input_range)
    # write mock bed to extract sequances
    bedfname = write_mock_bed(seqname, chrom, start, stop)
    chromosomes_fasta = retrieve_chromosomes_fasta(genome)  # retrieve chromosomes fasta
    # extract sequence
    chromfasta = None
    for c in chromosomes_fasta:
        if c.startswith(chrom):  # found match
            chromfasta = os.path.join(os.getcwd(), GENOMES_DIR, genome, c)
    if chromfasta is None or not os.path.isfile(chromfasta):
        raise FileNotFoundError(f"Fasta file not found for chromosome {chrom}")
    sequence = getfasta(bedfname, chromfasta)  # recover sequence
    os.remove(f"{chromfasta}.fai")  # remove fasta index
    os.remove(bedfname)  # remove mock bed
    return sequence


def generate_iupac_pam(pam: str):
    pam_chars = [PAMDICT[nt] for nt in pam]  # expand pam characters
    return ["".join(e) for e in product(*pam_chars)]


def search_pam_fwd(
    iupac_pam: List[str], sequence: str, pamstart: bool, guidelen: int, pamlen: int
) -> List[str]:
    guides = []  # list of retrieved guides
    for pam in iupac_pam:  # iterate over all possible pams
        for i in ([m.start() for m in re.finditer("(?=" + pam + ")", sequence)]):
            if pamstart:  # pam upstreeam the guide
                if i <= (len(sequence) - guidelen - pamlen):
                    guides.append(sequence[i + pamlen : i + pamlen + guidelen])
            elif i >= guidelen:  # pam downstream
                guides.append(sequence[i - guidelen : i])
    return guides

def search_pam_rev(
    iupac_pam: List[str], sequence: str, pamstart: bool, guidelen: int, pamlen: int
) -> List[str]:
    guides = []  # list of retrieved guides
    for pam in iupac_pam:  # iterate over all possible pams
        for i in ([m.start() for m in re.finditer("(?=" + pam + ")", sequence)]):
            if pamstart:  # pam upstreeam the guide
                if i >= guidelen:
                    guide = str(Seq(sequence[i - guidelen : i]).reverse_complement())
                    guides.append(guide)
            elif i <= (len(sequence) - guidelen - pamlen):  # pam downstream
                guide = str(Seq(sequence[i + pamlen : i + pamlen + guidelen]).reverse_complement())
                guides.append()
    return guides


def get_guides(sequence: str, pam: str, guidelen: str, pamstart: bool):
    pamlen = len(pam)  # pam length
    guidelen = int(guidelen)  # cast to int
    iupac_pam = generate_iupac_pam(pam)  # compute pam exploding iupac chars
    iupac_pam_rev = generate_iupac_pam(str(Seq(pam).reverse_complement()))
    sequence = sequence.upper()  # force uppercase
    # search guides on forward and reverse strands
    guides_fwd = search_pam_fwd(iupac_pam, sequence, pamstart, guidelen, pamlen)
    guides_rev = search_pam_rev(iupac_pam_rev, sequence, pamstart, guidelen, pamlen)
    return guides_fwd + guides_rev


if __name__ == "__main__":
    sequence = extract_sequence(">sequence1", "chrx:51-74", "test_genome")
    print(sequence)
    print(get_guides(sequence, "NGG", 10, False))
