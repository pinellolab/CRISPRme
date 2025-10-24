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
    """Clean an interval value string.

    This function removes commas, spaces, and periods from an interval value string
    and converts it to an integer.

    Args:
        value: The interval value string.

    Returns:
        The cleaned interval value as an integer.
    """
    return int(value.replace(",", "").replace(" ", "").replace(".", ""))


def retrieve_coordinates(input_range: str) -> Tuple[str, int, int]:
    """Retrieve coordinates from an input range string.

    This function parses an input range string in the format "chrom:start-stop"
    and returns the chromosome, start, and stop coordinates as a tuple.

    Args:
        input_range: The input range string.

    Returns:
        A tuple containing the chromosome, start coordinate, and stop coordinate.
    """
    chrom, coordinates = input_range.split(":")  # extract chromosome
    start = clean_interval_value(coordinates.split("-")[0])  # extract start position
    stop = clean_interval_value(coordinates.split("-")[1])  # extract stop position
    return chrom, start, stop


def write_mock_bed(seqname: str, chrom: str, start: int, stop: int) -> str:
    """Write a mock BED file.

    This function creates a temporary BED file containing a single interval
    defined by the provided sequence name, chromosome, start, and stop coordinates.

    Args:
        seqname: The name of the sequence.
        chrom: The chromosome name.
        start: The start coordinate.
        stop: The stop coordinate.

    Returns:
        The path to the created BED file.

    Raises:
        IOError: If an error occurs while writing the BED file.
    """
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
    """Retrieve chromosome FASTA files.

    This function retrieves a list of chromosome FASTA files from the specified
    genome directory, excluding FASTA index files (.fai).

    Args:
        genome: The name of the genome.

    Returns:
        A list of chromosome FASTA filenames.

    Raises:
        FileNotFoundError: If the genome directory is not found.
    """
    genomedir = os.path.join(os.getcwd(), GENOMES_DIR, genome)
    if not os.path.isdir(genomedir):
        raise FileNotFoundError(f"Cannot find Genome folder {genomedir}")
    return [
        f
        for f in os.listdir(genomedir)
        if os.path.isfile(os.path.join(genomedir, f)) and not f.endswith(".fai")
    ]


def getfasta(bedfname: str, chromfasta: str) -> str:
    """Extract sequence from FASTA using bedtools getfasta.

    This function extracts a sequence from a FASTA file based on coordinates
    specified in a BED file using the bedtools getfasta utility.

    Args:
        bedfname: Path to the BED file.
        chromfasta: Path to the chromosome FASTA file.

    Returns:
        The extracted sequence string.
    """
    # extract sequence from input coordinates
    sequence = subprocess.check_output(
        [f"{GETFASTA} -fi {chromfasta} -bed {bedfname}"], shell=True
    ).decode("utf-8")
    return sequence.split("\n")[1].strip()  # remove header


def extract_sequence(seqname: str, input_range: str, genome: str) -> str:
    """Extract a genomic sequence.

    This function extracts a genomic sequence from a specified FASTA file based on
    the provided sequence name, input range, and genome.

    Args:
        seqname: The name of the sequence.
        input_range: The genomic range in "chrom:start-stop" format.
        genome: The name of the genome.

    Returns:
        The extracted sequence string.

    Raises:
        FileNotFoundError: If the FASTA file for the specified chromosome is not found.
    """
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
    """Generate all possible PAM sequences from an IUPAC PAM string.

    This function takes an IUPAC PAM string as input and generates a list of all
    possible PAM sequences by expanding the IUPAC characters using the PAMDICT.

    Args:
        pam: The IUPAC PAM string.

    Returns:
        A list of all possible PAM sequences.
    """
    pam_chars = [PAMDICT[nt] for nt in pam]  # expand pam characters
    return ["".join(e) for e in product(*pam_chars)]


def search_pam_fwd(
    iupac_pam: List[str], sequence: str, pamstart: bool, guidelen: int, pamlen: int
) -> List[str]:
    """Search for PAMs in the forward strand.

    This function searches for all occurrences of PAM sequences from a given list
    of IUPAC PAMs in the forward strand of a given sequence.

    Args:
        iupac_pam: A list of IUPAC PAM strings.
        sequence: The input sequence string.
        pamstart: True if PAM is upstream of the guide, False otherwise.
        guidelen: The length of the guide sequence.
        pamlen: The length of the PAM sequence.

    Returns:
        A list of guide sequences found in the forward strand.
    """
    guides = []  # list of retrieved guides
    for pam in iupac_pam:  # iterate over all possible pams
        for i in [m.start() for m in re.finditer("(?=" + pam + ")", sequence)]:
            if pamstart:  # pam upstreeam the guide
                if i <= (len(sequence) - guidelen - pamlen):
                    guides.append(sequence[i + pamlen : i + pamlen + guidelen])
            elif i >= guidelen:  # pam downstream
                guides.append(sequence[i - guidelen : i])
    return guides


def search_pam_rev(
    iupac_pam: List[str], sequence: str, pamstart: bool, guidelen: int, pamlen: int
) -> List[str]:
    """Search for PAMs in the reverse strand.

    This function searches for all occurrences of PAM sequences from a given list
    of IUPAC PAMs in the reverse strand of a given sequence.

    Args:
        iupac_pam: A list of IUPAC PAM strings.
        sequence: The input sequence string.
        pamstart: True if PAM is upstream of the guide, False otherwise.
        guidelen: The length of the guide sequence.
        pamlen: The length of the PAM sequence.

    Returns:
        A list of guide sequences found in the reverse strand.
    """
    guides = []  # list of retrieved guides
    for pam in iupac_pam:  # iterate over all possible pams
        for i in [m.start() for m in re.finditer("(?=" + pam + ")", sequence)]:
            if pamstart:  # pam upstreeam the guide
                if i >= guidelen:
                    guide = str(Seq(sequence[i - guidelen : i]).reverse_complement())
                    guides.append(guide)
            elif i <= (len(sequence) - guidelen - pamlen):  # pam downstream
                guide = str(
                    Seq(
                        sequence[i + pamlen : i + pamlen + guidelen]
                    ).reverse_complement()
                )
                guides.append()
    return guides


def get_guides(sequence: str, pam: str, guidelen: str, pamstart: bool) -> List[str]:
    """Retrieve CRISPR guides from a given sequence.

    This function searches for CRISPR guides in both the forward and reverse strands
    of a given sequence, considering the provided PAM sequence, guide length, and
    PAM position relative to the guide.

    Args:
        sequence: The input sequence string.
        pam: The PAM sequence.
        guidelen: The length of the guide sequence.
        pamstart: True if PAM is upstream of the guide, False otherwise.

    Returns:
        A list of guide sequences found in both strands.
    """
    pamlen = len(pam)  # pam length
    guidelen = int(guidelen)  # cast to int
    iupac_pam = generate_iupac_pam(pam)  # compute pam exploding iupac chars
    iupac_pam_rev = generate_iupac_pam(str(Seq(pam).reverse_complement()))
    sequence = sequence.upper()  # force uppercase
    # search guides on forward and reverse strands
    guides_fwd = search_pam_fwd(iupac_pam, sequence, pamstart, guidelen, pamlen)
    guides_rev = search_pam_rev(iupac_pam_rev, sequence, pamstart, guidelen, pamlen)
    return guides_fwd + guides_rev
