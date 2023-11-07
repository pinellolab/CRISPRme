"""
This file contains the main functionality for designing CRISPR/Cas9 guides based on given input files.

Functions:
- `parse_commandline() -> argparse.Namespace`: Parses the command line arguments for the design_guides script.
- `check_input_args(args: argparse.Namespace) -> None`: Checks the validity of the input arguments.
- `design_guides(pamfile: str, genome: str, bedfile: str, mm_max: int, outname: str) -> None`: Designs guides for CRISPR/Cas9 genome editing.
- `main() -> None`: The main entry point for the design_guides script.

Example:
    ```python
    main()
    ```
"""


from utils import reverse_complement, read_pam, read_genome, read_coordinates
from encoder import encode_pam, encode_genome
from search import match_genome
from pam import PAM
from reporter import recover_guides

from typing import Tuple, List

import argparse
import pysam
import time
import sys
import os


def parse_commandline() -> argparse.Namespace:
    """
    Parses the command line arguments for the design_guides script.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.

    Example:
        ```python
        args = parse_commandline()
        print(args.fasta)  # Output: 'input.fasta'
        print(args.pam)  # Output: 'pam.txt'
        print(args.coordinates)  # Output: 'coordinates.bed'
        print(args.mm)  # Output: 1
        print(args.output_prefix)  # Output: 'guides_out'
        ```
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Given input fasta file, chromosome coordinates file (BED "
        "format [optional]) and PAM sequence, find all possible sgRNA guides in "
        "the fasta file (SUPPORTS IUPAC NOTATIONS)",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        metavar="FASTA-FILE",
        help="FASTA file contaning either a single chromosome or a target "
        "sequence extracted from genome",
    )
    parser.add_argument(
        "--pam",
        type=str,
        required=True,
        metavar="PAM-FILE",
        help="PAM file (CRISPRitz/CRISPRme format)",
    )
    parser.add_argument(
        "--coordinates",
        type=str,
        required=False,
        default="",
        metavar="BED-FILE",
        help="Coordinate file in BED format (Please Note: only one chromosome "
        "allowed at a time)",
    )
    parser.add_argument(
        "--mm",
        type=int,
        required=False,
        default=0,
        metavar="MM-NUM",
        help="Maximum allowed mismatches in the PAM sequence",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        required=False,
        default="guides_out",
        dest="output_prefix",
        help="Output file prefix: the generated files include the guides and a "
        "detailed report",
    )
    return parser.parse_args()


def check_input_args(args: argparse.Namespace) -> None:
    """
    Checks the validity of the input arguments.

    Args:
        args (argparse.Namespace): An object containing the parsed command line arguments.

    Returns:
        None

    Raises:
        FileNotFoundError: If any of the input files specified in the arguments are not found.
        ValueError: If any of the input files specified in the arguments are empty or if the number of mismatches is negative.

    Example:
        ```python
        args = parse_commandline()
        check_input_args(args)
        ```
    """

    if not os.path.isfile(args.fasta):
        raise FileNotFoundError(f"{args.fasta} not found!")
    if os.stat(args.fasta).st_size == 0:
        raise ValueError(f"{args.fasta} is empty!")
    if not os.path.isfile(args.pam):
        raise FileNotFoundError(f"{args.pam} not found!")
    if os.stat(args.pam).st_size == 0:
        raise ValueError(f"{args.pam} is empty!")
    if not os.path.isfile(args.coordinates):
        raise FileNotFoundError(f"{args.coordinates} not found!")
    if os.stat(args.coordinates).st_size == 0:
        raise ValueError(f"{args.coordinates} is empty")
    if args.mm < 0:
        raise ValueError(f"Forbidden number of mismatches ({args.mm})!")


def design_guides(
    pamfile: str, genome: str, bedfile: str, mm_max: int, outname: str
) -> None:
    """
    Designs guides for CRISPR/Cas9 genome editing based on the given PAM sequence, genome sequence, BED file, maximum number of mismatches, and output file name.

    Args:
        pamfile (str): The path to the PAM file.
        genome (str): The path to the genome file.
        bedfile (str): The path to the BED file (optional).
        mm_max (int): The maximum number of mismatches allowed.
        outname (str): The name of the output file.

    Returns:
        None

    Example:
        ```python
        pamfile = "path/to/pam.txt"
        genome = "path/to/genome.fasta"
        bedfile = "path/to/bedfile.bed"
        mm_max = 2
        outname = "output.txt"
        design_guides(pamfile, genome, bedfile, mm_max, outname)
        ```
    """

    pam = read_pam(pamfile)  # read PAM
    pam_bits = encode_pam(pam.pam), encode_pam(reverse_complement(pam.pam))
    sequence = read_genome(genome)
    if bedfile:  # extract subsequences
        coordinates = read_coordinates(bedfile)
        sequence_bits = {
            f"{coords[0]}:{coords[1]}-{coords[2]}": encode_genome(
                sequence.fetch(coords[0], coords[1], coords[2]).upper()
            )
            for coords in coordinates
        }
    else:
        sequence_bits = {
            seqname: encode_genome(sequence[seqname].upper())
            for seqname in sequence.references
        }
    # search PAM occurrences within the query sequences
    pam_positions = {
        seqname: match_genome(
            sequence_bits[seqname],
            pam_bits,
            pam.length,
            pam.full_length,
            mm_max,
            pam.pam_at_beginning,
        )
        for seqname in sequence_bits
    }
    # write report and guides file
    recover_guides(pam_positions, sequence, pam, outname)


def main():
    """
    The main entry point for the design_guides script.

    Returns:
        None

    Example:
        ```python
        main()
        ```
    """

    start = time.time()
    args = parse_commandline()  # parse command line args
    check_input_args(args)  # check input args consistency
    # run guide design
    sys.stderr.write(f"Searching guides in {args.fasta}...\n")
    design_guides(args.pam, args.fasta, args.coordinates, args.mm, args.output_prefix)
    sys.stderr.write(f"Elapsed time {(time.time() - start):.2f}s\n")


# entry point
if __name__ == "__main__":
    main()
