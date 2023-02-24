"""
"""

from crispritz import crispritz_index_genome
from utils import write

from typing import List, Tuple
from multiprocessing import Pool

import os


def index(
    indels_folder: str,
    gname: str,
    vname: str,
    pam: str,
    pamfile: str,
    bmax: int,
    threads: int,
    chrom: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Use the crispritz index-genome functionality to generate an index for the
    genome sequence that has been enriched with indels

    :param indels_folder: folder containing the genome enriched with indels
    :type indels_folder: str
    :param gname: genome name
    :type gname: str
    :param vname: variants dataset name
    :type vname: str
    :param pam: PAM sequence
    :type pam: str
    :param pamfile: PAM file
    :type pamfile: str
    :param bmax: max bulges
    :type bmax: int
    :param threads: threads
    :type threads: int
    :param chrom: chromosome
    :type chrom: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    # recover the folder storing the genome enriched with indels and the corresponding file
    genome = os.path.join(f"{gname}+{vname}_INDELS", f"{pam}_{bmax}_fake{chrom}")
    indels_file = os.path.join(indels_folder, f"fake_{vname}_{chrom}")
    if verbosity > 1:
        write(f"Indexing chr{chrom} indels")
    # run crispritz index-genome on the input data
    crispritz_index_genome(
        genome, indels_file, pamfile, bmax, threads, verbosity, debug
    )


def index_indels(
    genome: str,
    pamfile: str,
    pam: str,
    gname: str,
    vname: str,
    bmax: int,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Construct a ternary search tree (TST) to create an index of the genome
    sequence that has been enriched with indels present in the input variants
    dataset. The function runs in parallel creating a user-defined number of
    threads

    :param genome: enriched genome folder
    :type genome: str
    :param pamfile: PAM file
    :type pamfile: str
    :param pam: PAM sequence
    :type pam: str
    :param gname: genome name
    :type gname: str
    :param vname: variants dataset
    :type vname: str
    :param bmax: max bulges
    :type bmax: int
    :param threads: threads
    :type threads: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    # recover the chromosomes in the indels folder
    chroms = [f.split("_")[-1] for f in os.listdir(genome) if "chr" in f]
    # call crispritz index-genome in parallel
    args = [
        (genome, gname, vname, pam, pamfile, bmax, 1, chrom, verbosity, debug)
        for chrom in chroms
    ]  # force one thread
    with Pool(processes=threads) as pool:  # use threads processes
        pool.starmap(index, args)
