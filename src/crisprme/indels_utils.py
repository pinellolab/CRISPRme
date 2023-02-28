"""
"""

from crispritz import crispritz_index_genome, crispritz_search
from verbosity_handler import time_verbosity, write_verbosity
from utils import exception_handler, remove

from typing import List, Tuple
from multiprocessing import Pool
from glob import glob

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
    write_verbosity(f"indexing chr{chrom} indels", verbosity, 1, debug)
    # run crispritz index-genome on the input data
    crispritz_index_genome(
        genome, indels_file, pamfile, bmax, threads, 0, debug
    )  # force quiet run


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


def search(tst_variants_indels: str, pam: str, pamfile: str, bmax: int, chrom: str, guide: str, guidefile: str, mm: int, bDNA: int, bRNA: int, threads: int, verbosity: int, debug: bool):
    """Use crispritz search functionality to search potential off-targets on 
    variant genome indels in parallel

    :param tst_variants_indels: TST index
    :type tst_variants_indels: str
    :param pam: PAM
    :type pam: str
    :param pamfile: PAM file
    :type pamfile: str
    :param bmax: maximum bulge value
    :type bmax: int
    :param chrom: chromosome
    :type chrom: str
    :param guide: guide 
    :type guide: str
    :param guidefile: guide file
    :type guidefile: str
    :param mm: mismatches
    :type mm: int
    :param bDNA: DNA bulges
    :type bDNA: int
    :param bRNA: RNA bulges
    :type bRNA: int
    :param threads: threads
    :type threads: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    # recover chromosome associated TST 
    tst_variants_indels = os.path.join(tst_variants_indels, f"{pam}_{bmax}_fake{chrom}")
    reportname = f"fake{chrom}_{pam}_{guide}_{mm}_{bDNA}_{bRNA}"  # chrom report
    write_verbosity(f"Searching off-targets on {chrom} indels", verbosity, 1, debug)
    # call crispritz search
    crispritz_search(tst_variants_indels, pamfile, guidefile, reportname, mm, bDNA, bRNA, threads, 0, debug)  # force quiet run

def process_report_indels(gname: str, vname: str, pam: str, guide: str, mm: int, bDNA: int, bRNA: int, verbosity: int, debug: bool) -> None:
    """Merge single chromosome off-targets search reports

    :param gname: genome name
    :type gname: str
    :param vname: variants dataset
    :type vname: str
    :param pam: PAM
    :type pam: str
    :param guide: guide
    :type guide: str
    :param mm: mismatches
    :type mm: int
    :param bDNA: DNA bulges
    :type bDNA: int
    :param bRNA: RNA bulges
    :type bRNA: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    :raises OSError: raise on errors occurring while merging reports
    """
    write_verbosity("Constructing off-targets report on indels", verbosity, 2, debug)
    start_report = time_verbosity(verbosity, 2, debug)
    # define final report name
    reportfinal = f"indels_{gname}+{vname}_{pam}_{guide}_{mm}_{bDNA}_{bRNA}.targets.txt"
    try:  # start merging individual reports
        with open(reportfinal, mode="a+") as outfile:  # open final report
            header = ""  # header shared by all partial reports
            for report in glob(f"fakechr*_{pam}_{guide}_{mm}_{bDNA}_{bRNA}.targets.txt"):
                with open(report, mode="r") as infile:
                    lines = infile.readlines()
                    header = header or lines[0]  # store header if not empty
                    for line in lines[1:]:
                        outfile.write(line)
                # remove(report, debug)  # remove partial report
            outfile.seek(0, 0)  # move pointer to start
            content = outfile.read()
            outfile.seek(0, 0)
            outfile.write(header.rstrip("\r\n") + "\n" + content)  # add header to report
    except OSError:
        exception_handler(OSError, "An error occurred while constructing the summary report on indels", debug)
    stop_report = time_verbosity(verbosity, 2, debug)
    write_verbosity(f"Off-targets report on indels construction took {(stop_report - start_report):.2f}s", verbosity, 2, debug)


def search_indels(tst_variants_indels: str, gname: str, vname: str, pam: str, pamfile: str, mm: int, bDNA: int, bRNA: int, guide: str, guidefile: str, threads: int, verbosity: int, debug: bool) -> None:
    """Search variant genome indels for potential off-targets

    :param tst_variants_indels: TST index
    :type tst_variants_indels: str
    :param gname: genome name
    :type gname: str
    :param vname: variants dataset
    :type vname: str
    :param pam: PAM
    :type pam: str
    :param pamfile: PAM file
    :type pamfile: str
    :param mm: mismatches
    :type mm: int
    :param bDNA: DNA bulges
    :type bDNA: int
    :param bRNA: RNA bulges
    :type bRNA: int
    :param guide: guide
    :type guide: str
    :param guidefile: guide file
    :type guidefile: str
    :param threads: threads
    :type threads: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    """
    # recover chromosomes with indels and their TST index
    chroms = [d.split("_")[-1].replace("fake", "") for d in os.listdir(tst_variants_indels) if "chr" in d]
    bmax = max(bDNA, bRNA)
    args = [
        (tst_variants_indels, pam, pamfile, bmax, chrom, guide, guidefile, mm, bDNA, bRNA, 1, verbosity, debug)
        for chrom in chroms
    ]  # force one thread 
    with Pool(processes=threads) as pool:  # use threads processes
        pool.starmap(search, args)
    # construct the final off-targets report on indels
    process_report_indels(gname, vname, pam, guide, mm, bDNA, bRNA, verbosity, debug) 
