"""
"""

from crisprme_errors import TSTError
from crispritz import (
    ENRICHED_GENOME,
    GENOME_SNPS,
    crispritz_add_variants,
    crispritz_index_genome,
)
from indels_index import index_indels
from verbosity_handler import write_verbosity, time_verbosity
from utils import (
    CRISPRME_DIRS,
    LOG,
    exception_handler,
    move,
    raise_warning,
    remove_dir,
    write,
)

from typing import List, Optional
from glob import glob
from time import ctime, time

import sys
import re
import os


GENOMELIB = "genome_library"  # store indels indexes


def recover_fasta(genome: str, verbosity: bool, debug: bool) -> List[str]:
    """Recover FASTA files in the genome directory (assumes FA or FASTA
    extensions)

    :param genome: genome directory
    :type genome: str
    :param verbosity: verbosity level
    :type verbosity: bool
    :param debug: debug mode
    :type debug: bool
    :raises ValueError: raise if no FASTA file found in `genome`
    :return: chromosomes FASTA
    :rtype: List[str]
    """
    # try files with extension FA
    chroms_fasta = glob(os.path.join(genome, "*.fa"))
    if not chroms_fasta:  # try extension FASTA
        chroms_fasta = glob(os.path.join(genome, "*.fasta"))
    if not chroms_fasta:  # both extension didn't work
        exception_handler(
            ValueError,
            f"No FASTA file found in {genome} (extensions FA or FASTA)",
            debug,
        )
    chromosomes = [
        os.path.splitext(os.path.basename(chrom_fasta))[0]
        for chrom_fasta in chroms_fasta
    ]
    write_verbosity(
        f"Chromosomes found in {genome}: {','.join(chromosomes)}", verbosity, 2, debug
    )
    return chromosomes


def recover_vcf(vcf_dir: str, verbosity: bool, debug: bool) -> List[str]:
    """Recover VCF files in vcf_dir (assumes GZipped VCFs)

    :param vcf_dir: VCFs folder
    :type vcf_dir: str
    :param verbosity: verbosity level
    :type verbosity: bool
    :param debug: debug mode
    :type debug: bool
    :raises ValueError: raise if no VCF file found in vcf_dir
    :return: VCF files
    :rtype: List[str]
    """
    write_verbosity(
        f"Reading VCFs associated to each chromosome in {vcf_dir}", verbosity, 2, debug
    )
    chrom_vcfs = glob(os.path.join(vcf_dir, "*.vcf.gz"))  # extension VCF (GZipped)
    if not chrom_vcfs:  # no VCF found
        raise ValueError
    vcfs = [os.path.basename(chrom_vcf) for chrom_vcf in chrom_vcfs]
    write_verbosity(f"VCFs found in {vcf_dir}: {','.join(vcfs)}", verbosity, 2, debug)
    return vcfs


def recover_chrom_variants(vcfs: List[str], chromosomes: List[str]) -> List[str]:
    """Recover chromosome with FASTA and VCF files

    :param vcfs: VCF files
    :type vcfs: List[str]
    :param chromosomes: chromsome FASTA files
    :type chromosomes: List[str]
    :return: chromsomes with FASTA and VCF files
    :rtype: List[str]
    """
    chromosome_variants = [
        "fakechr" + re.findall(r"chr(.+?)\.", vcf)[0] for vcf in vcfs
    ]
    assert len(chromosome_variants) == len(vcfs)
    missing = list(
        set(chromosomes).difference(set([cv[4:] for cv in chromosome_variants]))
    )
    if missing:
        raise_warning(f"missing VCF files for {','.join(missing)}")
    return chromosome_variants


def enrich_genome(
    genome: str,
    vcf: str,
    gname: str,
    vname: str,
    pam: str,
    bmax: int,
    pamfile: str,
    threads: int,
    verbosity: bool,
    debug: bool,
) -> None:
    """Enrich the genome with SNVs and indels. Indels genome is indexed to
    speed-up the off-targets search

    :param genome: genome directory
    :type genome: str
    :param vcf: VCF directory
    :type vcf: str
    :param gname: genome name
    :type gname: str
    :param vname: VCF name
    :type vname: str
    :param pam: PAM
    :type pam: str
    :param bmax: maximum bulge value
    :type bmax: int
    :param pamfile: PAM file
    :type pamfile: str
    :param threads: threads
    :type threads: int
    :param verbosity: verbosity level
    :type verbosity: bool
    :param debug: debug mode
    :type debug: bool
    """
    variants_genome = f"{gname}+{vname}"  # variants genome name
    # if enriched genome not available, start enrichment
    # if True:
    if not os.path.isdir(variants_genome):
        write_verbosity(
            f"Started adding variants to the genome at {ctime()}", verbosity, 0, debug
        )
        start_addvariants = time_verbosity(verbosity, 1, debug)
        # add variants using crispritz
        crispritz_add_variants(vcf, genome, verbosity, debug)  # TODO: call in parallel
        write_verbosity(f"Adding variants completed at {ctime()}", verbosity, 0, debug)
        stop_addvariants = time_verbosity(verbosity, 1, debug)
        write_verbosity(
            f"Adding variants took {(stop_addvariants - start_addvariants):.2f}s",
            verbosity,
            1,
            debug,
        )
        # change name to variants genome
        source = os.path.join(GENOME_SNPS, f"{gname}_enriched")
        move(source, variants_genome, debug)
        # store snp dictionary
        variants_dictionary = os.path.join(CRISPRME_DIRS[2], f"dictionaries_{vname}")
        if not os.path.isdir(CRISPRME_DIRS[2]):  # Dictionaries folder
            os.mkdir(CRISPRME_DIRS[2])
        if not os.path.isdir(variants_dictionary):
            os.mkdir(variants_dictionary)
        move(os.path.join(GENOME_SNPS, "*.json"), variants_dictionary, debug)
        # store indels dictionaries
        assert os.path.isdir(CRISPRME_DIRS[2])
        indels_dictionaries = os.path.join(CRISPRME_DIRS[2], f"log_indels_{vname}")
        if not os.path.isdir(indels_dictionaries):
            os.mkdir(indels_dictionaries)
        move(os.path.join(GENOME_SNPS, "log*.txt"), indels_dictionaries, debug)
        # indels indexing
        write_verbosity(f"Start indels indexing at {ctime()}", verbosity, 0, debug)
        start_indel_index = time_verbosity(verbosity, 1, debug)
        if not os.path.isdir(GENOMELIB):
            os.mkdir(GENOMELIB)
        indels_folder = os.path.join(GENOMELIB, f"{pam}_{bmax}_{gname}+{vname}_INDELS")
        if not os.path.isdir(indels_folder):
            os.mkdir(indels_folder)
        # index indels using crispritz
        index_indels(
            ENRICHED_GENOME, pamfile, pam, gname, vname, bmax, threads, verbosity, debug
        )
        write_verbosity(f"Indels indexing completed at {ctime()}", verbosity, 0, debug)
        stop_indel_index = time_verbosity(verbosity, 1, debug)
        write_verbosity(
            f"Indels indexing took {(stop_indel_index - start_indel_index):.2f}s",
            verbosity,
            1,
            debug,
        )
        # store indels genome
        indels_genome = f"{variants_genome}_INDELS"
        if not os.path.isdir(indels_genome):
            os.mkdir(indels_genome)
        move(os.path.join(ENRICHED_GENOME, "fake*"), indels_genome, debug)
        remove_dir(ENRICHED_GENOME, debug)  # delete enriched genome original folder
    else:
        write_verbosity(
            f"Found genome {gname} enriched with variants {vname}", verbosity, 1, debug
        )
        raise_warning("skipping genome enrichment with variants")


def index_genome(
    pam: str,
    bmax: int,
    gname: str,
    genome: str,
    pamfile: str,
    threads: int,
    verbosity: int,
    debug: bool,
    vname: Optional[str] = "",
) -> str:
    """Build the TST index for the input genome

    :param pam: PAM
    :type pam: str
    :param bmax: maximum bulge
    :type bmax: int
    :param gname: genome name
    :type gname: str
    :param genome: genome directory
    :type genome: str
    :param pamfile: PAM file
    :type pamfile: str
    :param threads: threads
    :type threads: int
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    :param vname: VCF dataset name, defaults to ""
    :type vname: Optional[str], optional
    :return: TST index directory
    :rtype: str
    """
    index_folder = os.path.join(GENOMELIB, f"{pam}_{bmax}_{gname}")
    if vname:  # variants genome
        assert isinstance(vname, str)
        index_folder = f"{index_folder}+{vname}"
    if not os.path.isdir(index_folder):  # start indexing reference/snv genome
        message = "Start indexing {} genome at {}"
        message = (
            message.format("variants", ctime())
            if vname
            else message.format("reference", ctime())
        )
        write_verbosity(message, verbosity, 0, debug)
        start_index_genome = time_verbosity(verbosity, 1, debug)
        genome_name = f"{gname}+{vname}" if vname else gname
        genome_dir = f"{gname}+{vname}" if vname else genome
        # index snv/reference genome
        crispritz_index_genome(
            genome_name, genome_dir, pamfile, bmax, threads, verbosity, debug
        )
        message = "Indexing {} genome completed at {}"
        message = (
            message.format("variants", ctime())
            if vname
            else message.format("reference", ctime())
        )
        write_verbosity(message, verbosity, 0, debug)
        stop_index_genome = time_verbosity(verbosity, 1, debug)
        message = "Indexing {} genome took {:.2f}s"
        elapsed_time = stop_index_genome - start_index_genome
        message = (
            message.format("variants", elapsed_time)
            if vname
            else message.format("reference", elapsed_time)
        )
        write_verbosity(message, verbosity, 1, debug)
    # check if indexing worked
    if len(os.listdir(index_folder)) == 0:
        exception_handler(
            TSTError, f"An error occurred during TST indexing in {index_folder}", debug
        )
    return index_folder


def run_complete_search(
    genome: str,
    vcf: str,
    pam: str,
    pam_file: str,
    bmax: int,
    guides: List[str],
    output: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """_summary_

    :param genome: _description_
    :type genome: str
    :param vcf: _description_
    :type vcf: str
    :param guides: _description_
    :type guides: List[str]
    :param output: _description_
    :type output: str
    :param verbosity: _description_
    :type verbosity: int
    :param debug: _description_
    :type debug: bool
    """
    # TODO: base editor check in resultintegration
    # open the log file
    try:
        logfile = open(os.path.join(output, LOG), mode="w")
    except OSError:
        exception_handler(
            OSError,
            f"An error occurred while writing to {os.path.join(output, LOG)}",
            debug,
        )
    write_verbosity(f"Job start {ctime()}", verbosity, 0, debug)
    # TODO: remove queue file
    # parse VCF file
    write_verbosity(f"Opening VCF list file {vcf}", verbosity, 2, debug)
    # recover chromosomes FASTA files
    chromosomes = recover_fasta(genome, verbosity, debug)
    # recover genome directory basename
    gname = (
        os.path.basename(genome[:-1])
        if genome.endswith("/")
        else os.path.basename(genome)
    )
    try:
        with open(vcf, mode="r") as vcfinfile:
            start_vcf_tasks = time_verbosity(verbosity, 2, debug)
            write_verbosity(
                f"Start tasks for all datasets in {vcf}", verbosity, 2, debug
            )
            for line in vcfinfile:
                vcf_dataset = line.strip()
                if not vcf_dataset:  # skip empty lines
                    continue
                if not os.path.isdir(vcf_dataset):  # unable to find the VCF directory
                    exception_handler(
                        FileNotFoundError,
                        f"Unable to locate VCFs directory {vcf_dataset}",
                        debug,
                    )
                # recover VCF directory basename
                vname = (
                    os.path.basename(vcf_dataset[:-1])
                    if vcf_dataset.endswith("/")
                    else os.path.basename(vcf_dataset)
                )
                write_verbosity(
                    f"Starting analyis on {vcf_dataset}", verbosity, 0, debug
                )
                # recover VCFs and chromosomes with both FASTA and VCF files
                chrom_vcfs = recover_vcf(vcf_dataset, verbosity, debug)
                chromosomes_variants = recover_chrom_variants(chrom_vcfs, chromosomes)
                # TODO: create directories structure for web-site?
                # start genome snp and indels enrichment
                enrich_genome(
                    genome,
                    vcf_dataset,
                    gname,
                    vname,
                    pam,
                    bmax,
                    pam_file,
                    threads,
                    verbosity,
                    debug,
                )
                # index reference genome (guarenteed to be indexed only once)
                tst_reference = index_genome(
                    pam, bmax, gname, genome, pam_file, threads, verbosity, debug
                )
                # index snv genome
                tst_variants = index_genome(
                    pam,
                    bmax,
                    gname,
                    genome,
                    pam_file,
                    threads,
                    verbosity,
                    debug,
                    vname=vname,
                )

    except OSError:
        exception_handler(OSError, f"An error occurred while reading {vcf}", debug)
    stop_vcf_tasks = time_verbosity(verbosity, 2, debug)
    write_verbosity(
        f"Ended tasks for all datasets in {vcf} in {(stop_vcf_tasks - start_vcf_tasks):.2f}s",
        verbosity,
        2,
        debug,
    )
