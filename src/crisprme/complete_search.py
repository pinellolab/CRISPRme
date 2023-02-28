"""
"""

from crisprme_errors import TSTError
from crispritz import (
    ENRICHED_GENOME,
    GENOME_SNPS,
    crispritz_add_variants,
    crispritz_index_genome,
    crispritz_search,
)
from indels_utils import index_indels, search_indels
from sequence import write_guidefile
from verbosity_handler import write_verbosity, time_verbosity
from utils import (
    CRISPRME_DIRS,
    LOG,
    exception_handler,
    move,
    raise_warning,
    remove_dir,
    remove
)

from typing import List, Optional, Tuple
from glob import glob
from time import ctime, time

import sys
import re
import os


GENOMESENR = "genome_enriched"  # store enriched genomes
GENOMELIB = "genome_library"  # store TST indexes


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
    if missing := list(
        set(chromosomes).difference({cv[4:] for cv in chromosome_variants})
    ):
        raise_warning(f"missing VCF files for {','.join(missing)}") 
    return chromosome_variants


def recover_tst_index(genome_index: str, pam: str, bmax: int, gname: str, verbosity: int, debug: bool) -> Tuple[str, str]:
    """Recover reference and variant genomes TST indexes

    :param genome_index: directory storing the TST indexes
    :type genome_index: str
    :param pam: PAM
    :type pam: str
    :param bmax: max bulge value
    :type bmax: int
    :param gname: genome name
    :type gname: str
    :param verbosity: verbosity level
    :type verbosity: int
    :param debug: debug mode
    :type debug: bool
    :raises ValueError: raise on reference TST and input PAM mismatch
    :raises ValueError: raise on reference TST and input bmax mismatch
    :raises ValueError: raise on reference TST and input gname mismatch
    :raises ValueError: raise on variant TST and input PAM mismatch
    :raises ValueError: raise on variant TST and input bmax mismatch
    :raises ValueError: raise on variant TST and input gname mismatch 
    :return: reference and variant genome (if available) TST indexes
    :rtype: Tuple[str, str]
    """
    write_verbosity(f"Recovering TST index in {genome_index}", verbosity, 2, debug)
    tstindexes = os.listdir(genome_index)  # recover TST indexes
    # recover the reference TST
    tst_reference = [idx for idx in tstindexes if "+" not in idx]
    assert len(tst_reference) == 1  # should be only one
    tst_reference = tst_reference[0]
    # recover reference TST PAM, bmax and genome name
    pam_r, bmax_r, gname_r = tst_reference.split("_", 2)
    # check consistency with other input args
    if pam != pam_r:
        exception_handler(ValueError, f"Reference TST ({pam_r}) and input PAM ({pam}) mismatch", debug)
    if bmax != int(bmax_r):
        exception_handler(ValueError, f"Reference TST ({bmax_r}) input max bulge ({bmax}) mismatch", debug)
    if gname != gname_r:
        exception_handler(ValueError, f"Reference TST ({gname_r}) and input genome name ({gname}) mismatch", debug)
    tst_reference = os.path.join(genome_index, tst_reference)
    write_verbosity(f"Found reference genome TST: {tst_reference}", verbosity, 1, debug)
    # recover variants TST
    tst_variants = [idx for idx in tstindexes if "+" in idx and not idx.endswith("INDELS")]
    if tst_variants:
        assert len(tst_variants) == 1  # should be only one
        tst_variants = tst_variants[0]
        # recover variants PAM, bmax and genome name
        pam_v, bmax_v, gname_v = tst_variants.split("_", 2)
        gname_v = gname_v.split("+")[0]
        if pam != pam_v:
            exception_handler(ValueError, f"Variants TST ({pam_v}) and input PAM ({pam}) mismatch", debug)
        if bmax != int(bmax_v):
            exception_handler(ValueError, f"Variants TST ({bmax_v}) input max bulge ({bmax}) mismatch", debug)
        if gname != gname_v:
            exception_handler(ValueError, f"Variants TST ({gname_v}) and input genome name ({gname}) mismatch", debug)
        tst_variants = os.path.join(genome_index, tst_variants)
        write_verbosity(f"Found variants genome TST: {tst_variants}", verbosity, 1, debug)
    else:
        tst_variants = ""
    return tst_reference, tst_variants


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
    # check enriched genomes directory existence
    if not os.path.isdir(GENOMESENR):
        os.mkdir(GENOMESENR)
    variants_genome_snv = os.path.join(GENOMESENR, variants_genome)
    variants_genome_indels = os.path.join(GENOMESENR, f"{variants_genome}_INDELS")
    # if enriched genome not available, start enrichment
    if not os.path.isdir(variants_genome_snv) and not os.path.isdir(variants_genome_indels):
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
        move(source, variants_genome_snv, debug)
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
        move(os.path.join(ENRICHED_GENOME, "fake*"), variants_genome_indels, debug)
        remove_dir(ENRICHED_GENOME, debug)  # delete enriched genome original folder
    else:
        write_verbosity(
            f"Found genome {gname} enriched with variants {vname}", verbosity, 1, debug
        )
        raise_warning(f"skipping genome enrichment with variants because {variants_genome_snv} and {variants_genome_indels} were both found")


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
        genome_dir = os.path.join(GENOMESENR, f"{gname}+{vname}") if vname else genome
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
    else:
        raise_warning(f"skipping genome indexing because {index_folder} was found")
    # check if indexing worked
    if len(os.listdir(index_folder)) == 0:
        exception_handler(
            TSTError, f"An error occurred during TST indexing in {index_folder}", debug
        )
    return index_folder

def search(tst_index: str, gname: str, pam: str, pamfile: str, guide: str, guide_name: str, mm: int, bDNA: int, bRNA: int, threads: int, verbosity: int, debug: bool, variant: Optional[bool] = False, vname: Optional[str] = None) -> None:
    """Search potential off-targets for the input guide using the specified TST 
    index

    :param tst_index: TST index
    :type tst_index: str
    :param gname: genome name
    :type gname: str
    :param pam: PAM
    :type pam: str
    :param pamfile: PAM file
    :type pamfile: str
    :param guide: guide file
    :type guide: str
    :param guide_name: guide  
    :type guide_name: str
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
    :param variant: search variant genome, defaults to False
    :type variant: Optional[bool], optional
    :param vname: variants dataset, defaults to None
    :type vname: Optional[str], optional
    """
    if not variant:  # no variants dataset should be provided
        assert vname is None
    else:  # variants dataset should be provided
        assert (vname is not None and isinstance(vname, str) and bool(vname))
    message = "Start off-targets search on {} genome at {} for guide {}"
    message = message.format("variant", ctime(), guide_name) if variant else message.format("reference", ctime(), guide_name)
    write_verbosity(message, verbosity, 0, debug)
    start_search = time_verbosity(verbosity, 1, debug)
    reportname = f"{gname}_{pam}_{guide_name}_{mm}_{bDNA}_{bRNA}"
    if variant:
        reportname = reportname.replace(gname, f"{gname}+{vname}")
    # run crispritx search
    crispritz_search(tst_index, pamfile, guide, reportname, mm, bDNA, bRNA, threads, verbosity, debug)
    message = "Off-targets search on {} genome for guide {} completed at {}"
    message = message.format("variant", guide_name, ctime()) if variant else message.format("reference", guide_name, ctime())
    write_verbosity(message, verbosity, 0, debug)
    stop_search = time_verbosity(verbosity, 1, debug)
    write_verbosity(f"Off-targets search for {guide_name} took {(stop_search - start_search):.2f}s", verbosity, 1, debug)


def run_complete_search(
    genome: str,
    vcf: str,
    pam: str,
    pam_file: str,
    bmax: int,
    guides: List[str],
    mm: int,
    bDNA: int,
    bRNA: int,
    output: str,
    threads: int,
    verbosity: int,
    debug: bool,
    genome_index: Optional[str] = "",
) -> None:
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
    # recover genome directory basename
    gname = (
        os.path.basename(genome[:-1])
        if genome.endswith("/")
        else os.path.basename(genome)
    )
    # TODO: handle TST index and search
    if genome_index:  # precomputed TST-index as input
        tst_reference, tst_variants = recover_tst_index(genome_index, pam, bmax, gname, verbosity, debug)
    # TODO: remove queue file
    if vcf:  # parse VCF file
        write_verbosity(f"Opening VCF list file {vcf}", verbosity, 2, debug)
    # recover chromosomes FASTA files
    chromosomes = recover_fasta(genome, verbosity, debug)
    # index reference genome (indexed only once)
    tst_reference = index_genome(
        pam, bmax, gname, genome, pam_file, threads, verbosity, debug
    )
    if vcf:  # recover variants datasets to enrich the genome
        try:
            with open(vcf, mode="r") as vcfinfile:
                vcf_datasets = [line.strip() for line in vcfinfile if line.strip()]
        except OSError:
            exception_handler(OSError, f"An error occurred while reading {vcf}", debug) 
    tst_variants_dict = {}  # list of TST variants indexes
    if vcf:  # enrich reference genome with variants in vcf_datasets
        start_vcf_tasks = time_verbosity(verbosity, 2, debug)
        write_verbosity(
            f"Start tasks for all datasets in {vcf}", verbosity, 2, debug
        )
        for vcf_dataset in vcf_datasets:
            vcf_dataset = os.path.join(CRISPRME_DIRS[3], vcf_dataset) 
            # VCF dataset should be in VCFs
            if not os.path.isdir(vcf_dataset):
                exception_handler(FileNotFoundError, f"{os.path.basename(vcf_dataset)} not found in {CRISPRME_DIRS[3]}", debug)
                raise FileNotFoundError  # TODO: remove
            # recover VCF dataset name (directory name)
            vname = os.path.basename(vcf_dataset[:-1]) if vcf_dataset.endswith("/") else os.path.basename(vcf_dataset)
            write_verbosity(
                f"Starting genome enrichment using {vname}", verbosity, 0, debug
            )
            # recover VCFs and chromosomes with both FASTA and VCF files
            chrom_vcfs = recover_vcf(vcf_dataset, verbosity, debug)
            chromosomes_variants = recover_chrom_variants(chrom_vcfs, chromosomes)
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
            # index snv genome (indels already indexed)
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
            assert vname not in tst_variants_dict.keys()
            tst_variants_dict[vname] = tst_variants  # add the variants TST
        stop_vcf_tasks = time_verbosity(verbosity, 2, debug)
        write_verbosity(
            f"Ended tasks for all datasets in {vcf} in {(stop_vcf_tasks - start_vcf_tasks):.2f}s",
            verbosity,
            2,
            debug,
        )
    # search off-targets for each guide
    for guide in guides:
        # search the reference genome
        guidefile = write_guidefile(guide, verbosity, debug)  # create guide file
        search(tst_reference, gname, pam, pam_file, guidefile, guide, mm, bDNA, bRNA, threads, verbosity, debug)
        # search the variants genome
        for vname in tst_variants_dict:
            search(tst_variants_dict[vname], gname, pam, pam_file, guidefile, guide, mm, bDNA, bRNA, threads, verbosity, debug, variant=True, vname=vname)
            search_indels(f"{tst_variants_dict[vname]}_INDELS", gname, vname, pam, pam_file, mm, bDNA, bRNA, guide, guidefile, threads, verbosity, debug)
        remove(guidefile, debug)  # remove temporary guide file
        