"""
This module provides functionality for enriching genomic sequences by integrating
variant information from VCF files. It includes functions for parsing command line
arguments, retrieving chromosome data, loading FASTA and VCF files, encoding SNPs
using IUPAC notation, and saving enriched sequences along with variant and indel
dictionaries.

Key functionalities include:
- Parsing command line arguments to retrieve necessary directories and flags.
- Loading and processing genomic data from FASTA and VCF files.
- Enriching sequences by incorporating SNPs and indels.
- Saving the enriched sequences and associated variant information to specified
  directories.
"""

from variants_dictionary import (
    update_variant_dictionaries,
    save_variant_dictionary,
    save_indels_dictionary,
)
from utils import create_directory, copy

from glob import glob
from typing import Tuple, List, Dict, Union
from time import time

import shutil
import pysam
import sys
import os

IUPAC = {
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "a": "A",
    "t": "T",
    "c": "C",
    "g": "G",
    "N": "ATGC",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "r": "AG",
    "y": "CT",
    "s": "GC",
    "w": "AT",
    "k": "GT",
    "m": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "b": "CGT",
    "d": "AGT",
    "h": "ACT",
    "v": "ACG",
}
IUPAC_ENCODER = {
    "A": ["A"],
    "T": ["T"],
    "C": ["C"],
    "G": ["G"],
    "R": ["AG", "GA"],
    "Y": ["CT", "TC"],
    "S": ["GC", "CG"],
    "W": ["AT", "TA"],
    "K": ["GT", "TG"],
    "M": ["AC", "CA"],
    "B": ["CGT", "GCT", "TGC", "GTC", "CTG", "TCG"],
    "D": ["AGT", "GAT", "TAG", "ATG", "GTA", "TGA"],
    "H": ["ACT", "CAT", "TCA", "ATC", "CTA", "TAC"],
    "V": ["ACG", "CAG", "GAC", "AGC", "CGA", "GCA"],
    "N": [
        "ACGT",
        "CAGT",
        "GACT",
        "AGCT",
        "CGAT",
        "GCAT",
        "GCTA",
        "CGTA",
        "TGCA",
        "GTCA",
        "CTGA",
        "TCGA",
        "TAGC",
        "ATGC",
        "GTAC",
        "TGAC",
        "AGTC",
        "GATC",
        "CATG",
        "ACTG",
        "TCAG",
        "CTAG",
        "ATCG",
        "TACG",
    ],
}
IUPAC_ENCODER_REV = {v: k for k in IUPAC_ENCODER for v in IUPAC_ENCODER[k]}
VARIANTSGENOMEDIR = "variants_genome"
SNPSGENOMEDIR = "SNPs_genome"
INDELSGENOMEDIR = "INDELs_genome"  # NOTE only for retrocompatibility


def parse_commandline(args: List[str]) -> Tuple[str, str, bool]:
    """
    Parses command line arguments to retrieve directory paths and a boolean flag.
    This function ensures that the provided directories exist and returns their
    absolute paths along with a boolean indicating the presence of indels.

    Args:
        args (List[str]): A list of command line arguments containing two directory
            paths and a string indicating indels.

    Returns:
        Tuple[str, str, bool]: A tuple containing the absolute paths of the VCF
            directory and genome directory, and a boolean indicating if indels are
            present.

    Raises:
        ValueError: If the number of input arguments is not exactly three.
        FileNotFoundError: If the specified VCF or genome directory does not exist.
    """

    # parse commandline arguments
    if len(args) != 3:
        raise ValueError("Too many/few input arguments")
    vcfdir, genomedir, indels = args
    if not os.path.isdir(vcfdir):
        raise FileNotFoundError(f"Unable to find {vcfdir}")
    if not os.path.isdir(genomedir):
        raise FileNotFoundError(f"Unable to find {genomedir}")
    return os.path.abspath(vcfdir), os.path.abspath(genomedir), indels == "true"


def retrieve_chromosome_vcfs(vcfdir: str):
    """
    Retrieves VCF files corresponding to each chromosome from the specified directory.
    This function assumes a one-to-one mapping between VCF files and chromosomes,
    returning a dictionary that maps chromosome identifiers to their respective VCF
    file paths.

    Args:
        vcfdir (str): The directory containing VCF files.

    Returns:
        dict: A dictionary where keys are chromosome identifiers and values are the
            corresponding VCF file paths.
    """

    # recover VCF file for each chromosome within the input VCF directory
    # assume 1:1 mapping between VCFs and chromosomes
    chroms_vcfs = [
        f for f in glob(os.path.join(vcfdir, "*.vcf.gz")) if os.path.isfile(f)
    ]
    return {e: f for f in chroms_vcfs for e in f.split(".") if "chr" in e}


def retrieve_chromosome_fasta(genomedir: str):
    """
    Retrieves chromosome FASTA files from the specified genome directory.
    This function collects both `.fa` and `.fasta` files, returning a dictionary
    that maps chromosome names to their corresponding file paths, while handling
    enriched FASTA files appropriately.

    Args:
        genomedir (str): The directory containing chromosome FASTA files.

    Returns:
        dict: A dictionary where keys are chromosome names (without extensions)
            and values are the corresponding FASTA file paths.
    """

    # recover chromosome fasta files
    chroms_fasta = [
        f
        for f in glob(os.path.join(genomedir, "*.fa"))
        + glob(os.path.join(genomedir, "*.fasta"))
        if os.path.isfile(f)
    ]
    # handle cases where the enriched fasta has been put in the genome folder
    return {
        os.path.splitext(os.path.basename(chrom_fasta))[0].replace(
            ".enriched", ""
        ): chrom_fasta
        for chrom_fasta in chroms_fasta
    }


def load_fasta(fastafile: str) -> Tuple[str, List[str]]:
    """
    Loads a FASTA file and extracts the contig name and nucleotide sequence. This
    function reads the file, processes the sequence into a list for efficient
    nucleotide replacement, and handles any file loading errors.

    Args:
        fastafile (str): The path to the input FASTA file.

    Returns:
        Tuple[str, List[str]]: A tuple containing the contig name and a list of
            nucleotides in the sequence.

    Raises:
        OSError: If there is an error loading the specified FASTA file.
    """

    try:  # read the input FASTA file
        with open(fastafile, mode="r") as infile:
            contig = infile.readline().strip()[1:]
            fasta = infile.read().replace("\n", "").upper()
    except IOError as e:
        raise OSError(f"{fastafile} loading failed") from e
    # explode the sequence in list to allow efficient nt replacement during
    # enrichment -> allows to avoid calls to pysam.FastaFile (not efficient on
    # repetitve calls)
    return contig, list(fasta)


def index_vcf(vcffile: str) -> str:
    """
    Index a VCF file with tabix for faster parsing and access.

    Args:
        vcffile (str): The path to the VCF file to index.
        debug (bool): A flag indicating debug mode.

    Returns:
        str: The path to the generated index file.
    Raises:
        FileExistsError: If indexing the VCF file fails.
    """

    # vcf indexing enables faster vcf parsing and access
    sys.stdout.write(f"indexing {vcffile}\n")
    pysam.tabix_index(vcffile, preset="vcf")  # index input vcf with tabix
    tbi = f"{vcffile}.tbi"
    if not os.path.isfile(tbi) or os.stat(tbi).st_size <= 0:
        raise FileExistsError(f"Indexing {vcffile} failed")
    return tbi


def load_vcf(vcffile: str) -> Tuple[pysam.TabixFile, List[str]]:
    """
    Loads a VCF file using a Tabix index for efficient parsing. This function retrieves
    the VCF data and the associated sample names, handling the creation of the index
    if it does not already exist.

    Args:
        vcffile (str): The path to the input VCF file.

    Returns:
        Tuple[pysam.TabixFile, List[str]]: A tuple containing the loaded VCF as a
            TabixFile object and a list of sample names extracted from the VCF
            header.

    Raises:
        RuntimeError: If there is an error loading the specified VCF file.
    """

    # using tabix index allows efficient parsing of the VCF without requiring
    # file decompression and saturating memory while running enrichment in
    # parallel
    tbi = f"{vcffile}.tbi"
    # search for tabix index, if not found compute
    tbi_index = tbi if os.path.isfile(tbi) else index_vcf(vcffile)
    try:
        vcf = pysam.TabixFile(vcffile, index=tbi_index)  # load VCF as TabixFile
        # recover vcf samples (from 10th col in VCF)
        samples = vcf.header[-1].strip().split()[9:]
        return vcf, samples
    except IOError as e:
        raise RuntimeError(f"{vcffile} loading failed") from e


def encode_snp_iupac(
    fasta: List[str], ref: str, alt: str, pos: int
) -> Union[str, None]:
    """
    Encodes a single nucleotide polymorphism (SNP) using IUPAC notation based on
    the reference and alternate alleles.
    This function checks the validity of the reference and alternate alleles and
    returns the corresponding IUPAC symbol if the reference allele matches the FASTA
    sequence at the specified position.

    Args:
        fasta (List[str]): The FASTA sequence represented as a list of characters.
        ref (str): The reference allele.
        alt (str): The alternate allele(s) as a comma-separated string.
        pos (int): The position of the SNP in the FASTA sequence (1-based index).

    Returns:
        Union[str, None]: The IUPAC encoded string representing the SNP, or None
            if the input is invalid.

    Raises:
        ValueError: If there is a mismatch between the reference allele and the
            FASTA sequence at the specified position.
    """

    if len(ref) != 1:  # indel, likely  deletion
        return
    alleles_alt = alt.split(",")  # handle multiallelic sites
    if len(alleles_alt) == 1 and len(alleles_alt[0]) > 1:
        return
    refnt = fasta[pos - 1]  # 0-based position
    if refnt != ref:
        raise ValueError(
            f"Reference allele mismatch between FASTA and VCF file ({refnt} - {ref}, position {pos})"
        )
    iupac_string = "".join(
        set([refnt] + [aa for aa in alleles_alt if len(aa) == 1])
    )  # recover iupac symbol
    return IUPAC_ENCODER_REV[iupac_string]


def insert_snp(
    fasta: List[str], sequence_enr: List[str], snppos: int, ref: str, alt: str
) -> List[str]:
    """
    Inserts a single nucleotide polymorphism (SNP) into an enriched sequence based
    on the reference and alternate alleles.
    This function encodes the SNP using IUPAC notation and updates the sequence
    at the specified position if the SNP is valid.

    Args:
        fasta (List[str]): The FASTA sequence represented as a list of characters.
        sequence_enr (List[str]): The enriched sequence represented as a list of
            characters where the SNP will be inserted.
        snppos (int): The position in the sequence to insert the SNP (1-based index).
        ref (str): The reference allele.
        alt (str): The alternate allele.

    Returns:
        List[str]: The updated enriched sequence with the SNP inserted, or the
            original sequence if the SNP is invalid.
    """

    # recover iupac character to encode variants
    iupac_nt = encode_snp_iupac(fasta, ref, alt, snppos)
    if iupac_nt is None:  # likely indel
        return sequence_enr
    sequence_enr[snppos - 1] = iupac_nt
    return sequence_enr


def save_enriched_sequence(sequence_enr: str, header: str, outdir: str) -> None:
    """
    Saves an enriched nucleotide sequence to a specified output directory in FASTA
    format.

    Args:
        sequence_enr (str): The enriched nucleotide sequence to be saved.
        header (str): The header to be used for the FASTA file.
        outdir (str): The base directory where the enriched sequence will be saved.

    Raises:
        OSError: If there is an error writing the enriched sequence to the specified
            file.
    """

    fasta_enriched = os.path.join(outdir, f"{header}.enriched.fa")  # define output file
    try:  # write enriched sequence to fasta file
        with open(fasta_enriched, mode="w") as outfile:
            outfile.write(f">{header}\n{sequence_enr}\n")
    except IOError as e:
        raise OSError(f"Failed writing enriched sequence to {fasta_enriched}") from e


def enrich_sequence(
    fasta: List[str],
    contig: str,
    vcf: pysam.TabixFile,
    samples: List[str],
) -> Tuple[str, Dict[str, str], Dict[str, Tuple[str, List[str]]]]:
    """
    Enriches a nucleotide sequence by incorporating variants from a VCF file. This
    function updates the sequence using IUPAC notation for SNPs and maintains
    dictionaries to link samples to their respective variants and indels.

    Args:
        fasta (List[str]): The genomic sequence in FASTA format.
        contig (str): The name of the contig corresponding to the sequence.
        vcf (pysam.TabixFile): A Tabix file object containing variant information
            in VCF format.
        samples (List[str]): A list of sample identifiers.

    Returns:
        Tuple[str, Dict[str, str], Dict[str, Tuple[str, List[str]]]]: A tuple
            containing the enriched nucleotide sequence, a dictionary of variants,
            and a dictionary of indels.

    Raises:
        AssertionError: If the chromosome name in the VCF does not match the expected
            contig format.
    """

    # create the variant dictionary, used to link samples to the carried variants
    variants_dictionary = {}
    # create indels dictionary, used to link samples to the carried indels and
    # set local indel id (1-based)
    indels_dictionary, indelid = {}, 1
    # initialize the enriched sequence using iupac chars to denote variants
    sequence_enr = fasta.copy()
    # set sequence start position (0-based), indel start position is the same
    indelstart = 0
    for variant in vcf.fetch():
        variant = variant.strip().split()  # recover vcf fields for each variant
        chrom, snppos, snpid, ref, alt, _, vfilter, info = variant[:8]  # ignore FORMAT
        snppos = int(snppos)  # position is read as string
        # chrom must contain chr
        assert contig == (chrom := chrom if chrom[:3] == "chr" else f"chr{chrom}")
        if vfilter == "PASS":  # skip variants not labeled as pass
            genotypes = variant[9:]  # start samples idx
            if ":" in genotypes[0]:  # support format field with multiple fields
                genotypes = [gt.split(":")[0] for gt in genotypes]
            variants_dictionary, indels_dictionary, indelid, indelstart = (
                update_variant_dictionaries(
                    variants_dictionary,
                    indels_dictionary,
                    fasta,
                    chrom,
                    snppos,
                    snpid,
                    ref,
                    alt,
                    info,
                    genotypes,
                    samples,
                    indelstart,
                    indelid,
                )
            )  # update variants and indels dictionaries entry
            sequence_enr = insert_snp(
                fasta, sequence_enr, snppos, ref, alt
            )  # enrich the sequence
    return "".join(sequence_enr), variants_dictionary, indels_dictionary


def enrich(
    fastafile: str, vcffile: str, genomedir: str, vcfdir: str, outdir: str
) -> None:
    """
    Enriches a genomic sequence by integrating variant information from a VCF file.
    This function loads the FASTA and VCF data, processes the variants to update
    the sequence, and saves the enriched sequence along with the variant and indel
    dictionaries.

    Args:
        fastafile (str): The path to the input FASTA file.
        vcffile (str): The path to the input VCF file.
        out_seq (str): The directory where the enriched sequence will be saved.
        out_vdict (str): The directory where the snp dictionary will be saved.
        out_idict (str): The directory where the indel dictionary will be saved.
        debug (bool): A flag indicating whether to run in debug mode, which may
            include additional logging or checks.

    Returns:
        None
    """

    contig, fasta = load_fasta(fastafile)  # load input fasta
    vcf, samples = load_vcf(vcffile)  # load input vcf data
    sequence_enr, variants_dictionary, indels_dictionary = enrich_sequence(
        fasta, contig, vcf, samples
    )
    # save enriched sequence and variants dictionaries
    snps_genome_dir = create_directory(os.path.join(outdir, SNPSGENOMEDIR))
    enriched_genomedir = create_directory(
        os.path.join(snps_genome_dir, f"{os.path.basename(genomedir)}_enriched")
    )
    save_enriched_sequence(sequence_enr, contig, enriched_genomedir)
    # save variant dictionary in SNPs_genome
    save_variant_dictionary(variants_dictionary, contig, snps_genome_dir)
    # save indels fasta in fake_<vcf_folder>_<contig> and indels log in SNPs_genome
    indels_genome_dir = create_directory(
        os.path.join(outdir, INDELSGENOMEDIR)
    )  # NOTE only for retrocompatibility
    save_indels_dictionary(indels_dictionary, contig, vcfdir, outdir, snps_genome_dir)


def main():
    # parse command line arguments
    vcfdir, genomedir, indels = parse_commandline(sys.argv[1:])
    # recover chromosomes with associated VCFs
    chroms_vcf = retrieve_chromosome_vcfs(vcfdir)
    chroms_fasta = retrieve_chromosome_fasta(genomedir)
    chroms_wo_vcf = set(chroms_fasta.keys()) - set(
        chroms_vcf.keys()
    )  # recover chromosomes without vcf
    # create variants_genome directory -> contain the enriched genome sequence
    # snps enriched put in SNPs_genome
    # indels put in fake_<outdir>_<chrom> (INDELs_genome just for retrocompatibility)
    variants_genome_dir = create_directory(VARIANTSGENOMEDIR)
    for chrom, vcf in chroms_vcf.items():
        start = time()
        sys.stdout.write(f"Start genome enrichment with SNVs and indels on {chrom}\n")
        enrich(chroms_fasta[chrom], vcf, genomedir, vcfdir, variants_genome_dir)
        sys.stdout.write(
            f"Genome enrichment performed on {chrom} in {(time() - start):.2f}s\n"
        )
    # copy chromosomes without VCF in the enriched genome folder
    enriched_genomedir = os.path.join(
        variants_genome_dir, SNPSGENOMEDIR, f"{os.path.basename(genomedir)}_enriched"
    )
    assert os.path.isdir(enriched_genomedir)
    for chrom in chroms_wo_vcf:
        copy(
            chroms_fasta[chrom],
            os.path.join(enriched_genomedir, f"{chrom}.enriched.fa"),
        )


if __name__ == "__main__":
    main()
