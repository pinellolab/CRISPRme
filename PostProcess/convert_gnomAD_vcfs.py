"""
Module for converting gnomAD VCF files.

This module provides functionality to parse command line arguments, read sample IDs,
and convert gnomAD VCF files by updating their headers, filtering variants, and merging
alleles. It utilizes multiprocessing to handle multiple VCF files efficiently and
ensures robust error handling throughout the conversion process.

Key functions include:
- `parse_commandline`: Validates and parses command line arguments.
- `read_samples_ids`: Reads sample IDs from a specified file.
- `tabix_index`: Creates an index for a VCF file using tabix.
- `load_vcf`: Loads a VCF file and indexes it if necessary.
- `update_header`: Updates the VCF header with sample information.
- `variant_observed`: Checks if any allele count indicates a variant is observed.
- `format_variant_record`: Formats a variant record into a string.
- `convert_vcf`: Converts a VCF file by updating its header and filtering variants.
- `bcftools_merge`: Merges alleles in a VCF file using bcftools.
- `run_conversion_pipeline`: Runs the conversion pipeline for a single VCF file.
- `convert_gnomad_vcfs`: Main entry point for converting gnomAD VCF files based 
    on user parameters.

This module is designed to facilitate the processing of genomic data for research 
and analysis purposes.
"""

from utils import remove

from functools import partial
from typing import List, Tuple
from glob import glob
from io import TextIOWrapper

import multiprocessing
import subprocess
import itertools
import pysam
import gzip
import time
import sys
import os

GT = ["0/0", "0/1"]
GTLINE = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Sample Collapsed Genotype">'
BCFTOOLSNORM = "bcftools norm"


def parse_commandline(args: List[str]) -> Tuple[str, str, bool, bool, bool, int, frozenset, float, str]:
    """Parse and validate command line arguments for gnomAD VCF conversion.

    This function checks the provided arguments for correctness and extracts the
    necessary parameters for the gnomAD VCF conversion process. It ensures that
    the directory, joint processing flag, and thread count are valid before
    returning the parsed values.

    Positional args (6 required, 7th optional):
        gnomad_vcfs_dir samples_ids joint keep multiallelic threads [filter_pass_values]

    Optional keyword args (may appear anywhere after the positional args):
        --af-threshold FLOAT   Only include variants with AF_joint >= this value
                               (default: 0.0, i.e. no filtering).
        --output-dir PATH      Directory for converted output VCFs. Created if
                               absent. Default: same directory as each input VCF.

    Args:
        args (List[str]): A list of command line arguments.

    Returns:
        Tuple containing: gnomAD VCF directory, sample IDs, joint flag, keep flag,
        multiallelic flag, threads, filter_pass_values frozenset, af_threshold float,
        and output_dir string (empty string means same dir as input).

    Raises:
        ValueError: If the number of arguments is incorrect, if the specified
            directory is not valid, or if the number of threads is out of
            allowed range.
    """

    # Extract optional keyword args before counting positional args
    af_threshold = 0.0
    output_dir = ""
    positional = []
    i = 0
    while i < len(args):
        if args[i] == "--af-threshold":
            if i + 1 >= len(args):
                raise ValueError("Missing value for --af-threshold")
            af_threshold = float(args[i + 1])
            if not (0.0 <= af_threshold < 1.0):
                raise ValueError(f"--af-threshold must be in [0, 1), got {af_threshold}")
            i += 2
        elif args[i] == "--output-dir":
            if i + 1 >= len(args):
                raise ValueError("Missing value for --output-dir")
            output_dir = args[i + 1]
            i += 2
        else:
            positional.append(args[i])
            i += 1

    if len(positional) not in (6, 7):
        raise ValueError(
            "Wrong number of input arguments, cannot proceed with gnomAD VCF conversion"
        )
    gnomad_vcfs_dir, samples_ids, joint, keep, multiallelic, threads = positional[:6]
    filter_pass_values = frozenset(positional[6].split(",")) if len(positional) == 7 else frozenset({"PASS", "."})
    if not os.path.isdir(gnomad_vcfs_dir):
        raise ValueError(
            f"The specified gnomAD VCF directory is not a directory ({gnomad_vcfs_dir})"
        )
    threads = int(threads)
    if threads > multiprocessing.cpu_count() or threads < 0:
        raise ValueError(f"Forbidden number of threads selected ({threads})")
    joint = joint == "True"
    keep = keep == "True"
    multiallelic = multiallelic == "True"
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    return gnomad_vcfs_dir, samples_ids, joint, keep, multiallelic, threads, filter_pass_values, af_threshold, output_dir


def read_samples_ids(samples_ids: str):
    """
    Reads sample IDs from a file and returns a list of sample IDs excluding comments.

    Args:
        samples_ids (str): Path to the file containing sample IDs.

    Returns:
        list: List of sample IDs extracted from the file.

    Raises:
        FileNotFoundError: If the specified samples file is not found.
        IOError: If an error occurs while parsing the samples IDs file.
    """

    if not os.path.isfile(samples_ids):
        raise FileNotFoundError(f"Unable to locate {samples_ids}")
    try:
        with open(samples_ids, mode="r") as infile:  # parse sample ids file
            samplesdata = [line.strip() for line in infile]
    except IOError as e:
        raise IOError(
            f"An error occurred while parsing samples IDs file {samples_ids}"
        ) from e
    # create samples list
    return [sample.split()[0] for sample in samplesdata if not sample.startswith("#")]


def tabix_index(vcf_fname: str) -> str:
    """
    Creates an index for a VCF file using tabix and returns the path to the
    generated index file.

    Args:
        vcf_fname (str): Path to the VCF file to be indexed.

    Returns:
        str: Path to the generated tabix index file.

    Raises:
        FileExistsError: If indexing of the VCF file fails.
    """

    start = time.time()
    pysam.tabix_index(vcf_fname, preset="vcf")  # index input vcf with tabix
    tbi_index = f"{vcf_fname}.tbi"
    if not os.path.isfile(tbi_index) and os.stat(tbi_index).st_size <= 0:
        raise FileExistsError(f"Indexing {vcf_fname} failed")
    sys.stderr.write(
        f"Indexing {vcf_fname} completed in {(time.time() - start):.2f}s\n"
    )
    return tbi_index


def load_vcf(vcf_fname: str) -> pysam.VariantFile:
    """
    Loads a VCF file from the specified path, automatically indexing it if the
    tabix index is not found.

    Args:
        vcf_fname (str): Path to the VCF file to load.

    Returns:
        pysam.VariantFile: Loaded VCF file object.
    """

    # search for tabix index, if not found index the vcf
    tbi_index = (
        f"{vcf_fname}.tbi"
        if os.path.isfile(f"{vcf_fname}.tbi")
        else tabix_index(vcf_fname)
    )
    return pysam.VariantFile(vcf_fname, index_filename=tbi_index)


def update_header(header: pysam.VariantHeader, samples: List[str], joint: bool) -> str:
    """
    Updates the header of a VCF file with the specified samples and additional
    metadata fields.

    Args:
        header (pysam.VariantHeader): The original VCF header to be updated.
        samples (List[str]): List of sample names to be added to the header.

    Returns:
        str: Updated VCF header as a string.
    """

    header.add_line(GTLINE)  # add FORMAT metadata field
    header.add_samples(samples)  # add samples to header
    header = str(header).replace("<ID=AF_joint,", "<ID=AF,") if joint else str(header)  # type: ignore
    return header  # type: ignore


def variant_observed(allele_count: Tuple[int]) -> bool:
    """
    Checks if any allele count in the provided tuple is greater than zero, indicating
    the variant is observed.

    Args:
        allele_count (Tuple[int]): Tuple of allele counts for the variant.

    Returns:
        bool: True if any allele count is greater than zero, False otherwise.
    """

    return any((ac > 0 for ac in allele_count))


def format_variant_record(variant: pysam.VariantRecord, genotypes: str, af_field: str = "AF") -> str:
    """
    Formats a variant record into a string with specified genotypes and variant
    information.

    Args:
        variant (pysam.VariantRecord): The variant record to be formatted.
        genotypes (str): Genotypes information to be included in the format.
        af_field (str): INFO field name to read allele frequency from. Use
            "AF_joint" for gnomAD joint VCFs, "AF" otherwise. Default: "AF".

    Returns:
        str: Formatted variant record as a tab-separated string.
    """

    try:
        af = ",".join(list(map(str, variant.info[af_field])))
    except KeyError:  # catch potential AF missing in variant INFO field
        af = ",".join(["0.0" for _ in variant.alts])  # type: ignore
    variant_format = [
        variant.chrom,
        variant.pos,
        variant.id,
        variant.ref,
        ",".join(variant.alts),  # handle multiallelic sites # type: ignore
        variant.qual,
        ";".join(variant.filter.keys()),  # type: ignore
        f"AF={af}",  # keep allele frequencies
        "GT",  # add genotype to format
        genotypes,
    ]
    return "\t".join([f"{e}" if e is not None else "." for e in variant_format])


def convert_vcf(
    vcf_fname: str,
    samples: List[str],
    joint: bool,
    keep: bool,
    filter_pass_values: frozenset = frozenset({"PASS", "."}),
    af_threshold: float = 0.0,
    output_dir: str = "",
):
    """
    Converts a VCF file by updating the header, filtering variants, and creating
    a new compressed VCF file.

    Args:
        vcf_fname (str): Path to the input VCF file.
        samples (List[str]): List of sample names.
        keep (bool): Flag to determine whether to keep variants that do not pass
            the filter.
        filter_pass_values (frozenset): Set of VCF FILTER values to treat as
            passing. Default: frozenset({"PASS", "."}).
        af_threshold (float): Minimum AF_joint value required to include a variant.
            Default 0.0 (no filtering).
        output_dir (str): Directory for the output VCF. If empty, the output is
            written alongside the input file. Default: "".

    Returns:
        str: Path to the converted and compressed VCF file.

    Raises:
        OSError: If an error occurs while loading the input VCF file.
        IOError: If an error occurs during the conversion process.
    """

    # Determine output path: use output_dir if provided, else same dir as input
    out_basename = f"{os.path.splitext(os.path.basename(vcf_fname))[0]}.gz"
    if output_dir:
        vcf_outfname = os.path.join(output_dir, out_basename)
    else:
        vcf_outfname = os.path.join(os.path.dirname(vcf_fname), out_basename)
    try:
        vcf = load_vcf(vcf_fname)  # use pysam for optimized VCF loading
    except OSError as e:
        raise OSError(f"An error occurred while loading {vcf_fname}") from e
    # Pre-flight: warn if no accepted FILTER records found in the first 10,000 rows.
    # If every record has a hard-fail filter the output will be silently empty.
    _sample = list(itertools.islice(vcf.fetch(), 10_000))
    _n_ok = sum(1 for v in _sample if filter_pass_values & set(v.filter.keys()))
    if not keep and _n_ok == 0 and _sample:
        sys.stderr.write(
            f"WARNING: {vcf_fname}: sampled {len(_sample)} records, none have "
            f"FILTER in {filter_pass_values} — all variants will be skipped. "
            f"Re-run with keep=True or adjust --vcf-filter-pass-values.\n"
        )
    vcf.reset()  # rewind before the main loop below
    af_field = "AF_joint" if joint else "AF"  # field name for AF threshold filtering
    samples_ac = [
        f"AC_joint_{sample}" if joint else f"AC_{sample}" for sample in samples
    ]  # recover allele count field for each input sample
    try:
        with gzip.open(vcf_outfname, mode="wt") as outfile:
            # write the upated header to the converted vcf
            outfile.write(update_header(vcf.header.copy(), samples, joint))
            for variant in vcf:
                if not keep and not (filter_pass_values & set(variant.filter.keys())):
                    continue
                # AF threshold filter (skipped when af_threshold == 0.0)
                if af_threshold > 0.0:
                    try:
                        af_vals = variant.info[af_field]
                        af_max = max(af_vals) if hasattr(af_vals, "__len__") else af_vals
                    except KeyError:
                        af_max = 0.0
                    if af_max < af_threshold:
                        continue
                genotypes = "\t".join(
                    [
                        GT[1] if variant_observed(variant.info[sac]) else GT[0]
                        for sac in samples_ac
                    ]
                )
                outfile.write(
                    f"{format_variant_record(variant, genotypes, af_field)}\n"
                )  # write the formatted VCF line to the out vcf
    except IOError as e:
        raise IOError(f"An error occurred while converting {vcf_fname}") from e
    assert os.stat(vcf_outfname).st_size > 0
    return vcf_outfname


def bcftools_merge(vcf_fname: str, multiallelic: bool) -> str:
    """
    Merges alleles in a VCF file using bcftools and returns the path to the merged
    VCF file.

    Args:
        vcf_fname (str clean): Path to the input VCF file.
        multiallelic (bool): Flag indicating whether to handle multiallelic
            variants, determining the output file suffix.

    Returns:
        str: Path to the merged VCF file.

    Raises:
        subprocess.SubprocessError: If an error occurs during the allele merging
            process.
    """

    # recover file prefix
    vcf_outfname_prefix = os.path.splitext(os.path.splitext(vcf_fname)[0])[0]
    suffix = "multiallelic" if multiallelic else "biallelic"
    vcf_outfname = f"{vcf_outfname_prefix}.{suffix}.vcf.gz"
    m = "-m+" if multiallelic else "-m-"
    code = subprocess.call(
        f"{BCFTOOLSNORM} {m} -O z -o {vcf_outfname} {vcf_fname}",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    if code != 0:
        raise subprocess.SubprocessError(
            f"An error occurred while merging alleles in {vcf_fname}"
        )
    assert os.stat(vcf_outfname).st_size > 0
    remove(vcf_fname)  # remove intermediate vcf file
    return vcf_outfname


def run_conversion_pipeline(
    vcf_fname: str,
    samples: List[str],
    joint: bool,
    keep: bool,
    multiallelic: bool,
    filter_pass_values: frozenset = frozenset({"PASS", "."}),
    af_threshold: float = 0.0,
    output_dir: str = "",
) -> None:
    """
    Runs a conversion pipeline to process a VCF file by adding genotypes and merging
    variants based on specified options.

    Args:
        vcf_fname (str): Path to the input VCF file.
        samples (List[str]): List of sample names.
        keep (bool): Flag to determine whether to keep variants that do not pass
            the filter.
        multiallelic (bool): Flag indicating whether to handle multiallelic variants
            during merging.
        filter_pass_values (frozenset): Set of VCF FILTER values to treat as passing.
        af_threshold (float): Minimum AF_joint value required to include a variant.
            Default 0.0 (no filtering).
        output_dir (str): Directory for output VCFs. Default "" (same as input).

    Returns:
        None
    """

    # add genotypes to input VCF
    vcf_genotypes = convert_vcf(vcf_fname, samples, joint, keep, filter_pass_values, af_threshold, output_dir)
    # merge variants into mutlialleic/biallelic sites
    vcf_merged = bcftools_merge(vcf_genotypes, multiallelic)
    assert os.path.isfile(vcf_merged) and os.stat(vcf_merged).st_size > 0


def convert_gnomad_vcfs() -> None:
    """
    Converts gnomAD VCF files based on specified parameters and sample data.

    Returns:
        None

    Raises:
        ValueError: If the number of input arguments is incorrect, the gnomAD VCF
            directory is invalid, or an invalid number of threads is selected.
        FileNotFoundError: If no gnomAD VCF files are found in the specified directory.
        OSError: If an error occurs during the gnomAD VCF conversion process.
    """

    # read input arguments
    gnomad_vcfs_dir, samples_ids, joint, keep, multiallelic, threads, filter_pass_values, af_threshold, output_dir = (
        parse_commandline(sys.argv[1:])
    )
    start = time.time()
    # recover gnomAD file within the specified location (compressed with bgz extension)
    gnomad_vcfs = glob(os.path.join(gnomad_vcfs_dir, "*.vcf*bgz"))
    if not gnomad_vcfs:  # no gnomAD vcf found
        raise FileNotFoundError(f"No gnomAD VCF file found in {gnomad_vcfs_dir}")
    samples = read_samples_ids(samples_ids)  # recover samples data from sample file
    threads = multiprocessing.cpu_count() if threads == 0 else threads
    try:
        pool = multiprocessing.Pool(processes=threads)
        partial_run_conversion_pipeline = partial(
            run_conversion_pipeline,
            samples=samples,
            joint=joint,
            keep=keep,
            multiallelic=multiallelic,
            filter_pass_values=filter_pass_values,
            af_threshold=af_threshold,
            output_dir=output_dir,
        )
        pool.map(partial_run_conversion_pipeline, gnomad_vcfs)
        pool.close()
        pool.join()
    except OSError as e:
        raise OSError("An error occurred during gnomAD VCF conversion") from e
    sys.stderr.write(f"Elapsed time {(time.time() - start):.2f}s\n")


if __name__ == "__main__":
    convert_gnomad_vcfs()
