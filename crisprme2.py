""" """


from crisprme_version import __version__

from itertools import product
from typing import NoReturn, Tuple, List, Union
from argparse import (
    ArgumentParser,
    Namespace,
    RawDescriptionHelpFormatter,
    _SubParsersAction,
)
from colorama import Fore
from Bio.Seq import Seq

import multiprocessing
import subprocess
import contextlib
import pysam
import sys
import re
import os

# crisprme program name
CRISPRME = "crisprme"
# crisprme's functionalities
COMPLETESEARCH = "complete-search"
COMPLETETEST = "complete-test"
TARGETSINTEGRATION = "targets-integration"
GNOMADCONVERTER = "gnomad-converter"
GENERATEPERSONALCARD = "generate-personal-card"
WEBINTERFACE = "web-interface"
COMMANDS = [
    COMPLETESEARCH,
    COMPLETETEST,
    TARGETSINTEGRATION,
    GNOMADCONVERTER,
    GENERATEPERSONALCARD,
    WEBINTERFACE,
]
# crisprme's directories
CRISPRMEDIRS = [
    "Annotations",
    "Dictionaries",
    "Genomes",
    "PAMs",
    "Results",
    "samplesIDs",
    "VCFs",
]
# utils
DNA = ["A", "C", "G", "T"]
# IUPAC nucleotide code mapping - defining which bases each ambiguous code represents
IUPAC_CODE_MAP = {
    "A": "ARWMDHV",      # A and codes that include A
    "C": "CYSMBHV",      # C and codes that include C
    "G": "GRSKBDV",      # G and codes that include G
    "T": "TYWKBDH",      # T and codes that include T
    "R": "ARWMDHVSKBG",  # A or G (puRine)
    "Y": "CYSMBHVWKDT",  # C or T (pYrimidine)
    "S": "CYSMBHVKDRG",  # C or G (Strong)
    "W": "ARWMDHVYKBT",  # A or T (Weak)
    "K": "GRSKBDVYWHT",  # G or T (Keto)
    "M": "ARWMDHVYSBC",  # A or C (aMino)
    "B": "CYSMBHVRKDGWT", # C, G, or T (not A)
    "D": "ARWMDHVSKBGYT", # A, G, or T (not C)
    "H": "ARWMDHVYSBCKT", # A, C, or T (not G)
    "V": "ARWMDHVYSBCKG", # A, C, or G (not T)
    "N": "ACGTRYSWKMBDHV", # Any nucleotide
}
# script's path locations
SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))
CONDAPATH = "opt/crisprme/PostProcess"
WEBPATH = "opt/crisprme"  # webserver path location



class CRISPRmeArgumentParser(ArgumentParser):
    def error(self, message: str) -> NoReturn:
        command = f" {sys.argv[1]}" if sys.argv[1] in COMMANDS else ""
        errormsg = (
            f"{Fore.RED}\nERROR: {message}.{Fore.RESET}"
            + f"\n\nRun `{CRISPRME}{command} -h` for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error


def check_crisprme_directories() -> None:
    for directory in CRISPRMEDIRS:
        if not os.path.isdir(os.path.join(os.getcwd(), directory)):
            os.makedirs(os.path.join(os.getcwd(), directory))


def create_complete_search_parser(subparsers: _SubParsersAction) -> _SubParsersAction:
    parser_complete_search = subparsers.add_parser(
        COMPLETESEARCH,
        description="Automated end-to-end search pipeline that processes raw input "
        "data through off-target discovery, scoring, annotation, and post-analysis "
        "of results",
        help="perform a comprehensive off-target search across the reference "
        "genome and optionally variant-aware genomes. Includes CFD and CRISTA "
        "scoring, and automated target selection",
    )
    parser_complete_search.add_argument(
        "--genome",
        type=str,
        dest="genome",
        metavar="GENOME-DIR",
        required=True,
        help="path to the reference genome directory. Must contain one FASTA "
        "file per chromosome (unzipped)",
    )
    parser_complete_search.add_argument(
        "--vcf",
        type=str,
        dest="vcf_config",
        metavar="VCF-CONFIG-FILE",
        required=False,
        default="",
        help="path to the configuration file listing one or more VCF directories, "
        "containing one or more bgzipped VCFs (default: no variants used)",
    )
    group = parser_complete_search.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--guide",
        type=str,
        dest="guide",
        metavar="GUIDE-FILE",
        help="path to a file containing guide RNA sequences (one per line, no "
        "FASTA header). Cannot be used with --sequence or --coordinates",
    )
    group.add_argument(
        "--sequence",
        type=str,
        dest="sequence",
        metavar="FASTA-FILE",
        help="FASTA file containing guide sequences. Cannot be used with --guide or --coordinates",
    )
    group.add_argument(
        "--coordinates",
        type=str,
        dest="coordinates",
        metavar="BED-FILE",
        help="BED file with genomic coordinates for guide regions. Cannot be used with --guide or --sequence",
    )
    parser_complete_search.add_argument(
        "--pam",
        type=str,
        dest="pam",
        metavar="PAM",
        required=True,
        help="path to a file containing the PAMto use",
    )
    parser_complete_search.add_argument(
        "--be-window",
        type=str,
        dest="be_window",
        metavar="BASE-EDITING-WINDOW",
        required=False,
        default="",
        help="comma-separated range specifying the base editing window in the "
        "guide (e.g., 4,8)",
    )
    parser_complete_search.add_argument(
        "--be-base",
        type=str,
        dest="be_base",
        metavar="BASE-EDITING-BASE",
        required=False,
        default="",
        help="comma-separated list of nucleotide bases to check for editing "
        "susceptibility (e.g., A,C)",
    )
    parser_complete_search.add_argument(
        "--annotation",
        type=str,
        dest="annotation",
        metavar="ANNOTATION-BED-FILE",
        nargs="*",
        default=[],
        help="BED files with genome annotations (e.g., regulatory elements, "
        "enhancers). The fourth column must contain the annotation name. The "
        "flag accepts multiple input arguments (default: no annotations)",
    )
    parser_complete_search.add_argument(
        "--annotation-colnames",
        type=str,
        metavar="ANNOTATION-COLNAMES",
        dest="annotation_colnames",
        nargs="*",
        default=[],
        help="list of custom column names to use in the final report. Each name "
        "corresponds to one of the input BED files provided with '--annotation'. "
        "Must match the number and order of files in '--annotation' (default: "
        "annotation columns are named 'annotation_<i>')",
    )
    parser_complete_search.add_argument(
        "--gene-annotation",
        type=str,
        dest="gene_annotation",
        metavar="GENE-ANNOTATION-BED-FILE",
        required=False,
        default="",
        help="BED file with gene annotations (e.g., GENCODE). Used to identify "
        "the nearest gene for each target (default: skip nearest gene search)",
    )
    parser_complete_search.add_argument(
        "--samples-ids",
        type=str,
        dest="samples_ids",
        metavar="SAMPLESIDS-CONFIG-FILE",
        required=False,
        default="",
        help="path to config file listing files that contain sample IDs present "
        "in the VCF folders. Required only if --vcf is specified",
    )
    parser_complete_search.add_argument(
        "--mm",
        type=int,
        dest="mm",
        metavar="NUM-MISMATCHES",
        required=True,
        help="maximum number of mismatches allowed between the guide and off-targets",
    )
    parser_complete_search.add_argument(
        "--bDNA",
        type=int,
        dest="bdna",
        metavar="NUM-BULGE-DNA",
        required=False,
        default=0,
        help="maximum number of DNA bulges allowed in the search (default: 0)",
    )
    parser_complete_search.add_argument(
        "--bRNA",
        type=int,
        dest="brna",
        metavar="NUM-BULGE-RNA",
        required=False,
        default=0,
        help="maximum number of RNA bulges allowed in the search (default: 0)",
    )
    parser_complete_search.add_argument(
        "--merge",
        type=int,
        dest="merge_t",
        metavar="MERGE-WINDOW",
        required=False,
        default=3,
        help="nmber of nucleotides defining the merging window for candidate "
        "off-targets. Off-targets falling within this window are clustered, and "
        "the one with the highest score is selected as the representative (pivot) "
        "of the group (default: 3)",
    )
    parser_complete_search.add_argument(
        "--sorting-criteria",
        type=str,
        dest="sorting_criteria",
        metavar="SORTING-CRITERIA",
        required=False,
        default="mm+bulges,mm",
        help="comma-separated list of criteria to sort off-targets. Options: "
        "'mm', 'bulges', 'mm+bulges' (default: 'mm+bulges,mm')",
    )
    parser_complete_search.add_argument(
        "--sorting-criteria-scoring",
        type=str,
        dest="sorting_criteria_scoring",
        metavar="SORTING-CRITERIA-SCORING",
        required=False,
        default="mm+bulges",
        help="comma-separated list of criteria to sort scored off-targets. "
        "Score (CFD or CRISTA) has highest priority. Options: 'mm', 'bulges', "
        "or 'mm+bulges' (default: 'mm+bulges')",
    )
    parser_complete_search.add_argument(
        "--output",
        type=str,
        dest="outdir_complete_search",
        metavar="OUTPUT-DIR",
        required=True,
        help="name of the results directory (output saved in Results/<OUTPUT-DIR>)",
    )
    parser_complete_search.add_argument(
        "--threads",
        type=int,
        dest="threads_complete_search",
        metavar="THREADS-NUM",
        required=False,
        default=8,
        help="number of threads to use. Set to 0 to use all available cores (default: 8)",
    )
    parser_complete_search.add_argument(
        "--debug",
        dest="debug_complete_search",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return parser_complete_search


def create_complete_test_parser(subparsers: _SubParsersAction) -> _SubParsersAction:
    parser_complete_test = subparsers.add_parser(
        COMPLETETEST,
        description="Run an end-to-end test of CRISPRme's complete-search "
        "functionality using example data",
        help="run the full CRISPRme pipeline on example data. Supports testing "
        "on single chromosomes or entire genomes",
    )
    parser_complete_test.add_argument(
        "--chrom",
        type=str,
        dest="chrom",
        metavar="CHROMS",
        required=False,
        default="all",
        help="test the complete-search functionality on the specified chromosome "
        "(e.g., chr22). By default, the test is conducted on all chromosomes. "
        "(defult: all)",
    )
    parser_complete_test.add_argument(
        "--vcf-dataset",
        type=str,
        dest="vcf_dataset",
        metavar="VCF-DATASET",
        required=False,
        default="1000G",
        help="VCFs dataset to be used for CRISPRme's test. Available options "
        "include 1000 Genomes (1000G) and Human Genome Diversity Project (HGDP). "
        "To use the combined dataset type '1000G+HGDP' (default: 1000G)",
    )
    parser_complete_test.add_argument(
        "--threads",
        type=int,
        dest="threads_complete_test",
        metavar="THREADS-NUM",
        required=False,
        default=4,
        help="number of threads to use during test (default: 4)",
    )
    parser_complete_test.add_argument(
        "--debug",
        dest="debug_complete_test",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return subparsers


def create_targets_integration_parser(
    subparsers: _SubParsersAction,
) -> _SubParsersAction:
    parser_targets_integration = subparsers.add_parser(
        TARGETSINTEGRATION,
        description="Automates the integration of CRISPRme in-silico predictions "
        "with user-provided empirical data to generate a high-confidence target panel",
        help="integrate predicted targets with empirical data",
    )
    parser_targets_integration.add_argument(
        "--targets",
        type=str,
        dest="targets_file",
        metavar="TARGETS-FILE",
        required=True,
        help="CRISPRme off-targets report file to use as input for the "
        "integration process",
    )
    parser_targets_integration.add_argument(
        "--empirical-data",
        type=str,
        dest="emprical_data_file",
        metavar="EMPIRICAL-DATA-FILE",
        required=True,
        help="User-provided file containing empirical measurements (e.g., "
        "GUIDE-seq, CIRCLE-seq, or other wet-lab results) to assess and refine "
        "predicted off-targets",
    )
    parser_targets_integration.add_argument(
        "--output",
        type=str,
        dest="outdir_targets_integration",
        metavar="OUTPUT-DIR",
        required=True,
        help="ame of the output directory (integration results will be saved "
        "in Results/<OUTPUT-DIR>)",
    )
    parser_targets_integration.add_argument(
        "--debug",
        dest="debug_targets_integration",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return parser_targets_integration


def create_gnomad_converter_parser(subparsers: _SubParsersAction) -> _SubParsersAction:
    parser_gnomad_converter = subparsers.add_parser(
        GNOMADCONVERTER,
        description="Convert gnomAD VCF files (version >= 3.1) into CRISPRme-compatible "
        "format.\nThis utility ensures structural and content compatibility by "
        "preprocessing input VCFs and incorporating sample IDs for population-aware "
        "variant representation",
        help="convert gnomAD VCFs into CRISPRme-compatible format",
    )
    parser_gnomad_converter.add_argument(
        "--gnomad-vcf-dir",
        type=str,
        dest="gnomad_vcf_dir",
        metavar="GNOMAD-VCF-DIR",
        required=True,
        help="directory containing gnomAD VCF files (with .vcf.bgz extension). "
        "Subdirectories or symbolic links are not automatically followed",
    )
    parser_gnomad_converter.add_argument(
        "--samplesID",
        type=str,
        dest="samples_id",
        metavar="GNOMAD-SAMPLES-IDS-FILE",
        required=True,
        help="precomputed sample IDs file mapping gnomAD populations. This is "
        "required for generating population-stratified output VCFs.",
    )
    parser_gnomad_converter.add_argument(
        "--joint",
        action="store_true",
        dest="joint",
        help="specify this flag if input VCFs contain joint allele frequencies "
        "(e.g., gnomAD v4.1 joint exomes/genomes release).",
    )
    parser_gnomad_converter.add_argument(
        "--keep",
        action="store_true",
        dest="keep",
        help="retain all variants regardless of FILTER status. By default, only "
        "variants with FILTER=PASS are included.",
    )
    parser_gnomad_converter.add_argument(
        "--multiallelic",
        action="store_true",
        dest="multiallelic",
        help="enable merging of variants at the same genomic position to form "
        "multiallelic records. By default, each variant is treated as biallelic.",
    )
    parser_gnomad_converter.add_argument(
        "--threads",
        type=int,
        dest="threads_gnomad_converter",
        metavar="THREADS-NUM",
        required=False,
        default=8,
        help="number of threads to use (default: 8)",
    )
    parser_gnomad_converter.add_argument(
        "--debug",
        dest="debug_gnomad_converter",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return parser_gnomad_converter


def create_personal_card_parser(subparsers: _SubParsersAction) -> _SubParsersAction:
    parser_personal_card = subparsers.add_parser(
        GENERATEPERSONALCARD,
        description="Generate a personal card containing all private off-targets "
        "for a given sample.\nThis tool extracts private targets from CRISPRme "
        "integrated result files based on the provided guide and sample ID",
        help="generate a personal card for a specific sample by extracting all "
        "private targets",
    )
    parser_personal_card.add_argument(
        "--results-dir",
        type=str,
        dest="results_dir",
        metavar="RESULTS-DIR",
        required=True,
        help="path to the CRISPRme results directory from which to extract targets",
    )
    parser_personal_card.add_argument(
        "--guide",
        type=str,
        dest="guide_seq",
        metavar="GUIDE-SEQUENCE",
        required=True,
        help="guide RNA sequence used to extract matching private off-targets",
    )
    parser_personal_card.add_argument(
        "--sample-id",
        type=str,
        dest="sample_id",
        metavar="SAMPLE-ID",
        required=True,
        help="sample identifier for which the personal card should be generated",
    )
    parser_personal_card.add_argument(
        "--debug",
        dest="debug_personal_card",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return parser_personal_card


def create_web_interface_parser(subparsers: _SubParsersAction) -> _SubParsersAction:
    parser_web_interface = subparsers.add_parser(
        WEBINTERFACE,
        description="Start CRISPRme's local web interface server. This command "
        "launches a local server for interactive use of the CRISPRme web interface "
        "via a browser. No input files are required. Once launched, open your "
        "browser at http://localhost:8080 (or the assigned port)",
        help="start CRISPRme's local web interface for interactive use in the browser",
    )
    parser_web_interface.add_argument(
        "--debug",
        dest="debug_web_interface",
        action="store_true",
        help="run in debug mode (local execution vs conda/Docker)",
    )
    return parser_web_interface


def create_parser():
    parser = CRISPRmeArgumentParser(
        prog=CRISPRME,  # Replace with your actual tool name
        description="CRISPRme: A tool for variant- and haplotype-aware CRISPR "
        "off-targets analysis. Integrates robust functionalities for off-target "
        "detection, variant-aware search, and result analysis with support for "
        "genetic diversity assessment in CRISPR-based therapies.\n\n"
        "Important file format requirements:\n"
        f"- All FASTA files used by {CRISPRME} must be unzipped and separated "
        "by chromosome\n"
        f"- All VCF files used by {CRISPRME} must be gzipped with BGZIP and "
        "separated by chromosome",
        epilog=f"For more information on each command, use: {CRISPRME}.py "
        "<command> --help",
        formatter_class=RawDescriptionHelpFormatter,
    )
    # add version argument
    parser.add_argument("--version", action="version", version=__version__)
    # create subparsers for different functionalities
    subparsers = parser.add_subparsers(
        dest="command",
        title="Available commands",
        metavar="",  # needed for help formatting (avoid <command to be displayed>)
        description=None,
    )
    # complete search subcommand
    parser_complete_search = create_complete_search_parser(subparsers)
    # complete test subcommand
    parser_complete_test = create_complete_test_parser(subparsers)
    # targets integration subcommand
    parser_targets_integration = create_targets_integration_parser(subparsers)
    # gnomAD converter subcommand
    parser_gnomad_converter = create_gnomad_converter_parser(subparsers)
    # generate personal card subcommand
    parser_personal_card = create_personal_card_parser(subparsers)
    # web interface subcommand
    parser_web_interface = create_web_interface_parser(subparsers)
    return parser


def validate_file_exists(fname: str, flag: str, parser: CRISPRmeArgumentParser) -> str:
    if not os.path.isfile(fname):
        parser.error(f"The file specified for {flag} does not exist: {fname}")
    return os.path.abspath(fname)

def validate_directory_exists(dirname: str, flag: str, parser: CRISPRmeArgumentParser) -> str:
    if not os.path.isdir(dirname):
        parser.error(f"The directory specified for {flag} does not exist: {dirname}")
    return os.path.abspath(dirname)

def validate_be_window(be_window: str, parser: CRISPRmeArgumentParser) -> Tuple[int, int]:
    start, stop = 1, 0  # dummy window start & stop values
    if be_window:  # window selected, initialize start & stop
        window_coords = be_window.strip().split(",")  # comma-separated
        if len(window_coords) != 2:
            parser.error(
                "base editing window requires exactly a start and stop position. "
                f"Given {len(window_coords)} positions"
            )
        start, stop = window_coords
        if stop < start:
            parser.error(f"base editing window stop < start ({stop} < {start})")
    return int(start), int(stop)

def validate_be_base(be_base: str, parser: CRISPRmeArgumentParser) -> str:
    if not be_base:
        return "none"
    bents = be_base.strip().split(",")  # retrieve be bases
    if len(bents) != 2:
        parser.error(f"given more/less than 2 base editing bases ({len(bents)})") 
    if len(set(bents)) == 1:
        parser.error(f"base editing bases are identical ({be_base})")
    if any(nt not in DNA for nt in bents):
        parser.error("base editing bases appear to not be part of the DNA alphabet")
    return be_base

def validate_annotation(annotation: List[str], parser: CRISPRmeArgumentParser) -> List[str]:
    if not annotation:  # no input annotation file, assign dummy file
        return ["empty.txt"]
    # validate each annotation file
    return [validate_file_exists(fann, "--annotation", parser) for fann in annotation]  

def validate_annotation_colnames(annotation_colnames: List[str], annotations: List[str], parser: CRISPRmeArgumentParser) -> List[str]:
    if len(annotations) == 1 and annotations[0] == "empty.txt":  # no annotation
        return ["annotation_empty"]
    if not annotation_colnames:  # use default column names
        return [f"annotation_{i}" for i, _ in enumerate(annotations)]
    if len(annotation_colnames) != len(annotations):
        parser.error(f"Mismatching number of annotation column names ({len(annotation_colnames)}) and annotation datasets ({len(annotations)})")
    return annotation_colnames

def validate_gene_annotation(gene_annotation: str, parser: CRISPRmeArgumentParser) -> str:
    if not gene_annotation:  # no input gene annotation file, assign dummy file
        return "empty.txt"
    # validate gene annotation file
    return validate_file_exists(gene_annotation, "--gene-annotation", parser)

def validate_mismatches(mm: int, parser: CRISPRmeArgumentParser) -> int:
    if mm < 0:
        parser.error(f"negative number of mismatches selected ({mm})")
    return mm

def validate_bulges(bnum: int, btype: str, parser: CRISPRmeArgumentParser) -> int:
    if bnum < 0:
        parser.error(f"negative number of {btype} bulges selected ({bnum})")
    return bnum

def validate_merge_threshold(merge_t: int, parser: CRISPRmeArgumentParser) -> int:
    if merge_t < 0:
        parser.error(f"negative merge window selected ({merge_t})")
    return merge_t

def validate_sorting_criteria(sorting_criteria: str, flag: str, parser: CRISPRmeArgumentParser) -> str:
    criteria = sorting_criteria.split(",")
    if len(criteria) > len(set(criteria)):
        parser.error(f"repeated sorting criteria in {flag}")
    if len(criteria) > 3:
        parser.error(f"forbidden or repeated sorting criteria in {flag}")
    if any(c not in ["mm+bulges", "mm", "bulges"] for c in criteria):
        parser.error(f"forbidden sorting criteria selected in {flag}")
    return sorting_criteria

def validate_output_dir(outdir: str, parser: CRISPRmeArgumentParser) -> str:
    # output dir must be in Results folder
    outdir = os.path.join(os.getcwd(), CRISPRMEDIRS[4], outdir)
    if not os.path.isdir(outdir):  # output folder not present, create it
        os.makedirs(outdir)
    elif any(os.scandir(outdir)):  # output folder not empty, prevent run
        parser.error(
            f"Output folder {outdir} not empty! Select another output folder.\n"
            f"If the previous run threw an error please consider deleting {outdir}"
        )
    return validate_directory_exists(outdir, "--output", parser)

def validate_threads_num(threads: int, parser: CRISPRmeArgumentParser) -> int:
    if threads < 0:  # threads number cannot be negative 
        parser.error(f"Forbidden number of threads selected ({threads})")
    return multiprocessing.cpu_count() if threads == 0 else threads


def process_pam_file(pam_file: str) -> Tuple[str, int, bool]:
    try:
        with open(pam_file, mode="r") as infile:  # read pam and pam index
            guide_pam, pamidx = infile.readline().strip().split()
    except IOError as e:
        raise IOError(f"error while reading PAM file {pam_file}") from e
    pamidx_ = abs(int(pamidx)) 
    # pam may occur at the beginning or the end of the guide
    pam = guide_pam[:pamidx_] if int(pamidx) < 0 else guide_pam[-pamidx_:] 
    return pam, len(guide_pam), int(pamidx) < 0


def initialize_genome_index(vcf_config: str, pam: str, bmax: int, genomeref: str) -> Tuple[str, bool]:
    if vcf_config != "_":
        with open(vcf_config, mode="r") as infile:
            genome_indexes = ",".join([f"{pam}_{bmax}_{genomeref}_{os.path.basename(line)}" for line in infile if line.strip()])
        return genome_indexes, True
    return f"{pam}_{bmax}_{genomeref}", False

def write_command_line(argv: List[str], outdir: str) -> None:
    with open(os.path.join(outdir, ".crisprme_command_line"), mode="w") as outfile:
        command = " ".join(argv)
        outfile.write(f"CRISPRme input\n{command}\n")


def write_version(outdir: str) -> None:
    with open(os.path.join(outdir, ".crisprme_version"), mode="w") as outfile:
        outfile.write(f"CRISPRme version: {__version__}\n")

def write_params(outdir: str, genomeref: str, genomeidx: str, pam: str, mm: int, bdna: int, brna: int, annotations: List[str], refcmp: bool) -> None:
    with open(os.path.join(outdir, "Params.txt"), mode="w") as outfile:
        genomeref = genomeref.replace(" ", "_")
        outfile.write(f"Genome_selected\t{genomeref}\n")
        outfile.write(f"Genome_ref\t{genomeref}\n")
        outfile.write(f"Genome_idx\t{genomeidx}\n")
        outfile.write(f"Pam\t{pam}\n")
        outfile.write(f"Max_bulges\t{max(bdna, brna)}\n")
        outfile.write(f"Mismatches\t{mm}\n")
        outfile.write(f"DNA\t{bdna}\n")
        outfile.write(f"RNA\t{brna}\n")
        annotations_ = ",".join(annotations)
        outfile.write(f"Annotation\t{annotations_}\n")
        outfile.write(f"Ref_comp\t{refcmp}\n")

def initialize_search(argv: List[str], genome: str, vcf_config: str, pam: str, mm: int, bdna: int, brna: int, annotations: List[str], outdir: str) -> None:
    genomeref = os.path.basename(genome)  # reference genome identifier
    genomeidx, index = initialize_genome_index(vcf_config, pam, max(bdna, brna), genomeref)
    write_command_line(argv, outdir)
    write_version(outdir)
    write_params(outdir, genomeref, genomeidx if max(bdna, brna) > 0 else "None", pam, mm, bdna, brna, annotations, index)

def read_guide_file(guide_file: str) -> List[str]:
    with open(guide_file, mode="r") as infile:
        guides = {line.strip() for line in infile}
    return list(guides)


def expand_iupac_sequence(sequence: str) -> List[str]:
    sequence = sequence.upper()
    # validate sequence contains only valid IUPAC codes
    invalid_chars = set(sequence) - set(IUPAC_CODE_MAP.keys())
    if invalid_chars:
        raise ValueError(f"Invalid IUPAC codes in sequence '{sequence}': {', '.join(invalid_chars)}")
    # get all possible characters for each position
    nt_options = [IUPAC_CODE_MAP[nt] for nt in sequence]
    return ["".join(combination) for combination in product(*nt_options)]  # generate all combinations

def find_pam_positions(sequence: str, pam_variants: List[str]) -> List[int]:
    positions = []
    for pam in pam_variants:  # use lookahead to find overlapping matches
        matches = re.finditer(f"(?={re.escape(pam)})", sequence)
        positions.extend(match.start() for match in matches)
    return sorted(set(positions))  # remove duplicates and sort

def extract_guide_sequences(sequence: str, pam_positions: List[int], guide_length: int, pam_length: int, pam_at_start: bool, reverse_complement: bool = False) -> List[str]:
    guides = []
    sequence_length = len(sequence)
    for pos in pam_positions:
        if pam_at_start:  # PAM is at 5' end: PAM + guide
            start = pos + pam_length
            stop = start + guide_length
            if stop > sequence_length:  # check bounds
                continue
        else:  # PAM is at 3' end: guide + PAM
            stop = pos
            start = stop - guide_length
            if start < 0:   # check bounds
                continue
        guide_seq = sequence[start:stop]
        if guide_seq and len(guide_seq) == guide_length:
            if reverse_complement:
                guide_seq = str(Seq(guide_seq).reverse_complement())
            guides.append(guide_seq)
    return guides


def get_guides(sequence: str, pam: str, guidelen: int, pam_at_start: bool) -> List[str]:
    guides = set()  # avoid duplicates in guides list
    try:  # search guides on forward strand
        pam_variants = expand_iupac_sequence(pam)  # compute pams variants
        pam_positions = find_pam_positions(sequence, pam_variants)  # find pam occurrences
        # extract guide sequences from forward strand
        guides.update(extract_guide_sequences(sequence, pam_positions, guidelen, len(pam), pam_at_start))
    except ValueError as e:
        raise ValueError(f"Error processing forward PAM '{pam}': {e}") from e
    try:  # search guides on reverse strand
        reverse_pam = str(Seq(pam).reverse_complement())  # get reverse complement of PAM
        reverse_pam_variants = expand_iupac_sequence(reverse_pam)
        # find PAM positions on reverse strand (searching forward sequence for reverse PAM)
        reverse_positions = find_pam_positions(sequence, reverse_pam_variants)
        # extract guide sequences from reverse strand
        # for reverse strand, the logic is inverted:
        # - if original PAM is at beginning, reverse PAM guides come from before the PAM
        # - if original PAM is at end, reverse PAM guides come from after the PAM
        guides.update(extract_guide_sequences(sequence, reverse_positions, guidelen, len(pam), not pam_at_start, reverse_complement=True))
    except ValueError as e:
        raise ValueError(f"Error processing reverse PAM '{pam}': {e}") from e
    return list(guides)

def read_fasta_guide(fasta_guide: str, pam: str, guidelen: int, pam_at_start: bool) -> List[str]:
    with open(fasta_guide, mode="r") as infile:
        sequences = [line.strip() for line in infile if not line.startswith(">")]
    guides = []
    for sequence in sequences:
        guides.extend(get_guides(sequence, pam, guidelen, pam_at_start))
    return guides

def extract_guide_bed(genome: str, chrom: str, start: int, stop: int) -> str:
    chromfa = os.path.join(genome, f"{chrom}.fa")
    if not os.path.isfile(chromfa):  # try FASTA extension
        chromfa = os.path.join(genome, f"{chrom}.fasta")
    if not os.path.isfile(chromfa):  # both extension failed
        raise ValueError(f"guide extraction from bed file failed, unable to find {chromfa}")    
    with pysam.FastaFile(chromfa) as fasta:  # open FASTA and extract sequence
        # Note: pysam uses 0-based coordinates, bedtools uses 1-based
        sequence = fasta.fetch(chrom, start - 1, stop)
    subprocess.call(f"rm {chromfa}.fai")  # remove fasta index
    return sequence

def read_bed_guide(bed_guide: str, genome: str) -> List[str]:
    with open(bed_guide, mode="r") as infile:
        guides = {extract_guide_bed(genome, fields[0], int(fields[1].replace(",", "").replace(".", "")), int(fields[2].replace(",", "").replace(".", ""))) for line in infile for fields in [line.strip().split()[:3]]}
    return list(guides)

def format_guide(guide: str, pamlen: int, pam_at_start: bool) -> str:
    ns = "N" * pamlen  # add as many Ns as the pam length
    return f"{ns}{guide}" if pam_at_start else f"{guide}{ns}"

def initialize_input_guides(guide_file: Union[str, None], fasta_guide: Union[str, None], bed_guide: Union[str, None], genomedir: str, pam: str, guidelen: int, pam_at_start: bool, outdir: str) -> str:
    format_n = False
    if guide_file:  # guides given as guide file
        guides = read_guide_file(guide_file)
        format_n = True  # Ns already inserted at pam positions in guide
    elif fasta_guide:  # guides given as FASTA
        guides = read_fasta_guide(fasta_guide, pam, guidelen, pam_at_start)
    elif bed_guide:  # guides given as BED
        guides = read_bed_guide(bed_guide, genomedir)
    else:  # wtf happened?
        raise ValueError("We shouldn't be here")
    guides_file = os.path.join(outdir, "guides.txt")
    with open(guides_file, mode="w") as outfile:
        for guide in guides:
            guide_f = guide if format_n else format_guide(guide, len(pam), pam_at_start)
            outfile.write(f"{guide_f}\n")
    assert os.path.isfile(guides_file) and os.stat(guides_file).st_size > 0
    return guides_file

def establish_script_path_complete_search(debug: bool) -> str:
    if debug:  # run local installation
        sys.stdout.write("\nWarning: running in development mode\n")
        return os.path.join(os.getcwd(), "PostProcess")
    return os.path.join(SCRIPTPATH[:-3], CONDAPATH)  # run global (mamba)


def complete_search_crisprme(args: Namespace, parser: CRISPRmeArgumentParser) -> None:
    # check if all crisprme directories are available, if not create them
    check_crisprme_directories()
    # check input arguments validity 
    genome = validate_directory_exists(args.genome, "--genome", parser)
    if bool(args.vcf_config) and not bool(args.samples_ids):
        parser.error(
            "missing --samples-ids, samples ids file required to perform "
            "variant-aware off-targets search"
        )
    vcf_config = validate_file_exists(args.vcf_config, "--vcf", parser) if args.vcf_config else "_"
    guide_file = validate_file_exists(args.guide, "--guide", parser) if args.guide else None
    fasta_guide = validate_file_exists(args.sequence, "--sequence", parser) if args.sequence else None
    bed_guide = validate_file_exists(args.coordinates, "--coordinates", parser) if args.coordinates else None
    pam_file = validate_file_exists(args.pam, "--pam", parser)
    if bool(args.be_window) and not bool(args.be_base):
        parser.error(
            "base editing window selected, please indicate the bases to check "
            "for editing susceptibility"
        )
    if not bool(args.be_window) and bool(args.be_base):
        parser.error("base editing bases given, please indicate the window to check")
    be_start, be_stop = validate_be_window(args.be_window, parser)
    be_base = validate_be_base(args.be_base, parser)
    annotations = validate_annotation(args.annotation, parser)
    annotation_colnames = validate_annotation_colnames(args.annotation_colnames, annotations, parser)
    gene_annotation = validate_file_exists(args.gene_annotation, "--gene-annotation", parser) if args.gene_annotation else "empty.txt"
    if bool(args.samples_ids) and not bool(args.vcf_config):
        parser.error("--samples-ids selected, but no input arguments detected for --vcf")
    samples_ids = validate_file_exists(args.samples_ids, "--samples-ids", parser) if args.samples_ids else "empty.txt"
    mm = validate_mismatches(args.mm, parser)
    bdna = validate_bulges(args.bdna, "DNA", parser)
    brna = validate_bulges(args.brna, "RNA", parser)
    merge_t = validate_merge_threshold(args.merge_t, parser)
    sorting_criteria = validate_sorting_criteria(args.sorting_criteria, "--sorting-criteria", parser)
    sorting_criteria_score = validate_sorting_criteria(args.sorting_criteria_scoring, "--sorting-criteria-scoring", parser)
    outdir = validate_output_dir(args.outdir_complete_search, parser)
    threads = validate_threads_num(args.threads_complete_search, parser)
    # retrieve pam sequence pam + guide total length, pam position in guide
    pam, guide_pam_len, pam_at_start = process_pam_file(pam_file)  
    # create info files related to the current run
    initialize_search(sys.argv, genome, vcf_config, pam, mm, bdna, brna, annotations, outdir)
    # initialize guides file (guides.txt in output folder)
    guides_file = initialize_input_guides(guide_file, fasta_guide, bed_guide, genome, pam, (guide_pam_len - len(pam)), pam_at_start, outdir)
    log_verbose = os.path.join(outdir, "log_verbose.txt")
    log_error = os.path.join(outdir, "log_error.txt")
    void_mail = "_"
    sys.stdout.write(
        f"\nLaunching job {os.path.basename(outdir)}. The stdout is redirected in "
        f"{log_verbose} and stderr is redirected in {log_error}\n"
    )
    annotations = ",".join(annotations)
    annotation_colnames = ",".join(annotation_colnames)
    script_path = establish_script_path_complete_search(args.debug_complete_search)
    
    # TODO: remove 
    # annotations = os.path.join(script_path, "vuoto.txt")
    # gene_annotation = os.path.join(script_path, "vuoto.txt")
    
    # start search with set parameters
    with open(log_verbose, mode="w") as logv, open(log_error, mode="w") as loge:
            crisprme_run = (
                f"{os.path.join(script_path, 'submit_job_automated_new_multiple_vcfs.sh')} "
                f"{genome} {vcf_config} {guides_file} {pam_file} {annotations} "
                f"{samples_ids} {max(bdna, brna)} {mm} {bdna} {brna} {merge_t} "
                f"{outdir} {script_path} {threads} {os.getcwd()} {gene_annotation} "
                f"{void_mail} {be_start} {be_stop} {be_base} {sorting_criteria_score} "
                f"{sorting_criteria} {annotation_colnames}"
            )
            code = subprocess.call(crisprme_run, shell=True, stderr=loge, stdout=logv)
            if code != 0:
                raise OSError(f"CRISPRme run failed! See {log_error} for details\n")
    # subprocess.call(f"mv {guides_file} {outdir}/.guides.txt", shell=True)
    # subprocess.call(f"mv {outdir}/Params.txt {outdir}/.Params.txt", shell=True)


def complete_test_crisprme(args: Namespace) -> None:
    pass


def personal_card(args: Namespace, parser: CRISPRmeArgumentParser) -> None:
    # check targets directory existence
    targetsdir = validate_directory_exists(args.results_dir, "--results-dir", parser)
    guideseq = args.guide_seq  # guide sequence 
    sample_id = args.sample_id  # sample of interest
    script_path = establish_script_path_complete_search(args.debug_personal_card)
    crisprme_run = (
        f"{os.path.join(script_path, 'generate_sample_card.py')} {targetsdir} "
        f"{guideseq} {sample_id} {script_path}"
    )
    code = subprocess.call(crisprme_run, shell=True)
    if code != 0:
        parser.error(f"CRISPRme personal card generation failed!\n")


def establish_script_path_web_interface(debug: bool) -> str:
    if debug:  # run local installation
        sys.stdout.write("\nWarning: running in development mode\n")
        return os.getcwd()
    return os.path.join(SCRIPTPATH[:-3], WEBPATH)  # run global (mamba)


def web_interface(args: Namespace, parser: CRISPRmeArgumentParser) -> None:
    web_path = establish_script_path_web_interface(args.debug_web_interface)
    crisprme_run = os.path.join(web_path, "index.py")
    with contextlib.suppress(KeyboardInterrupt):
        code = subprocess.call(crisprme_run, shell=True)
        if code != 0:
            parser.error(f"CRISPRme web-interface failed!\n")


def main():
    parser = create_parser()  # initialize parser
    if len(sys.argv) == 1:  # if no arguments provided show help
        parser.print_help()
        sys.exit(os.EX_USAGE)
    args = parser.parse_args()  # parse input arguments
    # execute the appropriate crisprme function based on the command
    if args.command == COMPLETESEARCH:  # run complete-search command
        complete_search_crisprme(args, parser)
    elif args.command == COMPLETETEST:  # run complete-test command
        complete_test_crisprme(args)
    # elif args.command == TARGETSINTEGRATION:  # run targets-integration command
    #     # target_integration()
    #     pass
    # elif args.command == GNOMADCONVERTER:  # run gnomad-converter command
    #     # gnomAD_converter()
    #     pass
    # elif args.command == GENERATEPERSONALCARD:  # run generate-personal-card command
    #     # personal_card()
    #     pass
    elif args.command == WEBINTERFACE:  # run web-interface command
        web_interface(args, parser)
    else:
        # This shouldn't happen with argparse, but just in case
        parser.print_help()


if __name__ == "__main__":
    main()
