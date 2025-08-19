""" """

from crisprme_version import __version__

from typing import NoReturn
from argparse import (
    ArgumentParser,
    Namespace,
    RawDescriptionHelpFormatter,
    _SubParsersAction,
)
from colorama import Fore

import sys
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


class CRISPRmeArgumentParser(ArgumentParser):
    def error(self, message: str) -> NoReturn:
        command = f" {sys.argv[1]}" if sys.argv[1] in COMMANDS else ""
        errormsg = (
            f"{Fore.RED}\nERROR: {message}.{Fore.RESET}"
            + f"\n\nRun {CRISPRME}{command} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error


def error_argparse(error: str, command: str) -> None:
    errormsg = (
        f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
        + f"\n\nRun {CRISPRME} {command} -h for usage\n\n"
    )
    sys.stderr.write(errormsg)  # write error to stderr
    sys.exit(os.EX_USAGE)  # exit execution -> usage error


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
        default="_",
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
        default="none",
        help="comma-separated list of nucleotide bases to check for editing "
        "susceptibility (e.g., A,C)",
    )
    parser_complete_search.add_argument(
        "--annotation",
        type=str,
        dest="annotation",
        metavar="ANNOTATION-BED-FILE",
        nargs="?",
        default="empty.txt",
        help="BED file with genome annotations (e.g., regulatory elements, "
        "enhancers). The fourth column must contain the annotation name (default: "
        "no annotations)",
    )
    parser_complete_search.add_argument(
        "--gene-annotation",
        type=str,
        dest="gene_annotation",
        metavar="GENE-ANNOTATION-BED-FILE",
        required=False,
        default="empty.txt",
        help="BED file with gene annotations (e.g., GENCODE). Used to identify "
        "the nearest gene for each target (default: skip nearest gene search)",
    )
    parser_complete_search.add_argument(
        "--samplesID",
        type=str,
        dest="samples_id",
        metavar="SAMPLESIDS-CONFIG-FILE",
        required=False,
        default="empty.txt",
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
        dest="debug_gnomad_converter",
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


def complete_search_crisprme(args: Namespace) -> None:
    if len(sys.argv) <= 2:  # if no input arguments throw error
        error = f"No arguments provided for {COMPLETETEST} command"
        error_argparse(error, COMPLETETEST)
    print("parsing args...")


def complete_test_crisprme(args: Namespace) -> None:
    pass


def main():
    """Main function with argparse-based command handling."""
    parser = create_parser()

    # If no arguments provided, check directory and show help
    if len(sys.argv) == 1:
        # directoryCheck()
        parser.print_help()
        return

    # Parse arguments
    args = parser.parse_args()

    # Execute the appropriate function based on the command
    if args.command == COMPLETESEARCH:  # run complete-search command
        complete_search_crisprme(args)
    elif args.command == COMPLETETEST:  # run complete-test command
        complete_test_crisprme(args)
    elif args.command == TARGETSINTEGRATION:  # run targets-integration command
        # target_integration()
        pass
    elif args.command == GNOMADCONVERTER:  # run gnomad-converter command
        # gnomAD_converter()
        pass
    elif args.command == GENERATEPERSONALCARD:  # run generate-personal-card command
        # personal_card()
        pass
    elif args.command == WEBINTERFACE:  # run web-interface command
        # web_interface()
        pass
    else:
        # This shouldn't happen with argparse, but just in case
        parser.print_help()


if __name__ == "__main__":
    main()
