#!/usr/bin/env python

from typing import List, NoReturn, Tuple
from Bio.Seq import Seq

import subprocess
import itertools
import sys
import os
import re


version = "2.1.8"  #  CRISPRme version; TODO: update when required
__version__ = version

script_path = os.path.dirname(os.path.abspath(__file__))
origin_path = os.path.dirname(os.path.abspath(__file__))
# path where this file is located
# origin_path = os.path.dirname(os.path.realpath(__file__))
# conda path
conda_path = "opt/crisprme/PostProcess/"
# path corrected to use with conda
corrected_origin_path = script_path[:-3] + conda_path
corrected_web_path = f"{origin_path[:-3]}/opt/crisprme/"
# corrected_web_path = os.getcwd()

script_path = corrected_origin_path
current_working_directory = f"{os.getcwd()}/"
# script_path = corrected_web_path+"/PostProcess/"

input_args = sys.argv

if "--debug" in input_args:
    print("DEBUG MODE")
    script_path = current_working_directory + "PostProcess/"
    corrected_web_path = current_working_directory

VALID_CHARS = {
    "a",
    "A",
    "t",
    "T",
    "c",
    "C",
    "g",
    "G",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "r",
    "y",
    "s",
    "w",
    "k",
    "m",
    "b",
    "d",
    "h",
    "v",
}

CRISPRMEDIRS = [
    "Genomes", "Results", "Dictionaries", "VCFs", "Annotations", "PAMs", "samplesIDs",
]

def is_folder_empty(folder: str) -> bool:
    return any(os.scandir(folder))


# Input chr1:11,130,540-11,130,751
def extractSequence(name, input_range, genome_selected):
    name = "_".join(name.split())
    current_working_directory = os.getcwd() + "/"
    chrom = input_range.split(":")[0]
    start_position = (
        input_range.split(":")[1]
        .split("-")[0]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )
    end_position = (
        input_range.split(":")[1]
        .split("-")[1]
        .replace(",", "")
        .replace(".", "")
        .replace(" ", "")
    )

    list_chr = [
        f
        for f in os.listdir(current_working_directory + "Genomes/" + genome_selected)
        if os.path.isfile(
            os.path.join(current_working_directory + "Genomes/" + genome_selected, f)
        )
        and not f.endswith(".fai")
    ]
    add_ext = ".fa"
    if ".fasta" in list_chr[0]:
        add_ext = ".fasta"
    with open(current_working_directory + name + ".bed", "w") as b:
        b.write(chrom + "\t" + start_position + "\t" + end_position)

    output_extract = subprocess.check_output(
        [
            "bedtools getfasta -fi "
            + current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + add_ext
            + " -bed "
            + current_working_directory
            + name
            + ".bed"
        ],
        shell=True,
    ).decode("utf-8")
    try:
        os.remove(
            current_working_directory
            + "Genomes/"
            + genome_selected
            + "/"
            + chrom
            + ".fa.fai"
        )
    except:
        pass
    try:
        os.remove(current_working_directory + name + ".bed")
    except:
        pass
    ret_string = output_extract.split("\n")[1].strip()
    return ret_string


def getGuides(extracted_seq, pam, len_guide, pam_begin):
    len_pam = len(pam)
    # dict
    len_guide = int(len_guide)
    pam_dict = {
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
    list_prod = []
    for char in pam:
        list_prod.append(pam_dict[char])

    iupac_pam = []  # NNNNNNN NGG
    for element in itertools.product(*list_prod):
        iupac_pam.append("".join(element))

    rev_pam = str(Seq(pam).reverse_complement())
    list_prod = []
    for char in rev_pam:
        list_prod.append(pam_dict[char])

    # CCN NNNNNNN  -> results found with this pam must be reverse complemented
    iupac_pam_reverse = []
    for element in itertools.product(*list_prod):
        iupac_pam_reverse.append("".join(element))

    extracted_seq = extracted_seq.upper()
    len_sequence = len(extracted_seq)
    guides = []
    for pam in iupac_pam:
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    guides.append(extracted_seq[i + len_pam : i + len_pam + len_guide])
                else:
                    if i < len_guide:
                        continue
                    # guides.append(extracted_seq[i-len_guide:i+len_pam])           # i is position where first char of pam is found, eg the N char in NNNNNN NGG
                    # print('1 for:' , extracted_seq[i-len_guide:i])
                    guides.append(extracted_seq[i - len_guide : i])
    for pam in iupac_pam_reverse:  # Negative strand
        pos = [m.start() for m in re.finditer("(?=" + pam + ")", extracted_seq)]
        if pos:
            for i in pos:
                if pam_begin:
                    if i < len_guide:
                        continue
                    guides.append(
                        str(Seq(extracted_seq[i - len_guide : i]).reverse_complement())
                    )
                else:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    # guides.append(str(Seq(extracted_seq[i:i+len_pam+len_guide]).reverse_complement()))         # i is position where first char of pam is found, eg the first C char in CCN NNNNNN
                    # print('2 for:', str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
                    guides.append(
                        str(
                            Seq(
                                extracted_seq[i + len_pam : i + len_guide + len_pam]
                            ).reverse_complement()
                        )
                    )
    return guides
    # return guides for when adding to app.py


def check_crisprme_dirtree() -> None:
    """Ensures that the working directory contains the required CRISPRme folder 
    structure.
    
    Checks for the existence of each expected directory and creates any that are 
    missing.
    """
    # check that working folder respect crisprme's directory tree structure
    for directory in CRISPRMEDIRS:
        if not os.path.exists(os.path.join(current_working_directory, directory)):
            # expected folder not found, create it
            crisprmedir = os.path.join(current_working_directory, directory)
            os.makedirs(crisprmedir)

def print_help_complete_search() -> None:
    """Prints detailed help information for the complete-search functionality.

    Outputs a description of the pipeline and lists all available command-line 
    options to stderr, then exits the program.
    """
    # functionality description
    sys.stderr.write(
        "The complete-search functionality is an end-to-end automated pipeline "
        "that takes raw input files and performs the full workflow up to "
        "post-analysis. Starting from the user-provided genome, variants, guides, "
        "PAM and annotation files, it identifies potential CRISPR off-targets "
        "incorporating variant and haplotype information, and scores each "
        "candidate guide. The pipeline performs genome-wide searches, integrates "
        "annotation data, and generates comprehensive reports."
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--genome, specify the reference genome folder [REQUIRED]\n"
        "\t--vcf, specify a file listing VCF folders (one per line) [OPTIONAL]\n"
        "\t--guide, specify a file containing guide RNAs [REQUIRED if --sequence "
        "not provided]\n"
        "\t--sequence, specify a file with DNA sequences or BED coordinates to "
        "extract guides [REQUIRED if --guide not provided]\n"
        "\t--pam, specify a file containing the PAM sequence [REQUIRED]\n"
        "\t--be-window, specify the window to search for base editor "
        "susceptibility (e.g., --be-window 4,8) [OPTIONAL]\n"
        "\t--be-base, the base(s) for the chosen base editor (e.g., --be-base "
        "A,C) [OPTIONAL]\n"
        "\t--annotation, specify BED files with genome annotations (e.g., "
        "regulatory elements, enhancers). The fourth column must contain the "
        "annotation name. The input BED files must be compressed using bgzip "
        "[OPTIONAL]\n"
        "\t--personal_annotation, specify BED files with personal genomic "
        "annotations. The fourth column must contain the annotation name. The "
        "input BED files must be compressed using bgzip [OPTIONAL]\n"
        "\t--samplesID, specify a file listing sample files (one per line) "
        "present in samplesIDs folder [OPTIONAL]\n"
        "\t--gene_annotation, specify gene annotation (e.g., GENCODE) to find "
        "nearest gene for each target (must be bgzip-compressed) [OPTIONAL]\n"
        "\t--mm, number of mismatches allowed in the search [REQUIRED]\n"
        "\t--bDNA, number of DNA bulges allowed in the search [OPTIONAL]\n"
        "\t--bRNA, number of RNA bulges allowed in the search [OPTIONAL]\n"
        "\t--merge, window size (nucleotides) to merge candidate off-targets "
        "using the highest scoring as pivot [default: 3]\n"
        "\t--sorting-criteria-scoring, comma-separated list to sort targets by "
        "scoring criteria: 'mm', 'bulges', or 'mm+bulges' [default: 'mm+bulges']\n"
        "\t--sorting-criteria, comma-separated list to sort targets by 'mm', "
        "'bulges', or 'mm+bulges' [default: 'mm+bulges,mm']\n"
        "\t--output, specify the output folder name; results will be saved in "
        "Results/<name> [REQUIRED]\n"
        "\t--thread, set number of threads to use [default: 8]\n")
    sys.exit(1)


def error(msg: str) -> NoReturn:
    """Prints an error message to stderr and exits the program.

    This function is used to display error messages and terminate execution with 
    a non-zero status.
    
    Args:
        msg: The error message to display.
    """
    sys.stderr.write(f"Error: {msg}\n")
    sys.exit(1)

def _check_mandatory_args(args: List[str]) -> None:
    if "--genome" not in args:
        error("--genome is required")
    if "--guide" not in args and "--sequence" not in args:
        error("No input guide. One between --guide and --sequence must be specified")
    if "--pam" not in args:
        error("--pam is required")
    if "--mm" not in args:
        error("--mm is required")
    if "--output" not in args:
        error("--output is required")

def _check_genome(args: List[str]) -> str:
    """Retrieves and validates the reference genome directory from command-line 
    arguments.

    Ensures the --genome argument is provided and points to an existing directory. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the reference genome directory.

    Raises:
        SystemExit: If the --genome argument is missing or the specified directory 
            does not exist.
    """
    try:  # read genome input genome folder path
        genomedir = os.path.abspath(args[args.index("--genome") + 1])
    except IndexError:  # no argument for --genome
        error("Missing input for --genome. Reference genome folder must be specified")
    if not os.path.isdir(genomedir):
        error("The folder specified for --genome does not exist")
    return genomedir

def _check_vcf(args: List[str], variant: bool) -> str:
    """Retrieves and validates the VCF configuration file from command-line 
    arguments.

    Ensures the --vcf argument is provided and points to an existing file if 
    variant-aware search is enabled. Raises an error if the argument is missing 
    or invalid.

    Args:
        args: List of command-line arguments.
        variant: Boolean indicating if variant-aware search is enabled.

    Returns:
        The absolute path to the VCF configuration file, or "_" if variant is False.

    Raises:
        SystemExit: If the --vcf argument is missing or the specified file does 
            not exist.
    """
    if not variant:
        return "_"
    try:
        vcfdir = os.path.realpath(args[args.index("--vcf") + 1])
    except IndexError:
        error("Missing input for --vcf. VCF config file must be specified")
    if not os.path.isfile(vcfdir):
        error("The config file specified for --vcf does not exist")
    return vcfdir

def _check_guide(args: List[str], guide: bool) -> str:
    """Retrieves and validates the guide file from command-line arguments.

    Ensures the --guide argument is provided and points to an existing file, and 
    checks for conflicting input flags. Raises an error if the argument is missing,
    invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the guide file.

    Raises:
        SystemExit: If the --guide argument is missing, the specified file does 
            not exist, or there is a conflict with --sequence.
    """
    if not guide:
        return ""
    if "--guide" in args and "--sequence" in args:
        error("Error: Conflicting flags --guide and --sequence. Use only one")
    try:
        guidefile = os.path.abspath(args[args.index("--guide") + 1])
    except IndexError:
        error("Missing input for --guide. Guide file file must be specified")
    if not os.path.isfile(guidefile):
        error("The file specified for --guide does not exist")
    return guidefile

def _check_sequence(args: List[str], sequence: bool) -> str:
    """Retrieves and validates the sequence file from command-line arguments.

    Ensures the --sequence argument is provided and points to an existing file, 
    and checks for conflicting input flags. Raises an error if the argument is 
    missing, invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the sequence file.

    Raises:
        SystemExit: If the --sequence argument is missing, the specified file does 
            not exist, or there is a conflict with --guide.
    """
    if not sequence:
        return ""
    if "--guide" in args and "--sequence" in args:
        error("Error: Conflicting flags --guide and --sequence. Use only one")
    try:
        sequence_file = os.path.abspath(args[args.index("--sequence") + 1])
    except IndexError:
        error("Missing input for --sequence. Guide file file must be specified")
    if not os.path.isfile(sequence_file):
        error("The file specified for --sequence does not exist")
    return sequence_file

def _check_pam(args: List[str]) -> str:
    """Retrieves and validates the PAM file from command-line arguments.

    Ensures the --pam argument is provided and points to an existing file. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the PAM file.

    Raises:
        SystemExit: If the --pam argument is missing or the specified file does 
            not exist.
    """
    try:
        pamfile = os.path.abspath(args[args.index("--pam") + 1])
    except IndexError:
        error("Missing input for --pam. PAM file file must be specified")
        exit(1)
    if not os.path.isfile(pamfile):
        error("The file specified for --pam does not exist")
    return pamfile

def _check_be_window(args: List[str], be_window: bool) -> Tuple[int, int]:
    """Retrieves and validates the base editing window from command-line arguments.

    Ensures the --be-window argument is provided and correctly formatted, and 
    checks for required dependencies. Raises an error if the argument is missing, 
    invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        A tuple containing the start and end positions of the base editing window.

    Raises:
        SystemExit: If the --be-window argument is missing, incorrectly formatted, 
            or if --be-base is not provided when required.
    """
    if not be_window:
        return 1, 0
    if "--be-window" in args and "--be-base" not in args:
        error(
            "Missing --be-base argument. Please input the base editor to check "
            "in the specified window"
        )
    try:
        base_window = args[args.index("--be-window") + 1]
    except IndexError:
        error("Missing input for --be-window. Base editing window file must be specified")
    try:
        base_start, base_end = base_window.strip().split(",")
        base_start, base_end = int(base_start), int(base_end)
    except Exception:
        error("Invalid base editing window specified")
    if base_end < base_start:
        error("Invalid base editing window specified")
    return base_start, base_end

def _check_be_base(args: List[str], be_base: bool) -> str:
    """Retrieves and validates the base editor set from command-line arguments.

    Ensures the --be-base argument is provided and contains only valid nucleotide 
    characters, and checks for required dependencies. Raises an error if the argument 
    is missing, invalid, or in conflict.

    Args:
        args: List of command-line arguments.

    Returns:
        The base editor set as a string.

    Raises:
        SystemExit: If the --be-base argument is missing, contains invalid characters, 
            or --be-window is not provided when required.
    """
    if not be_base:
        return "none"
    if "--be-base" in args and "--be-window" not in args:
        error(
            "Missing --be-window argument. Please input the base editing window "
            "to check for the specified editor"
        )
    try:
        base_set = args[args.index("--be-base") + 1]
    except IndexError:
        error("Missing input for --be-base. Base editor file must be specified")
    if any(nt not in VALID_CHARS for nt in base_set.strip().split(",")):
        error("Invalid editor specified")
    return base_set 

def _decompress_file(fname: str, outfname: str) -> str:
    """Decompresses a bgzipped file to a specified output file.

    Uses gunzip to decompress the input file and writes the result to the output 
    file. Raises an error if decompression fails.

    Args:
        fname: Path to the gzipped input file.
        outfname: Path where the decompressed file will be written.

    Returns:
        The path to the decompressed output file.

    Raises:
        SystemExit: If decompression fails.
    """
    code = subprocess.call(f"gunzip -k -c {fname} > {outfname}", shell=True)
    if code != 0:
        error("Decompressing file failed")
    assert os.path.isfile(outfname)
    return outfname

def _compress_file(fname: str) -> str:
    """Compresses a file using bgzip and returns the path to the compressed file.

    Uses bgzip to compress the specified file. Raises an error if compression fails.

    Args:
        fname: Path to the file to be compressed.

    Returns:
        The path to the compressed file with a .gz extension.

    Raises:
        SystemExit: If compression fails.
    """
    code = subprocess.call(f"bgzip -f {fname}", shell=True)
    if code != 0:
        error("Compressing and indexing file failed")
    assert os.path.isfile(f"{fname}.gz")
    return f"{fname}.gz"

def _sort_bed(fname: str, outfname: str) -> str:
    """Sorts a BED file and writes the sorted output to a new file.

    Uses sort-bed to sort the input BED file and saves the result to the specified 
    output file. Raises an error if sorting fails.

    Args:
        fname: Path to the input BED file.
        outfname: Path where the sorted BED file will be written.

    Returns:
        The path to the sorted BED file.

    Raises:
        SystemExit: If sorting fails.
    """
    code = subprocess.call(f"sort-bed {fname} > {outfname}", shell=True)
    if code != 0:
        error("Sorting BED file failed")
    assert os.path.isfile(outfname)
    return outfname

def _cat_files(fname1: str, fname2: str, outfname: str) -> str:
    """Concatenates two files and writes the result to a new file.

    Uses the cat command to combine the contents of two files into a single output 
    file. Raises an error if concatenation fails.

    Args:
        fname1: Path to the first input file.
        fname2: Path to the second input file.
        outfname: Path where the concatenated file will be written.

    Returns:
        The path to the concatenated output file.

    Raises:
        SystemExit: If concatenation fails.
    """
    code = subprocess.call(f"cat {fname1} {fname2} > {outfname}", shell=True)
    if code != 0:
        error("Concatenating files failed")
    assert os.path.isfile(outfname)
    return outfname

def _mv_file(fname: str, outfname: str) -> str:
    """Renames or moves a file to a new location.

    Uses the mv command to move or rename the specified file. Raises an error if 
    the operation fails.

    Args:
        fname: Path to the source file.
        outfname: Path to the destination file.

    Returns:
        The path to the moved or renamed file.

    Raises:
        SystemExit: If the move or rename operation fails.
    """
    code = subprocess.call(f"mv {fname} {outfname}", shell=True)
    if code != 0:
        error("Renaming file failed")
    assert os.path.isfile(outfname)
    return outfname

def _rm_files(fnames: List[str]) -> None:
    """Removes a list of files from the filesystem.

    Iterates over the provided list of file paths and deletes each file. Raises 
    an error if any file cannot be removed.

    Args:
        fnames: List of file paths to remove.

    Raises:
        SystemExit: If removing any file fails.
    """
    for fname in fnames:
        code = subprocess.call(f"rm {fname}", shell=True)
        if code != 0:
            error("Failed removing file")


def _sort_annotation(annotationfile: str) -> str:
    """Sorts, compresses, and replaces a BED annotation file for downstream 
    analysis.

    Decompresses the input annotation file, sorts it, compresses it with bgzip, 
    and replaces the original file. Raises an error if any step fails.

    Args:
        annotationfile: Path to the input annotation file.

    Returns:
        The path to the sorted and compressed annotation file.

    Raises:
        SystemExit: If decompression, sorting, compression, or renaming fails.
    """
    annotationfile_decompressed = _decompress_file(annotationfile, f"{annotationfile}.tmp.bed")
    annotationfile_sorted = _sort_bed(annotationfile_decompressed, f"{annotationfile}.tmp.sorted.bed")
    annotationfile_sorted_bgzip = _compress_file(annotationfile_sorted)
    annfile = _mv_file(annotationfile_sorted_bgzip, annotationfile)
    _rm_files([annotationfile_decompressed])  # remove tmp files
    return annfile


def _check_annotation(args: List[str], annotation: bool) -> str:
    """Retrieves and validates the annotation file from command-line arguments.

    Ensures the --annotation argument is provided and points to an existing file, 
    or returns a mock file if not specified. Raises an error if the argument is 
    missing or invalid.

    Args:
        args: List of command-line arguments.
        annotation: Boolean indicating if annotation is required.

    Returns:
        The absolute path to the annotation file, or a mock file if annotation is 
            not required.

    Raises:
        SystemExit: If the --annotation argument is required but missing, or the 
            specified file does not exist.
    """
    if not annotation:
        return os.path.join(script_path, "vuoto.txt")  # mock annotation file
    try:
        annotationfile = os.path.abspath(args[args.index("--annotation") + 1])
    except IndexError:
        error("Missing input for --annotation. Annotation file must be specified")
    if not os.path.isfile(annotationfile):
        error("The file specified for --annotation does not exist")
    return _sort_annotation(annotationfile)  # sort input annotation file


def _tag_personal_annotation(fname: str, outfname) -> str:
    """Tags the fourth column of a BED file as personal and writes the result to 
    a new file.

    Modifies the annotation name in the fourth column by appending '_personal', 
    replaces spaces with tabs, and updates commas in the annotation name. Raises 
    an error if tagging fails.

    Args:
        fname: Path to the input BED file.
        outfname: Path where the tagged BED file will be written.

    Returns:
        The path to the tagged BED file.

    Raises:
        SystemExit: If tagging fails.
    """
    code = subprocess.call(
        f"awk '$4 = $4\"_personal\"' {fname} | sed \"s/ /\t/g\" | sed "
        f"\"s/,/_personal,/g\" > {outfname}", 
        shell=True,
    )
    if code != 0:
        error("Tagging personal annotation file failed")
    assert os.path.isfile(outfname)
    return outfname

def _process_personal_annotation(personal_annotationfile: str, annotationfile: str) -> str:
    """Integrates personal and reference annotation files into a single sorted 
    and compressed BED file.

    Decompresses and tags the personal annotation file, merges it with the reference 
    annotation file if present, sorts the combined file, compresses it, and returns 
    the path to the final file. Raises an error if any step fails.

    Args:
        personal_annotationfile: Path to the personal annotation BED file.
        annotationfile: Path to the reference annotation BED file.

    Returns:
        The path to the sorted and compressed combined annotation file.

    Raises:
        SystemExit: If decompression, tagging, concatenation, sorting, or 
            compression fails.
    """
    pannotation_decompressed = _decompress_file(
        personal_annotationfile, f"{personal_annotationfile}.tmp.bed"
    )
    pannotation_tag = f"{personal_annotationfile}.tmp.tag.bed"
    pannotation_tag = _tag_personal_annotation(pannotation_decompressed, pannotation_tag)
    _rm_files([pannotation_decompressed])  # remove tmp files
    concat_annotationfile = os.path.join(
        os.path.abspath(os.path.dirname(personal_annotationfile)),
        "annotation+personal.bed"
    )
    if annotationfile == os.path.join(script_path, "vuoto.txt"):
        concat_annotationfile = _mv_file(pannotation_tag, concat_annotationfile)
    else:  # concatenate personal and annotation file
        annotation_decompressed = _decompress_file(
            annotationfile, f"{annotationfile}.tmp.bed"
        )
        concat_annotationfile = _cat_files(
            annotation_decompressed, pannotation_tag, concat_annotationfile
        )
        _rm_files([annotation_decompressed, pannotation_tag])
    # sort concatenated annotation files
    concat_annotationfile_sorted = f"{concat_annotationfile}.sorted.bed"
    concat_annotationfile_sorted = _sort_bed(concat_annotationfile, concat_annotationfile_sorted)
    concat_annotationfile_bgzip = _compress_file(concat_annotationfile_sorted)
    _rm_files([concat_annotationfile, f"{concat_annotationfile_bgzip}.tbi"])
    assert os.path.isfile(concat_annotationfile_bgzip)
    return concat_annotationfile_bgzip
    

def _check_personal_annotation(args: List[str], annotationfile: str, personal_annotation: bool) -> str:
    """Retrieves and processes the personal annotation file if specified.

    Checks for the --personal_annotation argument, validates the file, and integrates 
    it with the reference annotation file. Returns the processed annotation file path, 
    or the original annotation file if no personal annotation is provided.

    Args:
        args: List of command-line arguments.
        annotationfile: Path to the reference annotation file.
        personal_annotation: Boolean indicating if personal annotation is required.

    Returns:
        The path to the processed annotation file.

    Raises:
        SystemExit: If the --personal_annotation argument is required but missing, 
            or the specified file does not exist.
    """
    if not personal_annotation:
        return annotationfile
    try:
        personal_annotationfile = os.path.abspath(args[args.index("--personal_annotation") + 1])
    except IndexError:
        error("Missing input for --personal_annotation. Annotation file must be specified")
    if not os.path.isfile(personal_annotationfile):
        error("The file specified for --personal_annotation does not exist")
    return _process_personal_annotation(personal_annotationfile, annotationfile)

def _check_samples_ids(args: List[str], variant: bool) -> str:
    """Retrieves and validates the samples ID file from command-line arguments.

    Ensures the --samplesID argument is provided and points to an existing file 
    if variant-aware search is enabled. Returns a mock file if variant is not used. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        variant: Boolean indicating if variant-aware search is enabled.

    Returns:
        The absolute path to the samples ID file, or a mock file if variant is False.

    Raises:
        SystemExit: If the --samplesID argument is missing or the specified file 
            does not exist.
    """
    if variant and "--samplesID" not in args:
        error("Missing --samplesID argument for variant-aware offtargets search")
    if not variant and "--samplesID" in args:
        error("Missing --samplesID selected, but missing --vcf argument")
    if not variant:  # use mock file for samples if variant not used
        return os.path.join(script_path, "vuoto.txt")
    try:
        samplefile = os.path.abspath(args[args.index("--samplesID") + 1])
    except IndexError:
        error("Missing input for --samplesID. Samples file must be specified")
        exit(1)
    if not os.path.isfile(samplefile):
        error("The file specified for --samplesID does not exist")
    return samplefile

def _check_gene_annotation(args: List[str], geneann: bool) -> str:
    """Retrieves and validates the gene annotation file from command-line arguments.

    Ensures the --gene_annotation argument is provided and points to an existing file, 
    or returns a mock file if not specified. Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        geneann: Boolean indicating if gene annotation is required.

    Returns:
        The absolute path to the gene annotation file, or a mock file if gene 
            annotation is not required.

    Raises:
        SystemExit: If the --gene_annotation argument is required but missing, 
            or the specified file does not exist.
    """
    if not geneann:
        return os.path.join(script_path, "vuoto.txt")
    try:
        gene_annotation = os.path.abspath(args[args.index("--gene_annotation") + 1])
    except IndexError:
        error("Missing input for --gene_annotation. Gene annotation file must be specified")
    if not os.path.isfile(gene_annotation):
        error("The file specified for --gene_annotation does not exist")
    return gene_annotation     

def _check_mm(args: List[str]) -> int:
    """Retrieves and validates the number of mismatches from command-line arguments.

    Ensures the --mm argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.

    Returns:
        The number of mismatches as an integer.

    Raises:
        SystemExit: If the --mm argument is missing or the value is negative.
    """
    try:
        mm = int(args[args.index("--mm") + 1])
    except IndexError:
        error("Missing input for --mm. Mismatches number must be specified")
    if mm < 0:
        error("Invalid number of mismatches specified")
    return mm

def _check_bdna(args: List[str], bdna: bool) -> int:
    """Retrieves and validates the number of DNA bulges from command-line arguments.

    Ensures the --bDNA argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        bdna: Boolean indicating if DNA bulges are required.

    Returns:
        The number of DNA bulges as an integer.

    Raises:
        SystemExit: If the --bDNA argument is missing or the value is negative.
    """
    if not bdna:
        return 0
    try:
        bDNA = int(args[args.index("--bDNA") + 1])
    except IndexError:
        error("Missing input for --bDNA. DNA bulges number must be specified")
    if bDNA < 0:
        error("Invalid number of DNA bulges specified")
    return bDNA

def _check_brna(args: List[str], brna: bool) -> int:
    """Retrieves and validates the number of RNA bulges from command-line arguments.

    Ensures the --bRNA argument is provided and is a non-negative integer. Raises 
    an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        brna: Boolean indicating if RNA bulges are required.

    Returns:
        The number of RNA bulges as an integer.

    Raises:
        SystemExit: If the --bRNA argument is missing or the value is negative.
    """
    if not brna:
        return 0
    try:
        bRNA = int(args[args.index("--bRNA") + 1])
    except IndexError:
        error("Missing input for --bRNA. RNA bulges number must be specified")
    if bRNA < 0:
        error("Invalid number of RNA bulges specified")
    return bRNA

def _check_merge(args: List[str], merge: bool) -> int:
    """Retrieves and validates the merge threshold from command-line arguments.

    Ensures the --merge argument is provided and is a non-negative integer. 
    Returns a default value if not specified. Raises an error if the argument is 
    missing or invalid.

    Args:
        args: List of command-line arguments.
        merge: Boolean indicating if merge threshold is required.

    Returns:
        The merge threshold as an integer.

    Raises:
        SystemExit: If the --merge argument is missing or the value is negative.
    """
    if not merge:
        return 3
    try:
        merge_t = int(args[args.index("--merge") + 1])
    except IndexError:
        error("Missing input for --merge. Merge threshold must be specified")
    if merge_t < 0:
        error("Invalid merge threshold specified")
    return merge_t

def _check_sorting_criteria_scoring(args: List[str], sorting_criteria: bool) -> str:
    """Retrieves and validates the sorting criteria for scoring from command-line 
    arguments.

    Ensures the --sorting-criteria-scoring argument is provided, contains valid 
    and non-repeated criteria, and does not exceed the allowed number of criteria. 
    Raises an error if the argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        sorting_criteria: Boolean indicating if sorting criteria for scoring is 
            required.

    Returns:
        The sorting criteria for scoring as a comma-separated string.

    Raises:
        SystemExit: If the --sorting-criteria-scoring argument is missing, contains 
            forbidden or repeated criteria, or exceeds the allowed number of criteria.
    """
    if not sorting_criteria:
        return "mm+bulges"
    try:
        sorting_criteria_scoring = args[args.index("--sorting-criteria-scoring") + 1]
    except IndexError:
        error(
            "Missing input for --sorting-criteria-scoring. Sorting criteria "
            "(scoring) must be specified"
        )
    if len(sorting_criteria_scoring.split(",")) > len(
        set(sorting_criteria_scoring.split(","))
    ):
        error("Repeated sorting criteria (scoring)\n")
    if len(sorting_criteria_scoring.split(",")) > 3:
        error("Forbidden or repeated sorting criteria (scoring)\n")
    if any(
        c not in ["mm+bulges", "mm", "bulges"]
        for c in sorting_criteria_scoring.split(",")
    ):
        error("Forbidden sorting criteria (scoring) selected\n")
    return sorting_criteria_scoring

def _check_sorting_criteria(args: List[str], sorting_criteria: bool) -> str:
    """Retrieves and validates the sorting criteria from command-line arguments.

    Ensures the --sorting-criteria argument is provided, contains valid and non-repeated 
    criteria, and does not exceed the allowed number of criteria. Raises an error if the 
    argument is missing or invalid.

    Args:
        args: List of command-line arguments.
        sorting_criteria: Boolean indicating if sorting criteria is required.

    Returns:
        The sorting criteria as a comma-separated string.

    Raises:
        SystemExit: If the --sorting-criteria argument is missing, contains forbidden 
            or repeated criteria, or exceeds the allowed number of criteria.
    """
    if not sorting_criteria:
        return "mm+bulges"
    try:
        sorting_criteria_fewest = args[args.index("--sorting-criteria") + 1]
    except IndexError:
        error(
            "Missing input for --sorting-criteria. Sorting criteria "
            "must be specified"
        )
    if len(sorting_criteria_fewest.split(",")) > len(
        set(sorting_criteria_fewest.split(","))
    ):
        error("Repeated sorting criteria\n")
    if len(sorting_criteria_fewest.split(",")) > 3:
        error("Forbidden or repeated sorting criteria\n")
    if any(
        c not in ["mm+bulges", "mm", "bulges"]
        for c in sorting_criteria_fewest.split(",")
    ):
        error("Forbidden sorting criteria selected\n")
    return sorting_criteria_fewest

def _check_output(args: List[str]) -> str:
    """Retrieves and validates the output folder from command-line arguments.

    Ensures the --output argument is provided and points to a valid directory. 
    Raises an error if the argument is missing, the folder does not exist, or is 
    not empty.

    Args:
        args: List of command-line arguments.

    Returns:
        The absolute path to the output folder.

    Raises:
        SystemExit: If the --output argument is missing, the folder does not exist, 
            or is not empty.
    """
    try:
        outputfolder = os.path.join(
            current_working_directory, CRISPRMEDIRS[1], args[args.index("--output") + 1]
        )
    except IndexError:
        error("Missing input for --output. Output folder must be specified")
    if os.path.isdir(outputfolder):  # check whether the folder is present or not
        if is_folder_empty(outputfolder):  # if present check if not empty
            error(
                f"Output folder {outputfolder} not empty!Select another "
                "output folder for the current CRISPRme run.If the previous "
                "run using the following folder threw an error, please delete "
                f"{outputfolder} before running a new CRISPRme search."
            )
    else:  # old folder doesn't exist, create it
        os.makedirs(outputfolder)    
    if not os.path.isdir(outputfolder):
        error("The folder specified for --output does not exist")
    return outputfolder

def _check_threads(args: List[str], threads: bool) -> int:
    """Retrieves and validates the number of threads from command-line arguments.

    Ensures the --thread argument is provided and is a positive integer. Returns 
    a default value if not specified. Raises an error if the argument is missing 
    or invalid.

    Args:
        args: List of command-line arguments.
        threads: Boolean indicating if the number of threads is specified.

    Returns:
        The number of threads as an integer.

    Raises:
        SystemExit: If the --thread argument is missing or the value is less than 1.
    """
    if not threads:
        return 8  # default use 8 threads
    try:
        thread = int(args[args.index("--thread") + 1])
    except IndexError:
        error("Missing input for --thread. Number of threads must be specified")
    if thread < 1:
        error("Invalid number of threads specified")
    return thread

def complete_search() -> None:
    args = input_args[2:]  # retrieve complete-search input arguments
    if "--help" in args or not args:  # print help
        print_help_complete_search()
    check_crisprme_dirtree()  # check crisprme directory tree structure
    _check_mandatory_args(args)  # check mandatory arguments
    genomedir = _check_genome(args)  # input genome folder
    variant = "--vcf" in args  # variant-aware search?
    vcfdir = _check_vcf(args, variant)  # input variants dataset
    guidefile = _check_guide(args, "--guide" in args)  # guide file
    sequence_file = _check_sequence(args, "--sequence" in args)  # sequence
    sequence_use = bool(sequence_file)  # use sequence file for guide
    assert sum([bool(guidefile), bool(sequence_file)]) == 1
    pamfile = _check_pam(args)  # pam file
    base_start, base_end = _check_be_window(args, "--be-window" in args)  # base editing window
    base_set = _check_be_base(args, "--be-base" in args)  # base editing bases
    annotationfile = _check_annotation(args, "--annotation" in args)  # annotation file
    annotationfile = _check_personal_annotation(args, annotationfile, "--personal_annotation" in args) # personal annotation file
    samplefile = _check_samples_ids(args, variant)  # samples ids file
    gene_annotation = _check_gene_annotation(args, "--gene_annotation" in args)  # gene annotation file
    mm = _check_mm(args)  # mismatches
    bDNA, bRNA = _check_bdna(args, "--bDNA" in args), _check_brna(args, "--bRNA" in args)  # bulges
    bMax = max(bDNA, bRNA)  # maximum number of bulges
    merge_t = _check_merge(args, "--merge" in args)  # merge threshold
    sorting_criteria_scoring = _check_sorting_criteria_scoring(args, "--sorting-criteria-scoring" in args)  # sorting criteria score columns
    sorting_criteria = _check_sorting_criteria(args, "--sorting-criteria" in args)  # sorting criteria (mm+bulges) columns
    outputfolder = _check_output(args)  # output folder
    thread = _check_threads(args, "--thread" in args)  # number of threads

    # extract pam seq from file
    pam_len = 0
    total_pam_len = 0
    with open(pamfile, "r") as pam_file:
        pam_char = pam_file.readline()
        total_pam_len = len(pam_char.split(" ")[0])
        index_pam_value = pam_char.split(" ")[-1]
        if int(pam_char.split(" ")[-1]) < 0:
            end_idx = int(pam_char.split(" ")[-1]) * (-1)
            pam_char = pam_char.split(" ")[0][0:end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(" ")[-1])
            pam_char = pam_char.split(" ")[0][end_idx * (-1) :]
            pam_len = end_idx
            pam_begin = False

    genome_ref = os.path.basename(genomedir)
    annotation_name = os.path.basename(annotationfile)
    nuclease = os.path.basename(pamfile).split(".")[0].split("-")[2]
    if bMax != 0:
        search_index = True
    else:
        search_index = False
    if variant:
        genome_idx_list = []
        with open(vcfdir, "r") as vcfs:
            for line in vcfs:
                if line.strip():
                    if line[-2] == "/":
                        line = line[:-2]
                    base_vcf = os.path.basename(line)
                    genome_idx_list.append(
                        pam_char
                        + "_"
                        + str(bMax)
                        + "_"
                        + genome_ref
                        + "+"
                        + base_vcf.strip()
                    )
        genome_idx = ",".join(genome_idx_list)
        ref_comparison = True
    else:
        genome_idx = pam_char + "_" + str(bMax) + "_" + genome_ref
        ref_comparison = False
    # os.chdir(script_path)
    # write crisprme version to file
    with open(outputfolder + "/.command_line.txt", "w") as p:
        p.write("input_command\t" + " ".join(sys.argv[:]))
        p.write("\n")
        p.close()
    with open(outputfolder + "/.version.txt", "w") as p:
        p.write("crisprme_version\t" + __version__)
        p.write("\n")
        p.close()
    # write parameters to file
    with open(outputfolder + "/Params.txt", "w") as p:
        p.write("Genome_selected\t" + genome_ref.replace(" ", "_") + "\n")
        p.write("Genome_ref\t" + genome_ref + "\n")
        if search_index:
            p.write("Genome_idx\t" + genome_idx + "\n")
        else:
            p.write("Genome_idx\t" + "None\n")
        p.write("Pam\t" + pam_char + "\n")
        p.write("Max_bulges\t" + str(bMax) + "\n")
        p.write("Mismatches\t" + str(mm) + "\n")
        p.write("DNA\t" + str(bDNA) + "\n")
        p.write("RNA\t" + str(bRNA) + "\n")
        p.write("Annotation\t" + str(annotation_name) + "\n")
        p.write("Nuclease\t" + str(nuclease) + "\n")
        # p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write("Ref_comp\t" + str(ref_comparison) + "\n")
        p.close()
    len_guide_sequence = total_pam_len - pam_len
    if sequence_use:
        guides = list()
        text_sequence = str()
        for line in open(sequence_file, "r"):
            text_sequence += line
        for name_and_seq in text_sequence.split(">"):
            if "" == name_and_seq:
                continue
            name = name_and_seq[: name_and_seq.find("\n")]
            seq = name_and_seq[name_and_seq.find("\n") :]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if "chr" in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split("\n"):
                    if "" == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = (
                        pieces_of_row[0]
                        + ":"
                        + pieces_of_row[1]
                        + "-"
                        + pieces_of_row[2]
                    )
                    extracted_seq = extractSequence(
                        name, seq_to_extract, genome_ref.replace(" ", "_")
                    )
                    guides.extend(
                        getGuides(
                            extracted_seq, pam_char, len_guide_sequence, pam_begin
                        )
                    )
            else:
                seq = seq.split()
                seq = "".join(seq)
                extracted_seq = seq.strip()
                guides.extend(
                    getGuides(extracted_seq, pam_char, len_guide_sequence, pam_begin)
                )
        temp_guides = list()
        for guide in guides:
            addN = "N" * pam_len
            if pam_begin:
                temp_guides.append(addN + guide)
            else:
                temp_guides.append(guide + addN)
        if len(temp_guides) > 1000000000:
            temp_guides = temp_guides[:1000000000]
        guides = temp_guides
        extracted_guides_file = open(outputfolder + "/guides.txt", "w")
        for guide in guides:
            extracted_guides_file.write(guide + "\n")
        extracted_guides_file.close()
    # print(guides)
    # exit(0)
    void_mail = "_"
    if sequence_use == False:
        os.system(f"cp {guidefile} {outputfolder}/guides.txt")
    print(
        f"Launching job {outputfolder}. The stdout is redirected in log_verbose.txt and stderr is redirected in log_error.txt"
    )
    # start search with set parameters
    with open(f"{outputfolder}/log_verbose.txt", "w") as log_verbose:
        with open(f"{outputfolder}/log_error.txt", "w") as log_error:
            crisprme_run = (
                f"{os.path.join(script_path, 'submit_job_automated_new_multiple_vcfs.sh')} "
                f"{genomedir} {vcfdir} {os.path.join(outputfolder, 'guides.txt')} "
                f"{pamfile} {annotationfile} {samplefile} {bMax} {mm} {bDNA} {bRNA} "
                f"{merge_t} {outputfolder} {script_path} {thread} {current_working_directory} "
                f"{gene_annotation} {void_mail} {base_start} {base_end} {base_set} "
                f"{sorting_criteria_scoring} {sorting_criteria}"
            )
            code = subprocess.call(
                crisprme_run, shell=True, stderr=log_error, stdout=log_verbose
            )
            if code != 0:
                raise OSError(
                    f"\nCRISPRme run failed! See {os.path.join(outputfolder, 'log_error.txt')} for details\n"
                )
            # subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), str(vcfdir), str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(
            #     samplefile), str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # else:
    #     with open(f"{outputfolder}/log_verbose.txt", 'w') as log_verbose:
    #         with open(f"{outputfolder}/log_error.txt", 'w') as log_error:
    #             subprocess.run([script_path+'./submit_job_automated_new_multiple_vcfs.sh', str(genomedir), '_', str(outputfolder)+"/guides.txt", str(pamfile), str(annotationfile), str(script_path+'vuoto.txt'),
    #                             str(bMax), str(mm), str(bDNA), str(bRNA), str(merge_t), str(outputfolder), str(script_path), str(thread), str(current_working_directory), str(gene_annotation),void_mail,str(base_start),str(base_end),str(base_set)], stdout=log_verbose, stderr=log_error)
    # change name of guide and param files to hidden
    os.system(f"mv {outputfolder}/guides.txt {outputfolder}/.guides.txt")
    os.system(f"mv {outputfolder}/Params.txt {outputfolder}/.Params.txt")


def target_integration():
    if "--help" in input_args:
        print(
            "This is the automated integration process that process the final result file to generate a usable target panel."
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--targets, used to specify the final result file to use in the panel creation process"
        )
        print(
            "\t--empirical_data, used to specify the file that contains empirical data provided by the user to assess in-silico targets"
        )
        print("\t--output, used to specify the output folder for the results")
        exit(0)

    if "--targets" not in input_args:
        print("--targets must be contained in the input")
        exit(1)
    else:
        try:
            target_file = os.path.abspath(input_args[input_args.index("--targets") + 1])
        except IndexError:
            print("Please input some parameter for flag --targets")
            exit(1)
        if not os.path.isfile(target_file):
            print("The file specified for --target_file does not exist")
            exit(1)

    # if "--vcf_dir" not in input_args:
    #     print("--vcf_dir non in input, multi-variant haplotype will not be calculated")
    #     vcf_dir = script_path+'vuota/'
    #     # exit(1)
    # else:
    #     try:
    #         vcf_dir = os.path.abspath(
    #             input_args[input_args.index("--vcf_dir")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --vcf_dir")
    #         exit(1)
    #     if not os.path.isdir(vcf_dir):
    #         print("The folder specified for --vcf_dir does not exist")
    #         exit(1)

    # if "--genome_version" not in input_args:
    #     print("--genome_version must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         genome_version = input_args[input_args.index(
    #             "--genome_version")+1]
    #     except IndexError:
    #         print("Please input some parameter for flag --genome")
    #         exit(1)

    # if "--guide" not in input_args:
    #     guidefile = script_path+'vuoto.txt'
    #     # print("--guide must be contained in the input")
    #     # exit(1)
    # else:
    #     try:
    #         guidefile = os.path.abspath(
    #             input_args[input_args.index("--guide")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --guide")
    #         exit(1)
    #     if not os.path.isfile(guidefile):
    #         print("The file specified for --guide does not exist")
    #         exit(1)

    if "--empirical_data" not in input_args:
        print("--empirical_data not in input, proceeding without empirical data")
        empiricalfile = script_path + "vuoto.txt"
        # exit(1)
    else:
        try:
            empiricalfile = os.path.abspath(
                input_args[input_args.index("--empirical_data") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --empirical_data")
            exit(1)
        if not os.path.isfile(empiricalfile):
            print("The file specified for --empirical_data does not exist")
            exit(1)

    # if "--gencode" not in input_args:
    #     print("--gencode must be contained in the input")
    #     exit(1)
    # else:
    #     try:
    #         gencode_file = os.path.abspath(
    #             input_args[input_args.index("--gencode")+1])
    #     except IndexError:
    #         print("Please input some parameter for flag --gencode")
    #         exit(1)
    #     if not os.path.isfile(gencode_file):
    #         print("The file specified for --gencode does not exist")
    #         exit(1)

    if "--output" not in input_args:
        print("--output must be contained in the input")
        exit(1)
    else:
        try:
            outputfolder = os.path.abspath(input_args[input_args.index("--output") + 1])
        except IndexError:
            print("Please input some parameter for flag --output")
            exit(1)
        if not os.path.isdir(outputfolder):
            print("The folder specified for --output does not exist")
            exit(1)

    os.system(
        f"{script_path}./empirical_integrator.py {target_file} {empiricalfile} {outputfolder}"
    )


def print_help_gnomad_converter():
    """
    Prints the help information for the gnomAD converter functionality, providing
    details on the conversion process from gnomAD VCFs to VCFs compatible with
    CRISPRme. It outlines the options available for specifying directories,
    sample IDs, variant filtering, multiallelic site handling, and thread usage
    during the conversion process.

    Raises:
        SystemExit: If the help information is displayed to guide users on using
        the gnomAD converter functionality.
    """

    # functionality description
    sys.stderr.write(
        "The gnomAD converter functionality simplifies the conversion process "
        "of gnomAD VCFs (versions 3.1 and 4.0) into VCFs supported by CRISPRme. "
        "It ensures a seamless transition while maintaining compatibility with "
        "CRISPRme's requirements, focusing on the structure and content of "
        "precomputed sample IDs file \n\n"
    )
    # options
    sys.stderr.write(
        "Options:\n"
        "\t--gnomAD_VCFdir, specifies the directory containing gnomAD VCFs. "
        "Files must have the BGZ extension\n"
        "\t--samplesID, specifies the precomputed sample IDs file necessary "
        "for incorporating population-specific information into the output "
        "VCFs\n"
        "\t--joint, optional flag to specify the input GnomAD VCF contain joint "
        "allele frequencies\n"
        "\t--keep, optional flag to retain all variants, regardless of their "
        "filter flag. By default, variants with a filter flag different from "
        "PASS are discarded\n"
        "\t--multiallelic, optional flag to merge variants mapped to the "
        "same position, creating multiallelic sites in the output VCFs. By "
        "default, each site remains biallelic\n"
        "\t--thread, used to set the number of thread used in the conversion "
        "process [default 8]\n"
    )
    sys.exit(1)


def gnomAD_converter():
    """
    Runs the gnomAD converter functionality based on specified arguments, converting
    gnomAD VCF files into formats compatible with CRISPRme.

    Raises:
        ValueError: If mandatory arguments are missing or have incorrect values.
        FileExistsError: If the specified gnomAD VCF directory cannot be located.
        FileNotFoundError: If the specified sample IDs file cannot be found.
        subprocess.SubprocessError: If an error occurs during the gnomAD VCF
            conversion process.
    """

    args = input_args[2:]  # recover gnomAD converter args
    if "--help" in args or not args:  # print help
        print_help_gnomad_converter()
        sys.exit(1)
    if "--gnomAD_VCFdir" not in args:
        raise ValueError(
            "--gnomAD_VCFdir is a mandatory argument required for the conversion "
            "process. Please specify the directory containing gnomAD VCFs using "
            "this option\n"
        )
    if "--samplesID" not in args:
        raise ValueError(
            "--samplesID is a mandatory argument required for the conversion "
            "process. Please specify the sample IDs file this option\n"
        )
    # read gnomAD directory arg
    try:
        gnomad_dir = args[args.index("--gnomAD_VCFdir") + 1]
        if gnomad_dir.startswith("--"):
            raise ValueError("Please input some parameter for flag --gnomAD_VCFdir\n")
        gnomad_dir = os.path.abspath(gnomad_dir)  # first sanity check passed
        if not os.path.isdir(gnomad_dir):
            raise FileExistsError(f"Unable to locate {gnomad_dir}")
    except IndexError as e:
        raise ValueError(
            "Please input some parameter for flag --gnomAD_VCFdir\n"
        ) from e
    # read samples ids arg
    try:
        samples_ids = args[args.index("--samplesID") + 1]
        if samples_ids.startswith("--"):
            raise ValueError("Please input some parameter for flag --samplesID")
        samples_ids = os.path.abspath(samples_ids)  # first sanity check passed
        if not os.path.isfile(samples_ids):
            raise FileNotFoundError(f"Unable to locate {samples_ids}")
    except IndexError as e:
        raise ValueError("Please input some parameter for flag --samplesID") from e
    # read joint gnomad vcf files
    joint = "--joint" in args
    # read keep arg
    keep = "--keep" in args  # keep all variants regardless of filter label
    # read multiallelic arg
    multiallelic = "--multiallelic" in args  # merge variants in multiallelic sites
    # read threads arg
    threads = 8
    if "--threads" in args:
        try:
            threads = int(args[args.index("--threads") + 1])
            if threads <= 0:
                raise ValueError(f"Forbidden number of threads ({threads})")
        except IndexError as e:
            raise ValueError("Missing or forbidden threads value") from e
    # run gnom AD converter
    gnomad_converter_script = os.path.join(script_path, "convert_gnomAD_vcfs.py")
    cmd = (
        f"python {gnomad_converter_script} {gnomad_dir} {samples_ids} {joint} "
        f"{keep} {multiallelic} {threads}"
    )
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise subprocess.SubprocessError(
            f"An error occurred while converting gnomAD VCFs in {gnomad_dir}"
        )


def personal_card():
    if "--help" in input_args:
        print(
            "This is the personal card generator that creates a files with all the private targets for the input sample"
        )
        print("These are the flags that must be used in order to run this function:")
        print(
            "\t--result_dir, directory containing the result from which extract the targets to generate the card"
        )
        print(
            "\t--guide_seq, sequence of the guide to use in order to exctract the targets"
        )
        print("\t--sample_id, ID of the sample to use in order to generate the card")
        exit(0)

    if "--result_dir" not in input_args:
        print("--result_dir not in input, please input a result directory")
        exit(1)
    else:
        try:
            result_dir = os.path.abspath(
                input_args[input_args.index("--result_dir") + 1]
            )
        except IndexError:
            print("Please input some parameter for flag --result_dir")
            exit(1)
        if not os.path.isdir(result_dir):
            print("The folder specified for --result_dir does not exist")
            exit(1)

    if "--guide_seq" not in input_args:
        print(
            "--guide_seq must be contained in the input, e.g. CTAACAGTTGCTTTTATCACNNN"
        )
        exit(1)
    else:
        try:
            guide = input_args[input_args.index("--guide_seq") + 1]
        except IndexError:
            print("Please input some parameter for flag --guide_seq")
            exit(1)
    if "--sample_id" not in input_args:
        print("--sample_id must be contained in the input, e.g. HG00001")
        exit(1)
    else:
        try:
            sample_id = input_args[input_args.index("--sample_id") + 1]
        except IndexError:
            print("Please input some parameter for flag --sample_id")
            exit(1)

    os.system(
        script_path
        + "./generate_sample_card.py "
        + result_dir
        + " "
        + guide
        + " "
        + sample_id
        + " "
        + script_path
    )


def web_interface():
    if "--help" in input_args:
        print(
            "This function must be launched without input, it starts a local server to use the web interface."
        )
        print(
            "Open your web-browser and write 127.0.0.1:8080 in the search bar if you are executing locally, if you are executing on an external server write <yourserverip>:8080 in search bar"
        )
        exit(0)
    subprocess.run(corrected_web_path + "/./index.py")


def crisprme_version():
    if len(input_args) != 2:
        sys.stderr.write("Wrong number of arguments for crisprme.py version\n")
        sys.exit(1)
    sys.stdout.write(f"v{__version__}\n")


def print_help_complete_test():
    """
    Prints the help information for executing comprehensive testing of the
    complete-search functionality provided by CRISPRme.

    Raises:
        SystemExit: If the help information is displayed to guide users on
            executing comprehensive testing.
    """

    # write intro message to stdout
    sys.stderr.write(
        "Execute comprehensive testing for complete-search functionality "
        "provided by CRISPRme\n"
    )
    # list functionality options
    sys.stderr.write(
        "Options:\n"
        "\t--chrom, test the complete-search functionality on the specified "
        "chromosome (e.g., chr22). By default, the test is conducted on all "
        "chromosomes\n"
        "\t--vcf_dataset, VCFs dataset to be used during CRISPRme testing. "
        "Available options include 1000 Genomes (1000G) and Human Genome "
        "Diversity Project (HGDP). To use the combined dataset type '1000G+HGDP' "
        "The default dataset is 1000 Genomes.\n"
        "\t--thread, number of threads.\n"
        "\t--debug, debug mode.\n"
    )
    sys.exit(1)


def complete_test_crisprme():
    """
    Executes comprehensive testing for the complete-search functionality provided
    by CRISPRme based on specified arguments.

    Raises:
        OSError: If the CRISPRme test fails, indicating an issue with the testing
            process.
    """

    if "--help" in input_args or len(input_args) < 3:
        print_help_complete_test()
        sys.exit(1)
    chrom = "all"
    if "--chrom" in input_args:  # individual chrom to test
        try:
            chrom = input_args[input_args.index("--chrom") + 1]
            if chrom.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --chrom\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --chrom\n")
            sys.exit(1)
    vcf_dataset = "1000G"
    if "--vcf_dataset" in input_args:  # specified variant dataset
        try:
            vcf_dataset = input_args[input_args.index("--vcf_dataset") + 1]
            if vcf_dataset.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some parameter for flag --vcf_dataset\n")
            sys.exit(1)
    threads = 4
    if "--thread" in input_args:  # number of threads to use during test
        try:
            threads = input_args[input_args.index("--thread") + 1]
            if threads.startswith("--"):
                sys.stderr.write("Please input some parameter for flag --thread\n")
                sys.exit(1)
        except IndexError:
            sys.stderr.write("Please input some value for flag --thread\n")
            sys.exit(1)
    debug = "--debug" in input_args  # run local or via conda/Docker
    # begin crisprme test
    script_test = os.path.join(script_path, "complete_test.py")
    code = subprocess.call(
        f"python {script_test} {chrom} {vcf_dataset} {threads} {debug}", shell=True
    )
    if code != 0:
        raise OSError(
            "\nCRISPRme test failed! See Results/crisprme-test-out/log_error.txt for details\n"
        )


# HELP FUNCTION
def crisprme_help() -> None:
    """
    Prints the general help information for CRISPRme, describing each available 
    functionality.

    Outputs usage instructions, requirements for input files, and a summary of 
    all main commands to stderr, then exits the program.
    """
    # print crisprme help; describe each functionality
    sys.stderr.write(
        "Help:\n\n"
        "- ALL FASTA FILEs USED BY THE SOFTWARE MUST BE UNZIPPED AND SEPARATED BY CHROMOSOME\n"
        "- ALL VCFs USED BY THE SOFTWARE MUST BE ZIPPED (WITH BGZIP) AND SEPARATED BY CHROMOSOME\n\n"
        "Functionalities:\n\n"
        "crisprme.py complete-search\n"
        "\tPerforms genome-wide off-targets search (reference and variant, if "
        "specified), including CFD and CRISTA analysis, and target selection\n\n"
        "crisprme.py complete-test\n"
        "\tTest the complete CRISPRme pipeline on single chromosomes or complete "
        "genomes\n\n"
        "crisprme.py targets-integration\n"
        "\tIntegrates in-silico targets with empirical data to generate a usable "
        "panel\n\n"
        "crisprme.py gnomAD-converter\n"
        "\tConverts gnomAD VCF files into CRISPRme compatible VCFs (supports "
        "gnomAD >= v3.1)\n\n"
        "crisprme.py generate-personal-card\n"
        "\tGenerates a personal card for specific samples by extracting all "
        "private targets\n\n"
        "crisprme.py web-interface\n"
        "\tActivates CRISPRme's web interface for local browser use\n\n"
        "crisprme.py --version\n"
        "\tPrints CRISPRme version to stdout and exit\n\n"
        "For additional information on each CRISPRme functionality type <function> "
        "--help (e.g. 'crisprme.py complete-search --help')\n"
    )
    sys.exit(1)  # stop execution


if len(sys.argv) < 2:
    check_crisprme_dirtree()  # check crisprme directory tree structure
    crisprme_help()  # no arg? print help
elif sys.argv[1] == "complete-search":  # run complete search
    complete_search()
elif sys.argv[1] == "complete-test":  # run complete test
    complete_test_crisprme()
elif sys.argv[1] == "targets-integration":  # run targets integration
    target_integration()
elif sys.argv[1] == "gnomAD-converter":  # run gnomad converter
    gnomAD_converter()
elif sys.argv[1] == "generate-personal-card":  # run create personal card
    personal_card()
elif sys.argv[1] == "web-interface":  # run web interface
    web_interface()
elif sys.argv[1] == "--version":  # print version
    crisprme_version()
else:
    sys.stderr.write(f"ERROR! {sys.argv[1]} is not an allowed!\n\n")
    crisprme_help()  # print help if invalid command is given
