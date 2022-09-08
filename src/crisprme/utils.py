"""Utils functions and variables used throughout CRISPRme code
"""


from argparse import Namespace
from contextlib import redirect_stdout, redirect_stderr
from typing import List, NoReturn, Optional, Tuple
from colorama import Fore, init
from Bio.Seq import Seq
from io import TextIOWrapper

import pybedtools
import subprocess
import itertools
import sys
import re
import os


# ------------------------------------------------------------------------------
# constant variables 
# SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
CONDA_PATH = "opt/crisprme/src/crisprme/PostProcess"
ORIGIN_PATH = os.path.dirname(
    os.path.join(os.path.abspath(__file__), "../..")
)[:-3] + CONDA_PATH
# WEB_PATH = ORIGIN_PATH[:-3] + "opt/crisprme"
CURRENT_WORKING_DIRECTORY = os.path.join(os.getcwd())
CRISPRme_COMMANDS = [
    "complete-search", 
    "gnomAD-converter", 
    "targets-integration", 
    "web-interface",
    "generate-personal-card"
]
CRISPRme_DIRS = [
    "Genomes", "Results", "Dictionaries", "VCFs", "Annotations", "PAMs", "samplesIDs"
]
IUPAC_ALPHABET = {
    "A",
    "T",
    "C",
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
    "a",
    "t",
    "c",
    "g",
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
PAM_DICT = {
        'A':  "ARWMDHV",
        'C':  "CYSMBHV",
        'G':  "GRSKBDV",
        'T':  "TYWKBDH",
        'R':  "ARWMDHVSKBG",
        'Y':  "CYSMBHVWKDT",
        'S':  "CYSMBHVKDRG",
        'W':  "ARWMDHVYKBT",
        'K':  "GRSKBDVYWHT",
        'M':  "ARWMDHVYSBC",
        'B':  "CYSMBHVRKDGWT",
        'D':  "ARWMDHVSKBGYT",
        'H':  "ARWMDHVYSBCKT",
        'V':  "ARWMDHVYSBCKG",
        'N':  "ACGTRYSWKMBDHV",
    }
LOG_FILE = "log.txt"
LOG_VERBOSE = "log_verbose.txt"
LOG_ERROR = "log_error.txt"
QUEUE_FILE = "queue.txt"


# ------------------------------------------------------------------------------
# utils functions
def sigint_handler() -> NoReturn:
    """Catch SIGINT and exit.
    """

    sys.stderr.write("\nSIGINT Caught. CRISPRme will exit.")
    sys.exit(3)  # ERRCODE 3 --> Keyboard interrupt


def exception_handler(
    exception_type: Exception, exception: str, debug: Optional[bool] = False
) -> NoReturn:
    """Handle exceptions (for debugging purposes).
    If `debug` is set to True, print the entire error stack to stderr, otherwise
    send the error message and gracefully exit.

    ...

    Parameters
    ----------
    exception_type : Exception
        Exception type
    exception : str
        Exception message
    debug : bool

    Returns
    -------
    NoReturn
    """

    init()  # initialize colorama
    if debug:
        raise exception_type(f"\n\n{exception}")
    else:
        sys.stderr.write(Fore.RED + f"\n\nERROR: {exception}\n" + Fore.RESET)
        sys.exit(1)  # ERRCODE 1 --> Exception raised


def check_directories_consistency(debug: bool) -> None:
    """Check CRISPRme directory tree status.
    If some directopry is missing, create it.

    ...

    Parameters
    ----------
    debug : bool

    Returns
    -------
    None
    """

    for d in CRISPRme_DIRS:
        # check CRISPRme directories existence in the current working dir
        if not os.path.exists(os.path.join(CURRENT_WORKING_DIRECTORY, d)):
            try:
                # try to create the missing directories
                os.makedirs(os.path.join(CURRENT_WORKING_DIRECTORY, d))
            except:
                exception_handler(
                    OSError, 
                    f"Directory {d} is missing. CRISPRme was unable to create it.",
                    debug
                )


def check_command_args(command: str, args: Namespace) -> bool:
    if command not in CRISPRme_COMMANDS:
        exception_handler(
            ValueError,
            f"Forbidden CRISPRme command found ({command})",
            True  # always trace this kind of errors
        )
    # check arguments consistency for each command
    args_dict = vars(args)
    if command != CRISPRme_COMMANDS[0]:  # complete-search
        # make sure that only complete-search has such args initialized
        if args_dict["genome"] != os.path.join(
            CURRENT_WORKING_DIRECTORY, "Genomes"
        ):
            return True
        if args_dict["vcf"] != "_":
            return True
        if bool(args_dict["sequence_guides"]):
            return True
        if bool(args_dict["pam"]):
            return True
        if bool(args_dict["be_window"]):
            return True
        if bool(args_dict["be_nucleotide"]):
            return True
        if bool(args_dict["annotation"]):
            return True
        if bool(args_dict["personal_annotation"]):
            return True
        if bool(args_dict["samples_file"]):
            return True
        if bool(args_dict["gene_annotation"]):
            return True
        if bool(args_dict["mm"]):
            return True
        if args_dict["bdna"] != 0:
            return True
        if args_dict["brna"] != 0:
            return True
        if args_dict["merge"] != 3:
            return True
        if bool(args_dict["output_name"]):
            return True
    if command != CRISPRme_COMMANDS[1]:  # gnomAD-converter
        # make sure that only gnomAD-converter has such args initialized
        if bool(args_dict["vcf_dir"]):
            return True
        if bool(args_dict["samples_file_gnomAD"]):
            return True
    if command != CRISPRme_COMMANDS[2]:  # targets-integration
        # make sure that only targets-integration has such args initialized
        if bool(args_dict["targets_file"]):
            return True
        if bool(args_dict["empirical_data"]):
            return True
        if args_dict["outdir"] != os.getcwd():
            return True
    if command != CRISPRme_COMMANDS[3]:  # web-interface
        # make sure that only web-interface has such args initialized
        if args_dict["help_web"]:
            return True
    if command != CRISPRme_COMMANDS[4]:  # generate-personal-card
        # make sure that only generate-personal-card has such args initialized
        if bool(args_dict["inputdir"]):
            return True
        if bool(args_dict["guide_seq"]):
            return True
        if bool(args_dict["sample_id"]):
            return True
    return False


def merge_annotation_files(
    func_annotation: str, personal_annotation: str, debug: bool
) -> str:
    """Merge functional and personal annotations in a single file.

    ...

    Parameters
    ----------
    func_annotation : str
        Functional annotation file
    personal_annotation : str
        Personal annotation file
    debug : bool

    Returns
    -------
    str
        Path to merged annotation file
    """
    
    # type and file existence checks already made
    cmd = f'awk \'$4 = $4\"_personal\"\' {personal_annotation} | sed "s/ /\t/g" | sed "s/,/_personal,/g" > {personal_annotation}.tmp'
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        exception_handler(
            subprocess.SubprocessError, 
            f"An error occurred while running {cmd}",
            debug
        )
    annotation_merge = f"{func_annotation}+personal.bed"
    cmd = f"cat {personal_annotation}.tmp {func_annotation} > {annotation_merge}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        exception_handler(
            subprocess.SubprocessError,
            f"An error occurred while running {cmd}",
            debug
        )
    cmd = f"rm -f {personal_annotation}.tmp"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        exception_handler(
            subprocess.SubprocessError,
            f"An error occurred while running {cmd}",
            debug
        )
    return annotation_merge


def parse_PAM_sequence_file(sequence_file: str, debug: bool) -> Tuple:
    """Parse the PAM file and recover the sequences.

    ...

    Parameters
    ----------
    sequence_file : str
        PAM sequence file
    debug : bool

    Returns
    -------
    Tuple
    """

    if not isinstance(sequence_file, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(sequence_file).__name__}",
            debug
        )
    if not os.path.isfile(sequence_file):
        exception_handler(
            FileNotFoundError,
            f"Unable to locate {sequence_file}",
            debug
        )
    try:
        # read the PAM file to recover the PAM sequence and start position
        with open(sequence_file, mode="r") as handle:
            full_pam, pam_value_idx = handle.readline().strip().split()
            full_pam_length = len(full_pam)  # full PAM length
            try:
                pam_value_idx = int(pam_value_idx)  # must be integer value
            except ValueError as verr:
                exception_handler(
                    ValueError,
                    f"Forbidden PAM character index found ({pam_value_idx})",
                    debug
                )
            except Exception as e:
                exception_handler(
                    Exception, "Unable to read PAM caharacter index", debug
                )
            # the PAM is placed at the end of the full sequence -> the start 
            # position is based on the end of the sequence
            # E.g. NNNNNNNGG -> PAM == NGG;  start == 3 
            if pam_value_idx < 0:
                pam_seq = full_pam[:-pam_value_idx]
                pam_length = -pam_value_idx
                pam_start = True
            else:
                pam_seq = full_pam[-pam_value_idx:]
                pam_length = pam_value_idx
                pam_start = False
    except:
        exception_handler(
            OSError, f"An error occurred while reading {sequence_file}", debug
        )
    return pam_seq, pam_length, full_pam, full_pam_length, pam_start


def parse_guide_sequences_file(
    sequences_file: str, 
    genome: str, 
    pam: str, 
    pam_length: int,
    guide_length: int, 
    pam_start: bool,
    debug: bool
) -> List:
    """Parse the file containing guide sequences.
    
    ...

    Parameters
    ----------
    sequences_file : str
        File containing guide sequences
    genome : str
        Genome location
    pam : str
        PAM sequence
    pam_length : int
        PAM sequence length
    guide_length : int
        Guide length
    pam_start : bool
        PAM start flag
    debug : bool

    Returns
    -------
    List
    """

    if not isinstance(sequences_file, str):
        exception_handler(
            TypeError,
            f"Exepcted {str.__name__}, got {type(sequences_file).__name__}",
            debug
        )    
    if not os.path.isfile(sequences_file):
        exception_handler(
            FileNotFoundError,
            f"Unable to locate {sequences_file}",
            debug
        )
    if not isinstance(genome, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(genome).__name__}",
            debug
        )
    if not os.path.exists(genome):
        exception_handler(
            FileNotFoundError, f"Unable to locate {genome}", debug
        )
    try:
        # read all sequence names and sequences in the file and create a 
        # single str storing them
        with open(sequences_file, mode="r") as handle:
            sequences_raw = "".join([line for line in handle])
    except OSError as oserr:
        exception_handler(
            OSError, f"Error while reading {sequences_file}", debug
        )
    guides = []  # guides list
    # recover seqnames and sequences
    for seqname_and_seq in sequences_raw.split(">"):
        if not seqname_and_seq:
            continue
        i = seqname_and_seq.find("\n")
        seqname = seqname_and_seq[:i]
        seq = seqname_and_seq[i:].strip()
        if "chr" in seq:  # BED file feature
            for line in seq.split("\n"):
                if not line:
                    continue
                line_split = line.strip().split()
                coords = f"{line_split[0]}:{line_split[1]}-{line_split[2]}"
                sequence = extract_sequence(
                    seqname, coords, genome.replace(" ", "_"), debug
                )
                guides.append(get_guides(sequence, pam, guide_length, pam_start))
        else:  # not BED file but simple guide sequences
            sequence = "".join(seq.split()).strip()
            guides.append(get_guides(sequence, pam, guide_length, pam_start))
    guides_tmp = []  # auxiliary guides list
    for guide in guides:
        n_prefix = "N" * pam_length  
        if pam_start:
            # add as many Ns as long is the PAM in front of the guide
            guides_tmp.append("".join([n_prefix, guide]))
        else:
            # append as many Ns as long is the PAM on the back of the guide
            guides_tmp.append("".join([guide, n_prefix]))
    # limit guides number 
    if len(guides_tmp) > 1000000000:
        guides_tmp = guide[:1000000000]
    guides = guides_tmp  # final guides considered during the search
    return guides


def extract_sequence(
    seqname: str, coords: str, genome: str, debug: bool
) -> str:
    """Extract genomic sequences in the genome region defined by the input
    coordinates.

    ...

    Parameters
    ----------
    seqname : str
        Sequence name
    coords : str
        Genomic coordinates
    genome : str
        Genome
    debug : bool
    
    Returns
    -------
    str
    """

    if not isinstance(seqname, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(seqname).__name__}",
            debug
        )
    if not isinstance(coords, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(coords).__name__}",
            debug
        )
    if not isinstance(genome, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(genome).__name__}",
            debug
        )
    if not os.path.exists(genome):
        exception_handler(FileNotFoundError, f"Unable to locate {genome}", debug)
    seqname = "_".join(seqname.split())
    chrom, positions = coords.split(":")
    start, stop = positions.split("-")
    start = start.replace(",", "").replace(".", "").replace(" ", "")
    stop = stop.replace(",", "").replace(".", "").replace(" ", "")
    chroms = [
        f for f in os.listdir(
            os.path.join(CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[0], genome)
        )
        if os.path.isfile(
            os.path.join(CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[0], genome, f)
        )
        and not f.startswith(".")  # skip hidden files
        and not f.endswith(".fai")
    ]
    fa_ext = "fasta" if ".fasta" in chroms[0] else "fa"
    # use pybedtools to extract sequences 
    seq = pybedtools.BedTool(f"{chrom} {start} {stop}", from_string=True)
    fastafile = os.path.join(
        CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[0], genome, f"{chrom}.{fa_ext}"
    )
    seq = seq.sequence(fi=fastafile)
    sequence = open(seq.seqfn).read().split("\n")[1].strip()
    # remove fai index file if it has been created 
    try:
        os.remove(
            os.path.join(
                CURRENT_WORKING_DIRECTORY, CRISPRme_DIRS[0], genome, f"{chrom}.fa.fai"
            )
        )
    except: 
        pass  # if fai not found, just skip
    return sequence


def get_guides(
    sequence: str, pam: str, guide_length: int, pam_start: int, debug: bool
) -> List:
    """Recover guide sequences (both forward and reverse).

    ...

    Parameters
    ----------
    sequence : str
        Base sequence
    pam : str
        PAM sequence
    guide_length : int
        Guide length
    pam_start : bool
        PAM start flag
    debug : bool

    Returns
    -------
    List
    """

    if not isinstance(sequence, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(sequence).__name__}",
            debug
        )
    if not isinstance(pam, str):
        exception_handler(
            TypeError,
            f"Expected {str.__name__}, got {type(pam).__name__}",
            debug
        )
    if not isinstance(guide_length, int):
        exception_handler(
            TypeError,
            f"Expected {int.__name__}, got {type(guide_length).__name__}",
            debug
        )
    if not isinstance(pam_start, int):
        exception_handler(
            TypeError,
            f"Expected {int.__name__}, got {type(pam_start).__name__}",
            debug
        )
    pam_length = len(pam)  # get PAM length
    # compute forward and reverse PAMs in IUPAC code
    pam_expanded = [PAM_DICT[c] for c in pam]
    pam_iupac = ["".join(e) for e in itertools.product(*pam_expanded)]
    pam_reverse = str(Seq(pam).reverse_complement())
    pam_expanded_rev = [PAM_DICT[c] for c in pam_reverse]
    pam_iupac_rev = ["".join(e) for e in itertools.product(*pam_expanded_rev)]
    sequence = sequence.upper()  # force sequeence nucleotides to be upper case
    sequence_len = len(sequence)  # get sequence length
    guides = []  # guides list
    # extract guides (on forward)
    for pam in pam_iupac:
        pos = ([m.start() for m in re.finditer('(?=' + pam + ')', sequence)])
        if pos:
            for i in pos:   
                # i is the position where the first char of pam is found
                # E.g. N in NNNNNN NGG
                if pam_start:
                    if i > (sequence_len - guide_length - pam_length):
                        continue
                    guides.append(
                        sequence[(i + pam_length):(i + pam_length + guide_length)]
                    )
                else:
                    if i < guide_length:
                        continue
                    guides.append(sequence[(i - guide_length):i])
    # extract guides (on reverse)
    for pam in pam_iupac_rev:
        pos = ([m.start() for m in re.finditer('(?=' + pam + ')', sequence)])
        if pos:
            for i in pos:
                if pam_start:
                    if i < guide_length:
                        continue
                    guides.append(
                        str(
                            Seq(
                                sequence[(i - guide_length):i]
                            ).reverse_complement()  # compute reverse guide
                        )
                    )
                else:
                    if i > (sequence_len - guide_length - pam_length):
                        continue
                    guides.append(
                        str(
                            Seq(
                                sequence[
                                    (i + pam_length):(i + guide_length + pam_length)
                                ]
                            ).reverse_complement()  # compute reverse guide
                        )
                    )
    assert bool(guides)  # check that at least one guide  has been found
    return guides  # return guides for when adding to app.py


def redir_stdout(stream: TextIOWrapper, msg: str) -> None:
    """Redirect the stout to the given stream channel.

    ...

    Parameters
    ----------
    stream : TextIOWrapper
        Stream to which the stdout is redicrected
    msg : str
        Message to redirect

    Returns
    -------
    None
    """

    assert isinstance(msg, str)
    try:
        with redirect_stdout(stream):
            sys.stdout.write("".join([msg, "\n"]))
    except:
        # always trace back such errors
        exception_handler(
            IOError, "An error occurred while redirecting stdout.", True
        )


def redir_stderr(stream: TextIOWrapper, msg: str) -> None:
    """Redirect the sterr to the given stream channel.

    ...

    Parameters
    ----------
    stream : TextIOWrapper
        Stream to which the stderr is redicrected
    msg : str
        Message to redirect

    Returns
    -------
    None
    """

    assert isinstance(msg, str)
    try:
        with redirect_stderr(stream):
            sys.stderr.write("".join([msg, "\n"]))
    except:
        # always trace back such kind of errors
        exception_handler(
            IOError, "An error occurred while redicrecting stderr.", True
        )
               
