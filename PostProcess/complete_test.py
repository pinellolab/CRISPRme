"""
This module provides functionality to execute the CRISPRme test workflow, including
downloading genomic and variant data, preparing input files, and running the CRISPRme
command-line tool. It includes functions for managing directories, downloading data
from various sources, and configuring test parameters.

Key functions include:
- `_assign_genome_directory_name`: Resolves the genome subdirectory name for the 
    requested chromosome selection.
- `download_genome_data`: Downloads genome data for a specified chromosome to the
    destination directory.
- `ensure_vcf_dataset_directory`: Ensures the existence of a directory for a specific
    VCF dataset.
- `download_vcf_data`: Downloads VCF data for a specific chromosome and variant dataset.
- `ensure_samplesids_directory`: Ensures the existence of the 'samplesIDs' directory.
- `download_samples_ids_data`: Downloads samples IDs data for a specific variant dataset.
- `ensure_annotation_directory`: Ensures the existence of the 'annotation' directory.
- `download_annotation_data`: Downloads gencode and encode annotation data to the
    'annotation' directory.
- `write_ngg_pamfile`: Writes a test PAM file containing the NGG sequence.
- `write_sg1617_guidefile`: Writes a test guide file containing the sg1617 guide
    sequence.
- `write_vcf_config`: Writes a test VCF list file for a specific variant dataset.
- `write_samplesids_config`: Writes a test samples ID list file for a specific
    variant dataset.
- `run_crisprme_test`: Executes the CRISPRme test workflow for a specified chromosome
    and dataset.
- `main`: The entry point of the module that orchestrates the test execution.

This module is designed to facilitate the testing and validation of the CRISPRme
tool, ensuring that all necessary data and configurations are correctly handled
before running the analysis.
"""

from utils import (
    check_crisprme_directory_tree,
    download,
    gunzip,
    untar,
    rename,
    compute_md5,
    CHROMS,
    CRISPRME_DIRS,
    MD5GENOME,
    MD51000G,
    MD5HGDP,
    MD5SAMPLES,
    MD5ANNOTATION,
)

from typing import Tuple

import subprocess
import json
import sys
import os

# benchmark registry (test/benchmark/benchmarks.json)
BENCHMARKS_JSON = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, "test", "benchmark",
    "benchmarks.json",
)

# define genome data url
HG38URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"

# define 1000G server and url
VCF1000GSERVER = "ftp.1000genomes.ebi.ac.uk"
VCF1000GURL = (
    "/vol1/ftp/data_collections/1000_genomes_project/release/"
    "20190312_biallelic_SNV_and_INDEL/"
    "ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)

# define hgdp server and url
VCFHGDPSERVER = "ngs.sanger.ac.uk"
VCFHGDPURL = "/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{}.vcf.gz"

# define test data url
TESTDATAURL = (
    "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/test/data/"
)

# define complete-test results folder name
COMPLETETESTRESDIR = "crisprme-test-out"


def check_output() -> None:
    """
    Check whether complete-test results already exist and prevent rerunning tests
    if previous output is detected.

    This function inspects the expected results directory and exits with an error
    message if a complete-test output folder is found, avoiding accidental
    overwriting or conflicting test runs.

    Raises:
        SystemExit: If a previous complete-test results folder is found.
        AssertionError: If the main results directory does not exist.
    """
    results_dir = os.path.abspath(os.path.join(os.getcwd(), CRISPRME_DIRS[1]))
    assert os.path.isdir(results_dir)
    # one output dir per registered benchmark (crisprme-test-out_<name>)
    for bench in load_benchmarks()["benchmarks"]:
        d = os.path.join(results_dir, f"{COMPLETETESTRESDIR}_{bench['name']}")
        if os.path.isdir(d):
            sys.stderr.write(
                "Complete-test already run once. Please delete the complete-test "
                f"results folder before running it again: {d}\n"
            )
            sys.exit(0)  # avoid throwing complete-search error on output folder


def _assign_genome_directory_name(chrom: str) -> str:
    """Return the genome subdirectory name for a chromosome selection.

    The whole-genome build (``chrom == "all"``) lives in ``hg38``; a
    single-chromosome build lives in ``hg38_<chrom>`` so per-chromosome test
    assets stay isolated from a full build.

    Args:
        chrom (str): Chromosome in UCSC format (e.g. "chr22"), or "all".

    Returns:
        str: Genome subdirectory name (e.g. "hg38" or "hg38_chr22").
    """
    return "hg38" if chrom == "all" else f"hg38_{chrom}"


def download_genome_data(chrom: str, dest: str) -> str:
    """
    Download genome data for the requested chromosome selection and return the
    directory containing the resulting FASTA file(s).

    The whole-genome build (``chrom == "all"``) is stored in ``<dest>/hg38``.
    A single-chromosome build is stored in ``<dest>/hg38_<chrom>``.

    Note on coordinates: this function only handles file placement; the FASTA
    payload is UCSC hg38 (1-based, closed genomic coordinates) and is not
    modified here.

    Args:
        chrom (str): Chromosome in UCSC format (e.g. "chr22"), or "all".
        dest (str): The CRISPRme "Genomes" directory in which to create the
            genome directory.

    Returns:
        str: Absolute path to the genome directory to pass to
            ``crisprme.py complete-search --genome``.

    Raises:
        ValueError: If the chromosome is invalid or an MD5 check fails.
        FileExistsError: If the destination directory does not exist.
    """
    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    sys.stderr.write(f"Downloading fasta file for chromosome(s) {chrom}\n")
    genome_dir = os.path.join(dest, _assign_genome_directory_name(chrom))
    if chrom == "all":
        chromstar = download(dest, http_url=f"{HG38URL}/bigZips/hg38.chromFa.tar.gz")
        # check genome md5
        if MD5GENOME[os.path.basename(chromstar)] != compute_md5(chromstar):
            raise ValueError(f"Download for {os.path.basename(chromstar)} failed")
        chromsdir = untar(chromstar, dest, "chroms")  # decompress archive
        # rename extracted chroms dir -> Genomes/hg38
        genome_dir = rename(chromsdir, genome_dir)
        assert os.path.isdir(genome_dir)  # was asserting the pre-rename path
    else:
        chromgz = download(dest, http_url=f"{HG38URL}/chromosomes/{chrom}.fa.gz")
        os.makedirs(genome_dir, exist_ok=True)  # create Genomes/hg38_{chrom}
        chromfa = gunzip(
            chromgz,
            os.path.join(genome_dir, os.path.splitext(os.path.basename(chromgz))[0]),
        )  # decompress chrom FASTA
        assert os.path.isfile(chromfa)
    return os.path.abspath(genome_dir)


def ensure_vcf_dataset_directory(dest: str, dataset: str) -> str:
    """
    Ensure the existence of a directory for a specific VCF dataset within the
    specified destination directory.

    Args:
        dest (str): The destination directory where the VCF dataset directory
            should be created.
        dataset (str): The name or identifier of the VCF dataset.

    Returns:
        str: The path to the VCF dataset directory.
    """

    vcf_dataset_dir = os.path.join(dest, f"hg38_{dataset}")
    if not os.path.exists(vcf_dataset_dir):
        os.mkdir(vcf_dataset_dir)
    return vcf_dataset_dir


def download_vcf_data(chrom: str, dest: str, dataset: str) -> None:
    """
    Download VCF data for a specific chromosome and variant dataset to the destination
    directory.

    Args:
        chrom (str): The chromosome identifier in UCSC format.
        dest (str): The destination directory to save the downloaded VCF data.
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP", or 1000G+HGDP).

    Returns:
        None

    Raises:
        ValueError: If the input chromosome or variant dataset is invalid.
        FileExistsError: If the destination directory does not exist.
    """

    # assume chromosomes given in UCSC format (chr1, chr2, etc.)
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Forbidden input chromosome ({chrom})")
    if not os.path.isdir(dest):  # check dest directory existence
        raise FileExistsError(f"Unable to locate {dest}")
    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # create VCF dataset directory within VCFs folder
    sys.stderr.write(f"Downloading VCF data for chromsome(s) {chrom}\n")
    for ds in dataset.split("+"):
        vcf_dataset_dir = ensure_vcf_dataset_directory(dest, ds)
        ftp_server = VCF1000GSERVER if ds == "1000G" else VCFHGDPSERVER
        vcf_url = VCF1000GURL if ds == "1000G" else VCFHGDPURL
        chroms = CHROMS if chrom == "all" else [chrom]
        md5data = MD51000G if ds == "1000G" else MD5HGDP
        for c in chroms:
            expected_fname = os.path.basename(vcf_url.format(c))
            dest_path = os.path.join(vcf_dataset_dir, expected_fname)
            # Resume: skip files already downloaded and checksum-verified in a
            # previous (interrupted) run, so a partially completed multi-hour
            # download can be restarted cheaply instead of starting over.
            if os.path.isfile(dest_path) and md5data.get(
                expected_fname
            ) == compute_md5(dest_path):
                sys.stderr.write(
                    f"{expected_fname} already present and verified; skipping\n"
                )
                continue
            # Retry on checksum mismatch (e.g. a truncated transfer from a slow or
            # flaky host) instead of aborting the whole run on the first failure.
            attempts = 3
            for attempt in range(1, attempts + 1):
                if ds == "1000G":
                    # the EBI 1000G host also serves HTTPS; prefer it, since FTP is
                    # frequently blocked on CI runners and institutional/cloud networks
                    vcf = download(
                        vcf_dataset_dir,
                        http_url=f"https://{ftp_server}{vcf_url.format(c)}",
                    )
                else:  # HGDP (Sanger) via FTP
                    vcf = download(
                        vcf_dataset_dir,
                        ftp_conn=True,
                        ftp_server=ftp_server,
                        ftp_path=vcf_url.format(c),
                    )
                if md5data[os.path.basename(vcf)] == compute_md5(vcf):
                    break  # download verified
                sys.stderr.write(
                    f"Checksum mismatch for {os.path.basename(vcf)} "
                    f"(attempt {attempt}/{attempts}); re-downloading\n"
                )
                if attempt == attempts:
                    raise ValueError(
                        f"Download for {os.path.basename(vcf)} failed after "
                        f"{attempts} attempts"
                    )


def ensure_samplesids_directory(dest: str) -> str:
    """
    Ensure the existence of the 'samplesIDs' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'samplesIDs' directory
            should be created.

    Returns:
        str: The path to the 'samplesIDs' directory.
    """

    samplesids_dir = os.path.join(dest, CRISPRME_DIRS[6])
    if not os.path.exists(samplesids_dir):
        os.mkdir(samplesids_dir)
    return samplesids_dir


def download_samples_ids_data(dataset: str) -> None:
    """
    Download samples IDs data for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g.,
            "1000G", "HGDP").

    Returns:
        None

    Raises:
        ValueError: If the variant dataset is unknown.
    """

    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # samples ids folder must be located within current directory
    # -- see check_crisprme_directory_tree() for details
    sys.stderr.write(f"Downloading sample ids for dataset(s) {dataset}\n")
    samplesids_dir = ensure_samplesids_directory(os.getcwd())
    for ds in dataset.split("+"):
        samplesid_fname = (
            "samplesIDs.1000G.txt" if ds == "1000G" else "samplesIDs.HGDP.txt"
        )
        samplesids = download(
            samplesids_dir, http_url=f"{TESTDATAURL}/samplesIDs/{samplesid_fname}"
        )
        if MD5SAMPLES[os.path.basename(samplesids)] != compute_md5(samplesids):
            raise ValueError(f"Download for {os.path.basename(samplesids)} failed")
        rename(samplesids, os.path.join(samplesids_dir, f"hg38_{ds}.samplesID.txt"))

def ensure_annotation_directory(dest: str) -> str:
    """
    Ensure the existence of the 'annotation' directory within the specified
    destination directory.

    Args:
        dest (str): The destination directory where the 'annotation' directory
            should be created.

    Returns:
        str: The path to the 'annotation' directory.
    """

    annotation_dir = os.path.join(dest, CRISPRME_DIRS[4])
    if not os.path.exists(annotation_dir):
        os.mkdir(annotation_dir)
    return annotation_dir


def download_annotation_data() -> Tuple[str, str]:
    """
    Download gencode and encode annotation data to the 'annotation' directory
    within the current working directory.

    Returns:
        Tuple[str, str]: Paths to the downloaded gencode and encode annotation
            files.
    """

    sys.stderr.write("Downloading ENCODE and GENCODE annotation data\n")
    annotation_dir = ensure_annotation_directory(os.getcwd())
    gencode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/gencode.protein_coding.bed.tar.gz",
        "gencode.protein_coding.bed",
    )
    encode = _retrieve_ann_data(
        annotation_dir,
        "Annotations/dhs+encode+gencode.hg38.bed.tar.gz",
        "dhs+encode+gencode.hg38.bed",
    )
    return gencode, encode


def _retrieve_ann_data(annotation_dir: str, url: str, fname: str) -> str:
    """
    Download and extract an annotation file, then compress it with bgzip.

    This function downloads an annotation archive, verifies its integrity, extracts
    the specified file, and compresses it using bgzip. It returns the path to the
    compressed annotation file.

    Args:
        annotation_dir (str): The directory to store the annotation data.
        url (str): The URL of the annotation archive to download.
        fname (str): The name of the file to extract and compress.

    Returns:
        str: The path to the compressed annotation file.

    Raises:
        ValueError: If the downloaded file fails the MD5 check.
    """

    # download gencode annotation
    annfile_tar = download(annotation_dir, http_url=os.path.join(TESTDATAURL, url))
    if MD5ANNOTATION[os.path.basename(annfile_tar)] != compute_md5(annfile_tar):
        raise ValueError(f"Download for {os.path.basename(annfile_tar)} failed")
    return os.path.join(untar(annfile_tar, annotation_dir), fname)


def ensure_pams_directory(dest: str) -> str:
    """
    Ensure the existence of the 'PAMs' directory within the specified destination
    directory.

    Args:
        dest (str): The destination directory where the 'PAMs' directory should
            be created.

    Returns:
        str: The path to the 'PAMs' directory.
    """

    pams_dir = os.path.join(dest, CRISPRME_DIRS[5])
    if not os.path.exists(pams_dir):
        os.mkdir(pams_dir)
    return pams_dir


def write_ngg_pamfile() -> str:
    """
    Write a test PAM file containing the NGG sequence to the 'PAMs' directory
    within the current working directory.

    Returns:
        str: The path to the created test PAM file.
    """

    sys.stderr.write("Creating PAM file\n")
    pams_dir = ensure_pams_directory(
        os.getcwd()
    )  # PAMs directory must be in current working dir
    pamfile = os.path.join(pams_dir, "20bp-NGG-SpCas9.txt")
    try:
        with open(pamfile, mode="w") as outfile:
            outfile.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")  # 20 + 3 bp (NGG)
    except IOError as e:
        raise IOError("An error occurred while writing the test PAM file") from e
    return pamfile


def write_sg1617_guidefile() -> str:
    """
    Write a test guide file containing the sg1617 guide sequence.

    Returns:
        str: The path to the created test guide file.
    """

    sys.stderr.write("Creating guide file\n")
    guidefile = "sg1617_test_guide.txt"
    try:
        with open(guidefile, mode="w") as outfile:
            outfile.write("CTAACAGTTGCTTTTATCACNNN\n")  # sg1617 guide
    except IOError as e:
        raise IOError("An error occerred while writing the test guide file") from e
    return guidefile


def load_benchmarks() -> dict:
    """Load the benchmark registry (falls back to the canonical Cas9 case)."""
    try:
        with open(BENCHMARKS_JSON) as fin:
            return json.load(fin)
    except (OSError, ValueError):
        return {
            "thresholds": {"mm": 4, "bDNA": 1, "bRNA": 1},
            "benchmarks": [{
                "name": "cas9_sg1617", "nuclease": "SpCas9",
                "pam_name": "20bp-NGG-SpCas9.txt",
                "pam_content": "NNNNNNNNNNNNNNNNNNNNNGG 3",
                "guide_file": "sg1617_test_guide.txt",
                "guide_crisprme": "CTAACAGTTGCTTTTATCACNNN",
            }],
        }


def write_pamfile(pam_name: str, pam_content: str) -> str:
    """Write a benchmark PAM file into the working-dir PAMs directory."""
    sys.stderr.write(f"Creating PAM file {pam_name}\n")
    pamfile = os.path.join(ensure_pams_directory(os.getcwd()), pam_name)
    with open(pamfile, mode="w") as outfile:
        outfile.write(pam_content.rstrip("\n") + "\n")
    return pamfile


def write_guidefile(guide_file: str, guide_seq: str) -> str:
    """Write a benchmark guide file into the working directory."""
    sys.stderr.write(f"Creating guide file {guide_file}\n")
    with open(guide_file, mode="w") as outfile:
        outfile.write(guide_seq.rstrip("\n") + "\n")
    return guide_file


def write_vcf_config(dataset: str) -> str:
    """
    Write a test VCF list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test VCF list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # config vcf list file
    sys.stderr.write(f"Creating VCF config file for dataset(s) {dataset}\n")
    vcf_config = "vcf.config.test.txt"
    try:
        with open(vcf_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                outfile.write(f"hg38_{ds}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return vcf_config


def write_samplesids_config(dataset: str) -> str:
    """
    Write a test samples ID list file for a specific variant dataset.

    Args:
        dataset (str): The name or identifier of the variant dataset (e.g., "1000G",
            "HGDP").

    Returns:
        str: The path to the created test samples ID list file.
    """

    # support for 1000 GP and HGDP datasets
    if dataset not in ["1000G", "HGDP", "1000G+HGDP"]:
        raise ValueError(f"Unknown variant dataset ({dataset})")
    # configure sample ids list
    sys.stderr.write(f"Creating samples config file for dataset(s) {dataset}\n")
    samples_config = "samplesIDs.config.test.txt"
    fname_map = {"1000G": "hg38_1000G.samplesID.txt", "HGDP": "hg38_HGDP.samplesID.txt"}
    try:
        with open(samples_config, mode="w") as outfile:
            for ds in dataset.split("+"):
                outfile.write(f"{fname_map[ds]}\n")
    except IOError as e:
        raise IOError("An error occurred while writing the test VCF list") from e
    return samples_config


def run_crisprme_test(chrom: str, dataset: str, threads: int, debug: bool) -> None:
    """Execute the CRISPRme test workflow for a specified chromosome and dataset.

    This function orchestrates the downloading of necessary genomic and VCF data,
    prepares input files, and runs the CRISPRme command-line tool to perform a
    complete search.

    Args:
        chrom (str): The chromosome to be analyzed.
        dataset (str): The dataset identifier for VCF data.
        threads (int): The number of threads to use for processing.
        debug (bool): A flag indicating whether to run in debug mode.

    Raises:
        Any exceptions raised by the called functions or subprocess.
    """

    check_crisprme_directory_tree(os.getcwd())  # check crisprme directory tree
    check_output()  # check complete-test output folder
    genome_dir = download_genome_data(chrom, CRISPRME_DIRS[0])  # download genome data
    download_vcf_data(chrom, CRISPRME_DIRS[3], dataset)  # download vcf data
    vcf = write_vcf_config(dataset)  # write test vcf list
    download_samples_ids_data(dataset)  # download vcf dataset samples ids
    samplesids = write_samplesids_config(dataset)  # write test samples ids list
    # download gencode and encode annotation data
    gencode, encode = download_annotation_data()
    debug_arg = "--debug" if debug else ""
    # Run one complete-search per registered benchmark. complete-search refuses
    # to run into a non-empty output folder, so each benchmark gets its OWN
    # output dir (crisprme-test-out_<name>); validate-test looks in each.
    registry = load_benchmarks()
    th = registry["thresholds"]
    bmax = max(th["bDNA"], th["bRNA"])
    for bench in registry["benchmarks"]:
        output_dir = f"{COMPLETETESTRESDIR}_{bench['name']}"
        pam = write_pamfile(bench["pam_name"], bench["pam_content"])
        guide = write_guidefile(bench["guide_file"], bench["guide_crisprme"])
        sys.stderr.write(
            f"Running complete-search for benchmark '{bench['name']}' "
            f"({bench.get('nuclease', '')}) -> {output_dir}\n"
        )
        crisprme_cmd = (
            f"crisprme.py complete-search --genome {genome_dir} "
            f"--bmax {bmax} --mm {th['mm']} --bDNA {th['bDNA']} --bRNA {th['bRNA']} "
            f"--merge 3 --pam {pam} --guide {guide} --vcf {vcf} "
            f"--samplesID {samplesids} --annotation {encode} "
            f"--gene_annotation {gencode} --output {output_dir} "
            f"--thread {threads} {debug_arg} --ci-cd-test"
        )
        returncode = subprocess.call(crisprme_cmd, shell=True)
        if returncode != 0:
            sys.stderr.write(
                "ERROR: complete-test failed during complete-search for benchmark "
                f"'{bench['name']}' (exit code {returncode}). See the log above.\n"
            )
            sys.exit(returncode)


def main():
    chrom, dataset, threads, debug = sys.argv[1:]  # read commandline args
    debug = debug == "True"
    run_crisprme_test(chrom, dataset, int(threads), debug)  # run crisprme test


if __name__ == "__main__":
    main()
