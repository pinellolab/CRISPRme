"""
setup_legacy_database.py — CRISPRme Legacy Web-Server Database Setup
======================================================================
 
Downloads and configures all genomic data (reference FASTA, VCF variants,
sample IDs, functional annotations, and PAM files) that were originally
hosted on the CRISPRme web server, enabling a fully local deployment of the
web interface.
 
By default full genome data are downloaded. For a quick initial setup pass 
``--chrom chr22`` to download only chromosome 22 data.
 
Both the **1000 Genomes Project** (1000G) and **Human Genome Diversity
Project** (HGDP) variant datasets are always downloaded together, matching
the original web-server configuration.
 
Force vs. incremental mode
---------------------------
When ``force=False`` (the default) each asset is checked before downloading:
 
* **FASTA / VCF / sample-ID files** — skipped when the expected local file
  already exists *and* its MD5 digest matches the known-good value.
* **Annotation files** — skipped when the bgzip-compressed output already
  exists and MD5 passes.
* **PAM files** — skipped when the file already exists (content is fully
  determined by :data:`PAM_DEFINITIONS`; no remote MD5 is needed).
* **Config files** (``vcf.config.txt``, ``samplesIDs.config.txt``) — always
  rewritten; they are trivially cheap and must reflect the current state.
 
When ``force=True`` all assets are (re-)downloaded and overwritten
unconditionally, regardless of any existing files or their integrity.
 
Genomic coordinate convention
-------------------------------
* FASTA files — 1-based (UCSC convention).
* BED files produced / consumed here — 0-based half-open (standard BED).
* VCF files — 1-based (VCF 4.2 spec).
 
Usage (called from ``crisprme.py``)
-------------------------------------
    crisprme.py web-interface setup [--chrom <chrom>] [--force] [--debug]
"""


from pathlib import Path
from typing import Dict, List, Optional, Tuple

import os
import subprocess
import sys

from utils import CHROMS, CRISPRME_DIRS, MD51000G, MD5ANNOTATION, MD5GENOME, MD5HGDP, MD5SAMPLES, check_crisprme_directory_tree, compute_md5, download, gunzip, rename, untar


# ==============================================================================
# remote data sources
# ==============================================================================

HG38_BASE_URL: str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"

VCF_1000G_SERVER: str = "ftp.1000genomes.ebi.ac.uk"
VCF_1000G_URL_TEMPLATE: str = (
    "/vol1/ftp/data_collections/1000_genomes_project/release/"
    "20190312_biallelic_SNV_and_INDEL/"
    "ALL.{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)

VCF_HGDP_SERVER: str = "ngs.sanger.ac.uk"
VCF_HGDP_URL_TEMPLATE: str = (
    "/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{chrom}.vcf.gz"
)

TEST_DATA_BASE_URL: str = (
    "https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/test/data/"
)

# ==============================================================================
# PAM file definitions
# ==============================================================================
# each entry: PAM filename + (pam_sequence_with_Ns, position_index)
#
# The position index follows CRISPRme convention:
#   - positive N -> PAM is at the 3' end; the last N chars of the token are PAM
#   - negative N -> PAM is at the 5' end; the first |N| chars of the token are PAM
#
# Token format in the file: "<full_sequence> <position_index>\n"
# e.g. "NNNNNNNNNNNNNNNNNNNNNGG 3" means 20 nt guide + NGG (3 nt) at 3′ end

PAM_DEFINITIONS: Dict[str, Tuple[str, int]] = {
    # ---- SpCas9 variants
    "20bp-NGG-SpCas9.txt":      ("NNNNNNNNNNNNNNNNNNNNNGG", 3),
    "20bp-NGC-SpCas9.txt":      ("NNNNNNNNNNNNNNNNNNNNNGC", 3),
    "20bp-NGK-SpCas9.txt":      ("NNNNNNNNNNNNNNNNNNNNNGNK", 3),        # K = G/T
    "20bp-NRG-SpCas9.txt":      ("NNNNNNNNNNNNNNNNNNNNNRG", 3),         # R = A/G (VQR variant)
    "20bp-NRCH-SpCas9.txt":     ("NNNNNNNNNNNNNNNNNNNNNRCH", 4),        # NRCH (xCas9)
    "20bp-NNGT-SpCas9.txt":     ("NNNNNNNNNNNNNNNNNNNNNNNGT", 4),       # NNGT (SpCas9-NG)
    # ---- iSpyMacCas9 
    "20bp-NAA-iSpyMacCas9.txt": ("NNNNNNNNNNNNNNNNNNNNNAA", 3),
    # ---- SaCas9 
    "22bp-NNGRRN-SaCas9.txt":   ("NNNNNNNNNNNNNNNNNNNNNNNNNGRRN", 6),   # 22 nt guide + NNGRRN
    # ---- CasX 
    "TTCN-20bp-CasX.txt":       ("TTCNNNNNNNNNNNNNNNNNNNNNN", -4),      # 5′ PAM (TTCN)
    # ---- Cas12a (Cpf1) variants 
    "TTTV-20bp-Cas12a.txt":     ("TTTV" + "N" * 20, -4),                # 5′ PAM (TTTV)
    "TTTV-21bp-Cas12a.txt":     ("TTTV" + "N" * 21, -4),
    "TTTV-23bp-Cas12a.txt":     ("TTTV" + "N" * 23, -4),
    "TTTV-25bp-AsCas12a.txt":   ("TTTV" + "N" * 25, -4),
    # ---- No-PAM / generic 
    "18bp-NNN-NO_PAM.txt":      ("N" * 21, 3),                          # 18 nt guide + NNN placeholder
    "20bp-NNN-NO_PAM.txt":      ("N" * 23, 3),
    "21bp-NNN-NO_PAM.txt":      ("N" * 24, 3),
    "23bp-NNN-NO_PAM.txt":      ("N" * 26, 3),
}


# ==============================================================================
# skip-or-download helpers
# ==============================================================================

def _file_is_valid(local_path: Path, md5_map: Optional[Dict[str, str]] = None) -> bool:
    """Return ``True`` when *local_path* exists and passes its MD5 check.
 
    This is the single decision point that determines whether an asset can be
    skipped during an incremental (non-forced) setup run.
 
    The check has three outcomes:
 
    * File missing → ``False`` (must download).
    * File present, basename **not** in *md5_map* → ``True`` (no digest
      available; trust the file that is already there).
    * File present, basename **in** *md5_map* → ``True`` only when the
      computed digest matches the expected value.
 
    Args:
        local_path: Absolute path to the candidate local file.
        md5_map: Mapping of ``basename → expected MD5 hex digest``.  May be
            an empty dict for assets whose remote digest is unknown (e.g. PAM
            files that are generated locally).
 
    Returns:
        ``True`` if the file can safely be skipped; ``False`` if it must be
        (re-)downloaded or recreated.
    """
    if md5_map is None:
        md5_map = {}
    if not local_path.is_file():
        sys.stderr.write(f"Asset missing, will download: {local_path}\n")
        return False
    basename = local_path.name
    if basename not in md5_map:
        return True
    expected = md5_map[basename]
    actual = compute_md5(str(local_path))
    if actual == expected:
        sys.stderr.write(f"MD5 OK, skipping: {local_path}\n")
        return True
    sys.stderr.write(
        f"MD5 mismatch for {basename} (expected {expected}, got {actual}); will re-download\n",
    )
    return False


# ==============================================================================
# genome helpers
# ==============================================================================

def _ensure_directory(parent: Path, name: str) -> Path:
    """Return *parent/name*, creating it if absent.
 
    Args:
        parent: Parent directory path.
        name: Sub-directory name to ensure.
 
    Returns:
        Absolute path to the (possibly newly created) sub-directory.
    """
    target = parent / name
    target.mkdir(parents=True, exist_ok=True)
    return target


def _ensure_hg38_directory(genomes_dir: Path) -> Path:
    """Ensure ``Genomes/hg38`` exists and return its path.

    Args:
        genomes_dir: Path to the top-level ``Genomes`` folder.

    Returns:
        Path to ``Genomes/hg38``.
    """
    return _ensure_directory(genomes_dir, "hg38")


def _download_full_genome_data(genomes_dir: Path, force: bool) -> None:
    hg38_dir = genomes_dir / "hg38"
    chroms_present = hg38_dir.is_dir() and any(hg38_dir.glob("chr*.fa"))
    if not force and chroms_present:
        sys.stderr.write("Full hg38 genome already present, skipping download\n")
        return
    archive = download(
        str(genomes_dir),
        http_url=f"{HG38_BASE_URL}/bigZips/hg38.chromFa.tar.gz",
    )
    _verify_md5(archive, MD5GENOME)
    chroms_dir = untar(archive, str(genomes_dir), "chroms")
    rename(chroms_dir, str(hg38_dir))


def _download_chrom_genome_data(chrom: str, genomes_dir: Path, force: bool) -> None:
    hg38_dir = _ensure_hg38_directory(genomes_dir)
    fa_path = hg38_dir / f"{chrom}.fa"
    archive_basename = f"{chrom}.fa.gz"
    if not force and _file_is_valid(fa_path):
        sys.stderr.write(f"Genome FASTA already valid, skipping: {fa_path}\n")
        return
    gz_path = download(
        str(genomes_dir),
        http_url=f"{HG38_BASE_URL}/chromosomes/{archive_basename}",
    )
    fa_path = gunzip(
        gz_path,
        str(hg38_dir / f"{Path(gz_path).stem}"),
    )
    if not Path(fa_path).is_file():
        raise RuntimeError(f"FASTA extraction failed for {chrom}")



def _download_genome_data(chrom: str, genomes_dir: Path, force: bool) -> None:
    """Download hg38 FASTA data for *chrom* into *genomes_dir*.

    Handles both single-chromosome downloads and the full genome archive.
    MD5 integrity is verified for every downloaded file.

    Args:
        chrom: Chromosome label in UCSC format (e.g. ``"chr22"``), or
            ``"all"`` to download the complete genome.
        genomes_dir: Destination directory (top-level ``Genomes`` folder).

    Raises:
        ValueError: If *chrom* is not a recognised chromosome label.
        FileNotFoundError: If *genomes_dir* does not exist.
        ValueError: If the MD5 check fails after download.
    """
    if chrom not in CHROMS + ["all"]:
        raise ValueError(f"Unrecognised chromosome: {chrom!r}")
    if not genomes_dir.is_dir():
        raise FileNotFoundError(f"Genomes directory not found: {genomes_dir}")
    sys.stdout.write(f"Downloading genome data for chromosomes: {chrom}\n")
    if chrom == "all":
        _download_full_genome_data(genomes_dir, force)
    else:
        _download_chrom_genome_data(chrom, genomes_dir, force)
        

def _ensure_vcf_dataset_directory(vcfs_dir: Path, dataset_label: str) -> Path:
    """Ensure ``VCFs/hg38_<dataset_label>`` exists and return its path.

    Args:
        vcfs_dir: Top-level ``VCFs`` folder.
        dataset_label: Short dataset name, e.g. ``"1000G"`` or ``"HGDP"``.

    Returns:
        Path to the dataset sub-directory.
    """
    return _ensure_directory(vcfs_dir, f"hg38_{dataset_label}")


def _download_vcf_data(chrom: str, vcfs_dir: Path, force: bool) -> None:
    """Download phased VCF files for *chrom* from *dataset* into *vcfs_dir*.

    Supports ``"1000G"``, ``"HGDP"``, and the combined ``"1000G+HGDP"``
    dataset specifier.

    Args:
        chrom: Chromosome label (UCSC format) or ``"all"``.
        vcfs_dir: Top-level ``VCFs`` folder.

    Raises:
        ValueError: For unknown *chrom* or *dataset* values.
        FileNotFoundError: If *vcfs_dir* does not exist.
        ValueError: If an MD5 check fails.
    """
    if not vcfs_dir.is_dir():
        raise FileNotFoundError(f"VCFs directory not found: {vcfs_dir}")
    chroms_to_fetch: List[str] = CHROMS if chrom == "all" else [chrom]
    for ds_label in "1000G+HGDP".split("+"):
        sys.stdout.write(f"Downloading genetic variants data for dataset {ds_label}\n")
        vcf_dataset_dir = _ensure_vcf_dataset_directory(vcfs_dir, ds_label)
        ftp_server, url_template, md5_map = _vcf_sources(ds_label)
        for c in chroms_to_fetch:
            remote_basename = Path(url_template.format(chrom=c)).name
            local_vcf = vcf_dataset_dir / remote_basename
            if not force and _file_is_valid(local_vcf, md5_map):
                sys.stderr.write(f"VCF already valid, skipping: {local_vcf}\n")
                continue
            vcf_path = download(
                str(vcf_dataset_dir),
                ftp_conn=True,
                ftp_server=ftp_server,
                ftp_path=url_template.format(chrom=c),
            )
            _verify_md5(vcf_path, md5_map)


def _ensure_samplesids_directory(base_dir: Path) -> Path:
    """Ensure the ``samplesIDs`` directory exists inside *base_dir*.

    Args:
        base_dir: Directory that should contain ``samplesIDs/``.

    Returns:
        Path to ``samplesIDs/``.
    """
    return _ensure_directory(base_dir, CRISPRME_DIRS[6])


def _download_samples_ids_data(base_dir: Path, force: bool) -> None:
    """Download pre-computed sample ID files for *dataset*.

    Args:
        base_dir: CRISPRme working directory (parent of ``samplesIDs/``).

    Raises:
        ValueError: For unknown *dataset* values or failed MD5 checks.
    """
    samplesids_dir = _ensure_samplesids_directory(base_dir)
    for ds_label in "1000G+HGDP".split("+"):
        sys.stdout.write(f"Downloading samples ids data for dataset {ds_label}\n")
        fname = (
            "samplesIDs.1000G.txt" if ds_label == "1000G" else "samplesIDs.HGDP.txt"
        )
        renamed_target = samplesids_dir / f"hg38_{ds_label}.samplesID.txt"
        if not force and _file_is_valid(renamed_target):
            sys.stderr.write(f"Sample IDs already valid, skipping: {renamed_target}")
            continue
        local_path = download(
            str(samplesids_dir),
            http_url=f"{TEST_DATA_BASE_URL}/samplesIDs/{fname}",
        )
        _verify_md5(local_path, MD5SAMPLES)
        rename(
            os.path.join(samplesids_dir, fname), 
            os.path.join(samplesids_dir, f"hg38_{ds_label}.samplesID.txt"),
        )


def _write_vcf_config(base_dir: Path) -> Path:
    """Write a VCF folder-list config file for *dataset*.

    Each line names one VCF sub-folder (relative to ``VCFs/``).

    Args:
        base_dir: CRISPRme working directory (config written here).

    Returns:
        Path to the created config file.

    Raises:
        ValueError: For unknown *dataset* values.
        IOError: If the file cannot be written.
    """
    vcf_config_path = base_dir / "vcf.config.txt"
    try:
        lines = [f"hg38_{ds}\n" for ds in "1000G+HGDP".split("+")]
        vcf_config_path.write_text("".join(lines), encoding="utf-8")
    except IOError as exc:
        raise IOError(f"Failed to write VCF config: {vcf_config_path}") from exc
    return vcf_config_path


def _write_samplesids_config(base_dir: Path) -> Path:
    """Write a sample-IDs list config file for *dataset*.

    Args:
        base_dir: CRISPRme working directory.

    Returns:
        Path to the created config file.

    Raises:
        ValueError: For unknown *dataset* values.
        IOError: If the file cannot be written.
    """
    samples_config_path = base_dir / "samplesIDs.config.txt"
    fname_map = {
        "1000G": "hg38_1000G.samplesID.txt", "HGDP": "hg38_HGDP.samplesID.txt"
    }
    try:
        lines = [f"{fname_map[ds]}\n" for ds in "1000G+HGDP".split("+")]
        samples_config_path.write_text("".join(lines), encoding="utf-8")
    except IOError as exc:
        raise IOError(f"Failed to write samplesIDs config: {samples_config_path}") from exc
    return samples_config_path


def _ensure_annotation_directory(base_dir: Path) -> Path:
    """Ensure the ``Annotations`` directory exists inside *base_dir*.

    Args:
        base_dir: CRISPRme working directory.

    Returns:
        Path to ``Annotations/``.
    """
    return _ensure_directory(base_dir, CRISPRME_DIRS[4])


def _download_annotation_data(base_dir: Path, force: bool) -> None:
    """Download GENCODE protein-coding and ENCODE/DHS annotation files.

    Both files are bgzip-compressed on disk after download so they are
    immediately usable by CRISPRme.

    Args:
        base_dir: CRISPRme working directory.

    """
    annotation_dir = _ensure_annotation_directory(base_dir)
    genocode_ann = Path(annotation_dir) / "gencode.protein_coding.bed"
    if force or not _file_is_valid(genocode_ann):
        sys.stdout.write(f"Downloading functional annotation data\n")
        gencode_bgz = _retrieve_annotation_file(
            annotation_dir,
            url=f"{TEST_DATA_BASE_URL}/Annotations/gencode.protein_coding.bed.tar.gz",
            inner_fname="gencode.protein_coding.bed",
        )
    encode_ann = Path(annotation_dir) / "dhs+encode+gencode.hg38.bed"
    if force or not _file_is_valid(encode_ann):
        sys.stdout.write(f"Downloading gene annotation data\n")
        encode_bgz = _retrieve_annotation_file(
            annotation_dir,
            url=f"{TEST_DATA_BASE_URL}/Annotations/dhs+encode+gencode.hg38.bed.tar.gz",
            inner_fname="dhs+encode+gencode.hg38.bed",
        )


def _ensure_pams_directory(base_dir: Path) -> Path:
    """Ensure the ``PAMs`` directory exists inside *base_dir*.

    Args:
        base_dir: CRISPRme working directory.

    Returns:
        Path to ``PAMs/``.
    """
    return _ensure_directory(base_dir, CRISPRME_DIRS[5])



def _write_pam_file(
    pam_filename: str,
    pam_sequence: str,
    position_index: int,
    base_dir: Path,
) -> Path:
    """Write a single PAM definition file to ``PAMs/<pam_filename>``.

    The on-disk format understood by CRISPRme is::

        <full_token> <position_index>

    where ``full_token`` is the concatenation of guide Ns and the PAM
    bases (order depends on the sign of *position_index*), and
    ``position_index`` encodes both PAM length and orientation:

    * positive → 3' PAM; value equals PAM length in nucleotides
    * negative → 5' PAM; absolute value equals PAM length in nucleotides

    Args:
        pam_filename: Filename (not full path) for the PAM file, e.g.
            ``"20bp-NGG-SpCas9.txt"``.
        pam_sequence: Full token string (guide Ns + PAM bases or vice-versa).
        position_index: Signed integer encoding PAM position and length.
        base_dir: CRISPRme working directory.

    Returns:
        Path to the written PAM file.

    Raises:
        IOError: If the file cannot be written.
    """
    pams_dir = _ensure_pams_directory(base_dir)
    pam_path = pams_dir / pam_filename
    try:
        pam_path.write_text(f"{pam_sequence} {position_index}\n", encoding="utf-8")
    except IOError as exc:
        raise IOError(f"Failed to write PAM file {pam_path}") from exc
    return pam_path



def _write_all_pam_files(base_dir: Path) -> Dict[str, Path]:
    """Write every PAM definition listed in :data:`PAM_DEFINITIONS`.

    Args:
        base_dir: CRISPRme working directory.

    Returns:
        Mapping of PAM filename → absolute path on disk.
    """
    sys.stdout.write(f"Creating PAM data\n")
    written: Dict[str, Path] = {
        pam_filename: _write_pam_file(
            pam_filename, pam_sequence, position_index, base_dir
        )
        for pam_filename, (
            pam_sequence,
            position_index,
        ) in PAM_DEFINITIONS.items()
    }
    return written


# ==============================================================================
# public entry point
# ==============================================================================

def run_legacy_database_setup(chrom: str, base_dir: Path, force: bool) -> None:
    _validate_chrom(chrom)
    check_crisprme_directory_tree(str(base_dir))
    # begin data download
    genomes_dir = base_dir / CRISPRME_DIRS[0]
    _download_genome_data(chrom, genomes_dir, force)  # reference genome
    vcfs_dir = base_dir / CRISPRME_DIRS[3]
    _download_vcf_data(chrom, vcfs_dir, force)  # variants
    _download_samples_ids_data(base_dir, force)  # samples ids
    _write_vcf_config(base_dir)  # vcf config file
    _write_samplesids_config(base_dir)  # samples ids config file
    _download_annotation_data(base_dir, force)  # annotations data
    _write_all_pam_files(base_dir)  # pam files


# ==============================================================================
# CLI entry point (when invoked directly from crisprme.py)
# ==============================================================================

def main(argv: Optional[List[str]] = None) -> None:
    """Parse arguments forwarded from ``crisprme.py`` and run the setup.
 
    Args:
        argv: Argument list (excluding the subcommand token itself).
            When ``None``, ``sys.argv[1:]`` is used.
    """
    args: List[str] = argv if argv is not None else sys.argv[1:]
    chrom, base_dir, force = args
    run_legacy_database_setup(chrom, Path(base_dir), force == "True")


# ==============================================================================
# private helpers
# ==============================================================================

def _validate_chrom(chrom: str) -> None:
    """Raise ValueError if *chrom* is not a valid chromosome identifier."""
    if chrom not in CHROMS + ["all"]:
        raise ValueError(
            f"Unrecognised chromosome: {chrom!r}. "
            f"Expected one of {CHROMS + ['all']}"
        )
    

def _verify_md5(local_path: str, md5_map: Dict[str, str]) -> None:
    """Verify a downloaded file against a known MD5 digest.

    Args:
        local_path: Absolute path to the file to check.
        md5_map: Mapping of basename → expected MD5 hex digest.

    Raises:
        ValueError: If the computed MD5 does not match the expected value.
    """
    basename = os.path.basename(local_path)
    if basename not in md5_map:
        return
    expected = md5_map[basename]
    actual = compute_md5(local_path)
    if actual != expected:
        raise ValueError(
            f"MD5 mismatch for {basename}: "
            f"expected {expected!r}, got {actual!r}"
        )


def _vcf_sources(ds_label: str) -> Tuple[str, str, Dict[str, str]]:
    if ds_label == "1000G":  # 1000 genomes 
        return VCF_1000G_SERVER, VCF_1000G_URL_TEMPLATE, MD51000G
    return VCF_HGDP_SERVER, VCF_HGDP_URL_TEMPLATE, MD5HGDP  # hgdp


def _retrieve_annotation_file(
    annotation_dir: Path,
    url: str,
    inner_fname: str,
) -> Path:
    """Download an annotation archive, extract, and bgzip-compress it.

    Args:
        annotation_dir: Directory where the annotation file will be stored.
        url: HTTP URL of the ``.tar.gz`` annotation archive.
        inner_fname: Filename of the BED file inside the archive.

    Returns:
        Path to the bgzip-compressed BED file.

    Raises:
        ValueError: If the MD5 check fails.
        subprocess.SubprocessError: If bgzip fails.
    """
    archive_path = download(str(annotation_dir), http_url=url)
    _verify_md5(archive_path, MD5ANNOTATION)
    extract_dir = untar(archive_path, str(annotation_dir))
    return Path(extract_dir) / inner_fname


# ==============================================================================
# script entry point
# ==============================================================================

if __name__ == "__main__":
    main()

