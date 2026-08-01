"""
This module implements CRISPRme's pre-flight input validator.

It performs a lightweight, always-on pass over the genome, VCF, PAM, guide, and
annotation inputs to `complete-search`, catching misconfigurations that would
otherwise only surface deep into a multi-hour run (or silently corrupt results
without any error at all). Every check here is grounded in a specific failure
mode traced through the installed CRISPRitz `enricher.py` / `crispritz.py`
source — see `validate_inputs_plan.md` at the repo root for the full mapping
from check to failure mode.

A slower, opt-in full-file scan (`run_full`, `--full_input_validate`) covers
issues that can only be caught by reading an entire VCF: structural variants
anywhere in the file, full chromosome coverage, INFO/AF field-position
consistency, FILTER PASS ratio, POS bounds, high-allele-count multiallelic
sites, breakend notation, duplicate records, mixed GT phasing, and indels near
a chromosome start.
"""

from typing import Dict, List, NamedTuple, Optional, Tuple

import gzip
import os
import sys

ERROR = "ERROR"
WARN = "WARN"

# valid IUPAC nucleotide codes (matches the pam_dict keys in crisprme.py's
# getGuides(), used to build all possible concrete PAM sequences)
IUPAC_BASES = set("ACGTRYSWKMBDHVN")
GZIP_MAGIC = b"\x1f\x8b"
# #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, + >=1 sample column
MIN_VCF_HEADER_FIELDS = 10
# #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
MIN_VCF_DATA_FIELDS = 8
# FILTER value the officially-released crispritz enricher.py hardcodes as
# "passing" (enricher.py:363: `if line[6] != 'PASS': continue`) — verified
# against GitHub pinellolab/CRISPRitz master. CRISPRme's own
# --vcf-filter-pass-values defaults to "PASS,." (intending to also accept
# '.', e.g. HPRC vcfwave output), but that flag isn't actually wired through
# to crispritz's add-variants call (see validate_inputs_plan.md's "Known
# separate bug" section), so today only literal 'PASS' passes regardless of
# the flag. (A local, not-yet-submitted CRISPRitz patch accepting '.' exists
# only in one developer's own conda env for HPRC testing — not part of any
# officially shipped crispritz release, so not assumed here.)
ENRICHER_PASS_VALUES = ("PASS",)
STRUCTURAL_VARIANT_LEN = 50
# enricher.py:302 builds indel context as genomeStr[pos-26 : pos+26+len(ref)];
# a negative start silently wraps via Python slicing instead of erroring
INDEL_START_PROXIMITY = 26
# enricher.py's allele-index matching does `str(value) in gt_string`, a
# substring test that misfires once two-digit allele indices appear
MULTIALLELIC_WARN_THRESHOLD = 10
# warn (not error) below this PASS ratio; a dataset legitimately restricted to
# one difficult chromosome can have a low ratio without being misconfigured
FILTER_PASS_WARN_RATIO = 0.10
# mock file CRISPRme substitutes for optional file arguments; skip content
# checks against it
MOCK_FILENAME = "vuoto.txt"


class Issue(NamedTuple):
    """A single validation finding.

    Attributes:
        severity: Either `ERROR` (aborts the run) or `WARN` (printed, non-fatal).
        message: Human-readable description of the finding.
    """

    severity: str
    message: str


class ValidationReport:
    """Accumulates validation issues and OK messages and renders a report.

    Check functions return plain lists of `Issue`; the report only exists to
    collect them alongside informational OK lines and produce the printed
    summary, so individual checks stay simple to unit test in isolation.
    """

    def __init__(self, title: str) -> None:
        self.title = title
        self._lines: List[str] = []
        self._error_count = 0
        self._warn_count = 0

    def ok(self, message: str) -> None:
        """Records an informational line for a check that passed cleanly."""
        self._lines.append(f"  [OK] {message}")

    def add(self, issues: List[Issue], ok_message: Optional[str] = None) -> None:
        """Records the outcome of a single check.

        Args:
            issues: Issues returned by a check function; empty if the check
                passed cleanly.
            ok_message: Message to print if `issues` is empty. Ignored if
                `issues` is non-empty.
        """
        if not issues:
            if ok_message:
                self.ok(ok_message)
            return
        for issue in issues:
            if issue.severity == ERROR:
                self._error_count += 1
            else:
                self._warn_count += 1
            self._lines.append(f"  [{issue.severity}] {issue.message}")

    def has_errors(self) -> bool:
        return self._error_count > 0

    def write(self) -> None:
        """Writes the rendered report to the appropriate stream.

        Goes to stderr only when the report contains an error (the run is
        about to abort); otherwise goes to stdout, since a clean or
        warning-only report doesn't block anything and shouldn't be mistaken
        for a failure signal by callers that treat stderr output as fatal
        (see issue #94).
        """
        stream = sys.stderr if self.has_errors() else sys.stdout
        stream.write(self.render())

    def render(self) -> str:
        """Renders the accumulated lines into the final report text."""
        divider = "=" * len(self.title)
        out = [self.title, divider, *self._lines, divider]
        if self.has_errors():
            out.append(
                f"Result: FAILED ({self._error_count} error"
                f"{'s' if self._error_count != 1 else ''}) — "
                "aborting before launching complete-search"
            )
        elif self._warn_count:
            out.append(
                f"Result: PASSED with {self._warn_count} warning"
                f"{'s' if self._warn_count != 1 else ''}"
            )
        else:
            out.append("Result: PASSED")
        return "\n".join(out) + "\n"


def _af_field(info_field: str) -> Tuple[Optional[int], Optional[str]]:
    """Locates the field `enricher.py` will treat as the site's AF.

    Mirrors `enricher.py:362` exactly (`if ele[0:2] == "AF"`) — a prefix
    match, not an exact-key match — so this returns the *same* field crispritz
    will actually use, including when it's the wrong one.

    Returns:
        `(position, key)` of the first INFO entry whose key starts with `AF`,
        or `(None, None)` if none is found.
    """
    for pos, entry in enumerate(info_field.split(";")):
        if entry[0:2] == "AF":
            return pos, entry.split("=", 1)[0]
    return None, None


def _af_value_parses_as_numeric(info_field: str, af_pos: int) -> bool:
    """Checks whether the AF field's value parses as one or more floats.

    Traced end-to-end (not assumed): `enricher.py:225` extracts this value
    with `line[7].split(";")[pos_AF][3:]` -- a hardcoded 3-character slice
    that only works correctly for an exact `AF=` key (see the ERROR case
    this check's caller already raises for a non-exact key). Even when the
    key *is* exactly `AF`, the value itself (e.g. `.`, empty, or otherwise
    non-numeric) is carried completely unvalidated: `enricher.py`'s dict
    entry -> `new_simple_analysis.py`'s `retrieveFromDict()`
    (`AF_list.append(split_entry[3].strip())`, no parsing) -> straight into
    the results file's `AF` column. `CRISPRme_plots.py` then does
    `pd.to_numeric(df["AF"])` with no `errors=` argument (defaults to
    `errors="raise"`) and no surrounding try/except -- a non-numeric value
    reaching that line crashes plot generation, near the very end of a
    multi-hour run.

    Args:
        info_field: The raw INFO column string.
        af_pos: Position of the AF entry within `info_field.split(";")`
            (from `_af_field`); only meaningful when that call's key was
            exactly `"AF"`.

    Returns:
        True if every comma-separated value (multi-allelic sites list one
        per ALT allele) parses as a float.
    """
    raw = info_field.split(";")[af_pos]
    if "=" not in raw:
        return False
    value_str = raw.split("=", 1)[1]
    if not value_str:
        return False
    for v in value_str.split(","):
        try:
            float(v)
        except ValueError:
            return False
    return True


def check_genome_fasta(genomedir: str) -> Tuple[List[Issue], List[str]]:
    """Checks that genome FASTA files are well-formed and collects chromosome
    identifiers.

    Chromosome identity for genome-VCF matching is filename-derived (matches
    how `submit_job_automated_new_multiple_vcfs.sh` and `crispritz.py`
    identify chromosomes when pairing enriched genomes with VCFs), not parsed
    from FASTA header content.

    Args:
        genomedir: Path to the reference genome directory (per-chromosome
            FASTA files).

    Returns:
        A tuple `(issues, chrom_names)` where `chrom_names` are FASTA
        filenames with their extension stripped (e.g. "chr1.fa" -> "chr1").
    """
    issues: List[Issue] = []
    chrom_names: List[str] = []
    if not os.path.isdir(genomedir):
        return [Issue(ERROR, f"Genome directory does not exist: {genomedir}")], []
    fasta_files = sorted(
        f for f in os.listdir(genomedir) if f.endswith((".fa", ".fasta"))
    )
    if not fasta_files:
        issues.append(
            Issue(ERROR, f"No .fa/.fasta files found in genome directory: {genomedir}")
        )
        return issues, chrom_names
    for fname in fasta_files:
        chrom_names.append(fname.rsplit(".", 1)[0])
        fpath = os.path.join(genomedir, fname)
        try:
            with open(fpath, "r") as fin:
                first_char = fin.read(1)
        except OSError as e:
            issues.append(Issue(ERROR, f"{fname}: could not read file ({e})"))
            continue
        if first_char != ">":
            issues.append(
                Issue(
                    ERROR,
                    f"{fname}: does not start with '>' — not a valid FASTA file",
                )
            )
    return issues, chrom_names


def check_vcf_filename_chrom(vcf_files: List[str]) -> Tuple[List[Issue], Dict[str, str]]:
    """Checks that each VCF filename contains a `chrN`-style token.

    Crispritz extracts chromosome identity by splitting the filename on `.`
    and looking for a token containing `"chr"` (`crispritz.py`
    `genomeEnrichment()`, lines ~609-617). A filename without this token
    silently produces an empty chromosome identifier, the dataset's genome
    path is built wrong, and the chromosome is dropped from the search with
    no error.

    Args:
        vcf_files: Absolute paths to `.vcf.gz` files in one dataset directory.

    Returns:
        A tuple `(issues, chrom_by_file)` mapping each filename (basename) to
        its extracted `chrN` token, or `""` if none was found.
    """
    issues: List[Issue] = []
    chrom_by_file: Dict[str, str] = {}
    for vf in vcf_files:
        fname = os.path.basename(vf)
        chrom = ""
        for token in fname.split("."):
            if "chr" in token:
                chrom = token
                break
        chrom_by_file[fname] = chrom
        if not chrom:
            issues.append(
                Issue(
                    ERROR,
                    f"{fname}: filename contains no chrN token — rename to e.g. "
                    "chr1.vcf.gz",
                )
            )
    return issues, chrom_by_file


def check_vcf_chrom_matches_genome(
    chrom_by_file: Dict[str, str], genome_chroms: List[str]
) -> List[Issue]:
    """Checks that each VCF's filename-derived chromosome exists in the genome
    directory.

    A mismatch (e.g. `1.vcf.gz` vs. a genome directory using `chr1.fa`) puts
    that chromosome in crispritz's "no vcf" bucket — it's silently treated as
    reference-only with no warning (`crispritz.py` `genomeEnrichment()`,
    `chr_wihtout_vcf` computation, lines ~625-634). This is the AoU
    `chr`-prefix bug in its general form.

    One aggregated `WARN` per dataset, not one `ERROR` per mismatched file
    (changed 2026-07-31): the historical bug this guards against was a
    systematic, whole-dataset naming mismatch, but a single legitimately-
    absent chromosome (e.g. no `chrY.fa` in the genome build, while the VCF
    dataset happens to include a `chrY.vcf.gz`) shouldn't block an entire
    multi-dataset run. Informing the user which chromosomes will be
    silently reference-only is the actual goal, not blocking on it.

    Args:
        chrom_by_file: Mapping produced by `check_vcf_filename_chrom`.
        genome_chroms: Chromosome identifiers from the genome directory
            (from `check_genome_fasta`).

    Returns:
        A single WARN issue listing the count and example mismatched
        files/chromosomes, or an empty list if every VCF's chromosome token
        matches a genome chromosome.
    """
    genome_chrom_set = set(genome_chroms)
    mismatches = [
        (fname, chrom)
        for fname, chrom in chrom_by_file.items()
        if chrom and chrom not in genome_chrom_set
    ]
    if not mismatches:
        return []
    examples = ", ".join(f"{fname} ('{chrom}')" for fname, chrom in mismatches[:5])
    return [
        Issue(
            WARN,
            f"{len(mismatches)} of {len(chrom_by_file)} VCF file(s) have a "
            f"chromosome token with no matching genome FASTA ({examples}"
            f"{', ...' if len(mismatches) > 5 else ''}) — these will be "
            "silently treated as reference-only for that chromosome",
        )
    ]


def check_vcf_content(vcf_path: str) -> List[Issue]:
    """Checks a single VCF's compression, header, and first data record.

    Reads only the header block plus the first data record — bounded, cheap,
    safe to run on every invocation. Mirrors what `enricher.py` actually does
    when it opens a VCF:

    - `gzip.open()` is called unconditionally (`enricher.py:19`) — a non-gzip
      file crashes immediately.
    - The header is found by scanning for the literal substring `#CHROM`
      (`enricher.py:351-355`); if absent, the header variable is never set and
      a later reference crashes (indel/haplotype mode only).
    - The header must have >=1 sample genotype column, or every variant
      silently produces zero sample associations (a sites-only VCF looks like
      it works, but nothing is ever enriched per-sample).
    - The first data record after the header is used, regardless of its
      FILTER value, to locate the `AF=` field's position in the INFO column
      (`enricher.py:356-364`); missing INFO entirely, or no `AF`-prefixed
      field, either IndexErrors or (in indel/haplotype mode) UnboundLocalErrors
      downstream.

    Args:
        vcf_path: Path to a single `.vcf.gz` file.

    Returns:
        A list of issues; empty if the file passes all checks in this
        function.
    """
    fname = os.path.basename(vcf_path)
    issues: List[Issue] = []
    try:
        with open(vcf_path, "rb") as fin:
            magic = fin.read(2)
    except OSError as e:
        return [Issue(ERROR, f"{fname}: could not read file ({e})")]
    if magic != GZIP_MAGIC:
        issues.append(
            Issue(ERROR, f"{fname}: not gzip-compressed (expected .vcf.gz)")
        )
        return issues  # can't safely read further as gzip

    try:
        with gzip.open(vcf_path, "rt") as fin:
            header_fields: Optional[List[str]] = None
            for line in fin:
                if "#CHROM" in line:
                    header_fields = line.strip().split("\t")
                    break
            if header_fields is None:
                issues.append(
                    Issue(ERROR, f"{fname}: no #CHROM header line found")
                )
                return issues
            if len(header_fields) < MIN_VCF_HEADER_FIELDS:
                issues.append(
                    Issue(
                        ERROR,
                        f"{fname}: #CHROM header has {len(header_fields)} columns, "
                        f"expected >={MIN_VCF_HEADER_FIELDS} (no sample genotype "
                        "columns — this looks like a sites-only VCF; every "
                        "variant will silently get zero sample associations)",
                    )
                )
            # first data record, regardless of FILTER (matches enricher.py
            # behavior — pos_AF is computed before the FILTER check runs)
            first_record: Optional[List[str]] = None
            for line in fin:
                if line.strip():
                    first_record = line.rstrip("\n").split("\t")
                    break
            if first_record is None:
                return issues  # empty VCF (header only); nothing more to check
            if len(first_record) < MIN_VCF_DATA_FIELDS:
                issues.append(
                    Issue(
                        ERROR,
                        f"{fname}: first data record has {len(first_record)} "
                        f"columns, expected >={MIN_VCF_DATA_FIELDS} (malformed row)",
                    )
                )
                return issues
            info_field = first_record[7]
            af_pos, af_key = _af_field(info_field)
            if af_pos is None:
                issues.append(
                    Issue(
                        ERROR,
                        f"{fname}: first data record's INFO field has no "
                        "AF-prefixed entry (e.g. 'AF=0.1') — enricher.py "
                        "locates the AF field's position exactly once, from "
                        "this same first record; if nothing matches, that "
                        "position is never assigned, and the SNP/indel "
                        "dictionary-creation step crashes immediately with a "
                        "NameError on the first PASS record, before any "
                        "genome enrichment happens",
                    )
                )
            elif af_key != "AF":
                issues.append(
                    Issue(
                        ERROR,
                        f"{fname}: first data record's INFO field matches on "
                        f"'{af_key}=', not an exact 'AF=' — crispritz locates "
                        "the AF field by prefix, not exact key, and extracts "
                        "its value with a fixed 3-character slice that only "
                        "works for a true 'AF=' key. For any longer key like "
                        f"'{af_key}', this produces a garbled, non-numeric "
                        "value for every record in this file, which crashes "
                        "CRISPRme_plots.py's pd.to_numeric() at the very end "
                        "of the run. Common with raw gnomAD-style VCFs where "
                        "a population-specific field (e.g. AF_afr) appears "
                        "before the true AF field",
                    )
                )
            elif not _af_value_parses_as_numeric(info_field, af_pos):
                issues.append(
                    Issue(
                        ERROR,
                        f"{fname}: first data record's AF field has a "
                        "non-numeric or missing value (e.g. '.', empty) — "
                        "this value is carried through with no validation "
                        "anywhere in the pipeline and crashes "
                        "CRISPRme_plots.py's pd.to_numeric() at the very end "
                        "of the run",
                    )
                )
    except OSError as e:
        issues.append(Issue(ERROR, f"{fname}: error reading VCF content ({e})"))
    return issues


def check_tbi_files(vcf_dir: str) -> List[Issue]:
    """Warns if `.tbi` index files are present in a VCF dataset directory.

    Crispritz's `.tbi` removal loop mutates the file list while iterating it
    (`crispritz.py` `genomeEnrichment()`, lines ~600-604), so some `.tbi`
    files slip through the removal and `enricher.py`'s unconditional
    `gzip.open()` then crashes or produces garbage on them.

    Args:
        vcf_dir: Path to one VCF dataset directory.

    Returns:
        A single warning `Issue` listing the stray files, or an empty list.
    """
    try:
        tbi_files = sorted(f for f in os.listdir(vcf_dir) if f.endswith(".tbi"))
    except OSError:
        return []
    if not tbi_files:
        return []
    return [
        Issue(
            WARN,
            f"{vcf_dir}: {len(tbi_files)} .tbi file(s) found "
            f"({', '.join(tbi_files[:5])}{'...' if len(tbi_files) > 5 else ''}) "
            "— remove before running; crispritz's .tbi removal has a known "
            "bug that can let some slip through",
        )
    ]


def check_pam_file(pamfile: str) -> List[Issue]:
    """Checks that the PAM file has a valid sequence + length-offset format.

    `crispritz.py` `indexGenome()` does `PAM.split()[1]` to get the signed
    length offset (line ~81) — an IndexError if there's no second token, a
    ValueError if it isn't numeric. No character validation is done at all,
    so an invalid IUPAC character is only caught much later inside the C
    search binary with no clear message.

    Args:
        pamfile: Path to the PAM file.

    Returns:
        A list of issues; empty if the file passes.
    """
    fname = os.path.basename(pamfile)
    try:
        with open(pamfile, "r") as fin:
            content = fin.read()
    except OSError as e:
        return [Issue(ERROR, f"{fname}: could not read file ({e})")]
    tokens = content.split()
    if len(tokens) < 2:
        return [
            Issue(
                ERROR,
                f"{fname}: expected '<PAM_SEQUENCE> <offset>', found "
                f"{len(tokens)} whitespace-separated token(s)",
            )
        ]
    sequence, offset = tokens[0], tokens[1]
    issues: List[Issue] = []
    try:
        int(offset)
    except ValueError:
        issues.append(
            Issue(ERROR, f"{fname}: offset '{offset}' is not a valid integer")
        )
    invalid_chars = sorted({c for c in sequence.upper() if c not in IUPAC_BASES})
    if invalid_chars:
        issues.append(
            Issue(
                ERROR,
                f"{fname}: sequence '{sequence}' contains invalid IUPAC "
                f"character(s): {', '.join(invalid_chars)}",
            )
        )
    if not sequence:
        issues.append(Issue(ERROR, f"{fname}: PAM sequence is empty"))
    return issues


def check_guide_file(guidefile: str) -> List[Issue]:
    """Checks that all guides in the guide file have the same length.

    No length-consistency check exists anywhere in the Python layer; guides
    of mixed length are fed straight to the C search binary with undefined
    behavior.

    Args:
        guidefile: Path to the guide RNA file.

    Returns:
        A list of issues; empty if the file passes.
    """
    fname = os.path.basename(guidefile)
    try:
        with open(guidefile, "r") as fin:
            guides = [line.strip() for line in fin if line.strip()]
    except OSError as e:
        return [Issue(ERROR, f"{fname}: could not read file ({e})")]
    if not guides:
        return [Issue(ERROR, f"{fname}: guide file is empty")]
    lengths = {len(g) for g in guides}
    if len(lengths) > 1:
        return [
            Issue(
                ERROR,
                f"{fname}: guides have inconsistent lengths {sorted(lengths)} "
                "— all guides must be the same length",
            )
        ]
    return []


def check_gzip_compressed(fname_path: str, label: str) -> List[Issue]:
    """Checks that a file is gzip-compressed.

    Only applies to `--gene_annotation`: `post_process.sh` unconditionally
    runs `gunzip -k -c` on it (no fallback for an already-plain-text file,
    unlike `--annotation`'s `_sort_annotation`, which now handles either).
    `gunzip` decompresses true BGZF and plain gzip identically, so checking
    for gzip specifically -- not full BGZF block structure -- is the
    complete, correct requirement here, not an approximation of one.

    Args:
        fname_path: Path to the file to check.
        label: Human-readable label for error messages (e.g. "gene
            annotation file").

    Returns:
        A list of issues; empty if the file passes.
    """
    fname = os.path.basename(fname_path)
    try:
        with open(fname_path, "rb") as fin:
            magic = fin.read(2)
    except OSError as e:
        return [Issue(ERROR, f"{label} {fname}: could not read file ({e})")]
    if magic != GZIP_MAGIC:
        return [
            Issue(
                ERROR,
                f"{label} {fname}: not gzip-compressed — compress with "
                f"'bgzip {fname}' before running",
            )
        ]
    return []


def _read_samples_file(samplesfile: str) -> List[str]:
    """Reads sample IDs (column 0) from a `--samplesID` file.

    Mirrors `associateSample.py:47-58` (installed crispritz): an optional `#`
    header line is skipped, remaining lines are tab-separated with the sample
    ID in column 0.
    """
    sample_ids: List[str] = []
    with open(samplesfile, "r") as fin:
        for line in fin:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            sample_ids.append(line.strip().split("\t")[0])
    return sample_ids


def check_samples_in_vcf_header(vcf_path: str, sample_ids: List[str]) -> List[Issue]:
    """Checks that every VCF sample column exists in the `--samplesID` file.

    `annotator.py:241` does `dict_sample_to_pop[sample]` with no fallback,
    where `dict_sample_to_pop` is built from `--samplesID`
    (`associateSample.py`). Any VCF sample name absent from that file raises
    an uncaught `KeyError` — but only inside the annotation step, i.e. after
    the full multi-hour search has already completed. Cheap to catch here:
    the VCF's `#CHROM` header and the samples file are both already read in
    full by other lightweight checks.

    Args:
        vcf_path: Path to a single `.vcf.gz` file.
        sample_ids: Sample IDs from `--samplesID` (from `_read_samples_file`).

    Returns:
        A single error `Issue` listing missing sample names, or an empty list.
    """
    fname = os.path.basename(vcf_path)
    try:
        with gzip.open(vcf_path, "rt") as fin:
            header_fields: Optional[List[str]] = None
            for line in fin:
                if "#CHROM" in line:
                    header_fields = line.strip().split("\t")
                    break
    except OSError:
        return []  # already reported by check_vcf_content
    if header_fields is None or len(header_fields) < MIN_VCF_HEADER_FIELDS:
        return []  # already reported by check_vcf_content
    vcf_samples = header_fields[9:]
    sample_id_set = set(sample_ids)
    missing = sorted(s for s in vcf_samples if s not in sample_id_set)
    if not missing:
        return []
    return [
        Issue(
            ERROR,
            f"{fname}: {len(missing)} VCF sample(s) not found in --samplesID "
            f"({', '.join(missing[:5])}{', ...' if len(missing) > 5 else ''}) "
            "— annotation will crash with a KeyError on these samples after "
            "the search completes",
        )
    ]


def run_lightweight(
    genomedir: str,
    vcf_dataset_dirs: List[str],
    guidefile: str,
    pamfile: str,
    gene_annotation_file: Optional[str] = None,
    samplesfile: Optional[str] = None,
) -> ValidationReport:
    """Runs all lightweight (always-on) input checks and returns the report.

    Reads only file headers / first few records — bounded and cheap enough to
    run unconditionally on every `complete-search` invocation, before the
    pipeline subprocess launches.

    `--annotation` is not checked here -- see the comment above the
    `gene_annotation_file` check for why (in short: `_sort_annotation` now
    accepts either compressed or plain input, so there's nothing to validate
    upfront without risking a false error).

    Args:
        genomedir: Path to the reference genome directory.
        vcf_dataset_dirs: Resolved paths to VCF dataset directories (one per
            line in the `--vcf` config file, resolved the same way
            `submit_job_automated_new_multiple_vcfs.sh` does:
            `<current_working_directory>/VCFs/<dataset_name>`). Empty if the
            search is not variant-aware.
        guidefile: Path to the guide RNA file.
        pamfile: Path to the PAM file.
        gene_annotation_file: Path to the gene annotation file, or the
            `vuoto.txt` mock path if not used.
        samplesfile: Path to the `--samplesID` file, or the `vuoto.txt` mock
            path / `None` if not used.

    Returns:
        A `ValidationReport` with every check's outcome recorded. Call
        `.has_errors()` to decide whether to abort, and `.render()` to get the
        printable report.
    """
    report = ValidationReport("CRISPRme input validator (lightweight)")

    genome_issues, genome_chroms = check_genome_fasta(genomedir)
    report.add(
        genome_issues,
        ok_message=f"Genome directory: {len(genome_chroms)} FASTA file(s), all start with '>'",
    )

    has_samplesfile = bool(samplesfile) and os.path.basename(samplesfile) != MOCK_FILENAME
    sample_ids = _read_samples_file(samplesfile) if has_samplesfile else []

    for vcf_dir in vcf_dataset_dirs:
        if not os.path.isdir(vcf_dir):
            report.add([Issue(ERROR, f"VCF dataset directory does not exist: {vcf_dir}")])
            continue
        vcf_files = [
            os.path.join(vcf_dir, f)
            for f in sorted(os.listdir(vcf_dir))
            if f.endswith(".vcf.gz")
        ]
        if not vcf_files:
            report.add(
                [Issue(ERROR, f"No .vcf.gz files found in VCF dataset directory: {vcf_dir}")]
            )
            continue

        chrom_issues, chrom_by_file = check_vcf_filename_chrom(vcf_files)
        report.add(
            chrom_issues,
            ok_message=f"{os.path.basename(vcf_dir)}: all filenames contain a chrN token",
        )
        report.add(
            check_vcf_chrom_matches_genome(chrom_by_file, genome_chroms),
            ok_message=f"{os.path.basename(vcf_dir)}: chromosome tokens match genome directory",
        )
        for vf in vcf_files:
            report.add(
                check_vcf_content(vf),
                ok_message=f"{os.path.basename(vf)}: gzip-compressed, valid #CHROM "
                "header with samples, first record has an AF field",
            )
            if has_samplesfile:
                report.add(
                    check_samples_in_vcf_header(vf, sample_ids),
                    ok_message=f"{os.path.basename(vf)}: all VCF samples found in --samplesID",
                )
        report.add(check_tbi_files(vcf_dir))

    report.add(check_pam_file(pamfile), ok_message=f"PAM file: {os.path.basename(pamfile)}")
    report.add(
        check_guide_file(guidefile), ok_message=f"Guide file: {os.path.basename(guidefile)}"
    )
    # --annotation is not checked here: _sort_annotation (crisprme.py) treats
    # a file not ending in .gz as already-plain-text, no decompression
    # attempted -- so an uncompressed --annotation file is valid input, and a
    # gzip-required check would produce a false error on it. --gene_annotation
    # has no such fallback (post_process.sh always runs gunzip on it), so it
    # still needs this check.
    if gene_annotation_file and os.path.basename(gene_annotation_file) != MOCK_FILENAME:
        report.add(
            check_gzip_compressed(gene_annotation_file, "gene annotation file"),
            ok_message=f"Gene annotation file is gzip-compressed: {os.path.basename(gene_annotation_file)}",
        )
    return report


# ---------------------------------------------------------------------------
# Full-tier checks (--full_input_validate): require reading every record in
# every VCF, so they're opt-in rather than run on every invocation.
# ---------------------------------------------------------------------------


def _normalize_chrom(chrom: str) -> str:
    """Adds a `chr` prefix if missing (e.g. '1' -> 'chr1').

    Matches the normalization already added to `pool_post_analisi_indel.py`
    for the same reason: some VCF datasets (e.g. 1000G GRCh38) store the
    `CHROM` column without a `chr` prefix even though genome indices and
    filenames use `chr`-prefixed names.
    """
    return chrom if chrom.startswith("chr") else "chr" + chrom


def _genome_chrom_length(genomedir: str, chrom: str) -> Optional[int]:
    """Reads a single-chromosome FASTA fully to compute its sequence length.

    Full-file read — only used by the opt-in full-tier scan, for the POS
    bounds check (`enricher.py:146-162` indexes `genomeStr[int(line[1])]`
    with no bounds check).

    Returns:
        Total sequence length in bases, or `None` if no matching FASTA file
        is found.
    """
    for ext in (".fa", ".fasta"):
        fpath = os.path.join(genomedir, chrom + ext)
        if os.path.isfile(fpath):
            break
    else:
        return None
    length = 0
    with open(fpath, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                continue
            length += len(line.strip())
    return length


def check_vcf_full_scan(
    vcf_path: str, expected_chrom: str, chrom_length: Optional[int]
) -> List[Issue]:
    """Reads an entire VCF once and runs every whole-file check together.

    Bundled into a single pass (rather than one read per check) since this is
    already the slow, opt-in tier — no reason to re-read a potentially huge
    file multiple times. Covers:

    - REF/ALT length survey: variants >50bp anywhere in the file slow the
      search unexpectedly (structural variants aren't only in the first
      record).
    - Chromosome consistency: every record's `CHROM` column should match the
      chromosome implied by the filename (`chrN` token) — a stray record for
      a different chromosome, e.g. from a bad VCF subset/split, wouldn't be
      caught by the filename-only lightweight check.
    - AF field-position consistency: `enricher.py` caches the AF field's
      position from the *first* record only and reuses it for every later
      record (`enricher.py:356-364`). If INFO field ordering isn't constant
      record-to-record, later records silently get the wrong AF value.
    - FILTER PASS ratio: warns if effectively no variants will survive
      enrichment, against the values the installed crispritz actually
      hardcodes (`ENRICHER_PASS_VALUES`), not `--vcf-filter-pass-values`
      (that flag isn't wired through to crispritz today).
    - POS bounds vs. chromosome length: a build/genome mismatch causes an
      IndexError in `enricher.py`, only reachable by scanning far enough into
      the file.
    - Multiallelic sites with >=10 ALT alleles: `enricher.py`'s allele-index
      matching is a substring test (`str(value) in gt_string`) that misfires
      once two-digit allele indices appear.
    - Breakend/SV notation in ALT (`]`/`[` without a leading `<`): only
      `<DEL>`-style symbolic ALT is filtered out; breakend text is spliced as
      literal sequence into the enriched genome.
    - Symbolic ALT notation (`<DEL>`, `<DUP>`, etc.): a different, milder
      failure mode than breakend -- `enricher.py` correctly filters these
      out of the indel/SV path, but does so silently, so these records are
      never searched at all with zero signal to the user (not corrupted,
      just invisibly omitted).
    - Duplicate CHR+POS records: the per-variant dict is keyed by
      `(chrom, pos)` only, so a duplicate record silently overwrites the
      first with no warning.
    - Mixed GT phasing separators (`|` vs `/`): `new_simple_analysis.py`
      decides phasing for the *entire chromosome* from whichever dict entry
      its iteration happens to hit first — a mixed phased/unphased VCF
      (common with merged cohorts) causes silently wrong haplotype/sample
      attribution. The validator can only warn here; the real fix belongs in
      `new_simple_analysis.py`.
    - Indels within 26bp of a chromosome start: `enricher.py:302`'s
      `pos-26` context window goes negative, and Python silently wraps a
      negative slice index to the end of the string instead of erroring.
    - AF/ALT allele-count consistency for multiallelic SNP sites:
      `enricher.py`'s `add_to_dict_snps` indexes the comma-separated AF list
      by each single-character ALT allele's 1-indexed position in the full
      ALT list, with no bounds check -- fewer AF values than that highest
      position raises an uncaught `IndexError` that crashes enrichment
      outright for the whole run.

    Args:
        vcf_path: Path to a single `.vcf.gz` file.
        expected_chrom: Chromosome token derived from the filename (from
            `check_vcf_filename_chrom`).
        chrom_length: Total length of the matching genome FASTA, or `None` if
            it couldn't be determined (skips the POS-bounds check).

    Returns:
        A list of issues; empty if the file passes every check.
    """
    fname = os.path.basename(vcf_path)
    issues: List[Issue] = []

    long_variant_count = 0
    rogue_chroms: Dict[str, int] = {}
    af_positions_seen: set = set()
    pass_count = 0
    total_count = 0
    multiallelic_count = 0
    breakend_count = 0
    symbolic_alt_count = 0
    seen_positions: Dict[Tuple[str, str], int] = {}
    duplicate_examples: List[str] = []
    phased_seen = False
    unphased_seen = False
    indel_near_start_count = 0
    out_of_bounds_count = 0
    af_index_overflow_count = 0
    expected_chrom_norm = _normalize_chrom(expected_chrom) if expected_chrom else None

    try:
        with gzip.open(vcf_path, "rt") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < MIN_VCF_DATA_FIELDS:
                    continue
                chrom, pos_str, _id, ref, alt, _qual, filt, info = fields[:8]
                total_count += 1

                # chromosome consistency
                if expected_chrom_norm and _normalize_chrom(chrom) != expected_chrom_norm:
                    rogue_chroms[chrom] = rogue_chroms.get(chrom, 0) + 1

                # duplicate CHR+POS
                key = (chrom, pos_str)
                if key in seen_positions:
                    seen_positions[key] += 1
                    if len(duplicate_examples) < 5:
                        duplicate_examples.append(f"{chrom}:{pos_str}")
                else:
                    seen_positions[key] = 1

                # REF/ALT length survey + breakend + multiallelic + indel proximity
                alt_alleles = alt.split(",")
                if len(alt_alleles) >= MULTIALLELIC_WARN_THRESHOLD:
                    multiallelic_count += 1
                is_breakend = any(
                    ("[" in a or "]" in a) and not a.startswith("<") for a in alt_alleles
                )
                if is_breakend:
                    breakend_count += 1
                if any(a.startswith("<") for a in alt_alleles):
                    symbolic_alt_count += 1
                if len(ref) > STRUCTURAL_VARIANT_LEN or any(
                    len(a) > STRUCTURAL_VARIANT_LEN for a in alt_alleles
                ):
                    long_variant_count += 1
                is_indel = any(len(a) != len(ref) for a in alt_alleles)
                if is_indel:
                    try:
                        if int(pos_str) < INDEL_START_PROXIMITY:
                            indel_near_start_count += 1
                    except ValueError:
                        pass

                # POS bounds
                if chrom_length is not None:
                    try:
                        if int(pos_str) < 1 or int(pos_str) > chrom_length:
                            out_of_bounds_count += 1
                    except ValueError:
                        pass

                # AF field position + FILTER pass ratio
                af_pos, _af_key = _af_field(info)
                if af_pos is not None:
                    af_positions_seen.add(af_pos)
                if filt in ENRICHER_PASS_VALUES:
                    pass_count += 1

                # AF/ALT allele-count consistency (enricher.py's multiallelic
                # SNP branch, add_to_dict_snps): a comma-separated ALT with
                # >=1 single-character allele indexes the comma-separated AF
                # list by each such allele's 1-indexed position in the full
                # ALT list, with no bounds check -- fewer AF values than the
                # highest such position raises an uncaught IndexError
                if len(ref) == 1 and len(alt_alleles) > 1 and af_pos is not None:
                    snp_positions = [
                        i + 1 for i, a in enumerate(alt_alleles) if len(a) == 1
                    ]
                    if snp_positions:
                        af_values = info.split(";")[af_pos][3:].split(",")
                        if len(af_values) < max(snp_positions):
                            af_index_overflow_count += 1

                # GT phasing separator survey
                for sample in fields[9:]:
                    gt = sample.split(":", 1)[0]
                    if "|" in gt:
                        phased_seen = True
                    elif "/" in gt:
                        unphased_seen = True
    except OSError as e:
        return [Issue(ERROR, f"{fname}: error reading VCF content ({e})")]

    if total_count == 0:
        return issues  # empty VCF (header only); nothing more to check

    if long_variant_count:
        issues.append(
            Issue(
                WARN,
                f"{fname}: {long_variant_count} variant(s) have REF/ALT >"
                f"{STRUCTURAL_VARIANT_LEN}bp — will slow the search",
            )
        )
    if rogue_chroms:
        examples = ", ".join(f"{c} ({n}x)" for c, n in list(rogue_chroms.items())[:5])
        issues.append(
            Issue(
                ERROR,
                f"{fname}: {sum(rogue_chroms.values())} record(s) have a CHROM "
                f"value other than the filename-implied '{expected_chrom}' "
                f"({examples}) — likely a bad VCF subset/split",
            )
        )
    if len(af_positions_seen) > 1:
        issues.append(
            Issue(
                WARN,
                f"{fname}: AF field position varies across records "
                f"({sorted(af_positions_seen)}) — crispritz caches the "
                "position from the first record only, so later records with "
                "a different INFO field order will get the wrong AF value",
            )
        )
    if pass_count == 0:
        # unlike the <10% case below (a subset excluded), 0% means this
        # whole dataset silently contributes nothing to enrichment -- a
        # total, guaranteed void, same severity class as the sites-only-VCF
        # lightweight check, not a partial-impact WARN
        issues.append(
            Issue(
                ERROR,
                f"{fname}: 0 of {total_count} variants have FILTER == "
                f"'{ENRICHER_PASS_VALUES[0]}' — every variant will be "
                "excluded from enrichment",
            )
        )
    elif pass_count / total_count < FILTER_PASS_WARN_RATIO:
        issues.append(
            Issue(
                WARN,
                f"{fname}: only {pass_count}/{total_count} "
                f"({100 * pass_count / total_count:.1f}%) variants have "
                f"FILTER == '{ENRICHER_PASS_VALUES[0]}' — most variants "
                "will be excluded from enrichment",
            )
        )
    if out_of_bounds_count:
        issues.append(
            Issue(
                ERROR,
                f"{fname}: {out_of_bounds_count} variant(s) have a POS outside "
                f"the chromosome bounds (length {chrom_length}) — will crash "
                "enrichment with an IndexError",
            )
        )
    if multiallelic_count:
        issues.append(
            Issue(
                WARN,
                f"{fname}: {multiallelic_count} site(s) have >="
                f"{MULTIALLELIC_WARN_THRESHOLD} ALT alleles — crispritz matches "
                "allele indices by substring, which misattributes samples once "
                "two-digit allele indices appear",
            )
        )
    if af_index_overflow_count:
        issues.append(
            Issue(
                ERROR,
                f"{fname}: {af_index_overflow_count} multiallelic SNP site(s) "
                "have fewer comma-separated AF values than crispritz's "
                "allele-index lookup needs — enricher.py indexes the AF list "
                "by each SNP ALT allele's position with no bounds check, "
                "raising an uncaught IndexError that crashes the entire "
                "enrichment step outright",
            )
        )
    if breakend_count:
        issues.append(
            Issue(
                ERROR,
                f"{fname}: {breakend_count} record(s) use breakend/SV notation "
                "in ALT (e.g. 'G]17:1234]') — only <DEL>-style symbolic ALT is "
                "filtered out, so this will be spliced as literal text into "
                "the enriched genome, corrupting it",
            )
        )
    if symbolic_alt_count:
        issues.append(
            Issue(
                WARN,
                f"{fname}: {symbolic_alt_count} record(s) use symbolic ALT "
                "notation (e.g. '<DEL>', '<DUP>') — enricher.py silently "
                "excludes these from the search entirely (not corrupted, "
                "just never searched), with no other signal that they were "
                "dropped",
            )
        )
    if duplicate_examples:
        dup_count = sum(1 for n in seen_positions.values() if n > 1)
        issues.append(
            Issue(
                WARN,
                f"{fname}: {dup_count} duplicate CHR+POS record(s) found "
                f"(e.g. {', '.join(duplicate_examples)}) — the per-variant "
                "dict is keyed by CHR+POS only, so all but the last duplicate "
                "are silently dropped",
            )
        )
    if phased_seen and unphased_seen:
        issues.append(
            Issue(
                WARN,
                f"{fname}: both phased ('|') and unphased ('/') genotypes "
                "found — new_simple_analysis.py decides haplotype-mode "
                "processing for the whole chromosome from whichever record "
                "it scans first. If that locks in phased mode, unphased "
                "samples are then incorrectly counted as variant carriers "
                "on both haplotypes regardless of their actual genotype "
                "(including homozygous-reference samples)",
            )
        )
    if indel_near_start_count:
        issues.append(
            Issue(
                WARN,
                f"{fname}: {indel_near_start_count} indel(s) within "
                f"{INDEL_START_PROXIMITY}bp of the chromosome start — "
                "enricher.py's context window goes negative here, and Python "
                "silently wraps the slice to the end of the sequence instead "
                "of erroring, corrupting that variant's FASTA record",
            )
        )
    return issues


def run_full(
    genomedir: str,
    vcf_dataset_dirs: List[str],
) -> ValidationReport:
    """Runs all full-tier (`--full_input_validate`) checks and returns the report.

    Reads every record in every VCF — slow on large datasets, so this only
    runs when explicitly requested. Meant to be called in addition to (after)
    `run_lightweight`, not instead of it.

    Args:
        genomedir: Path to the reference genome directory.
        vcf_dataset_dirs: Resolved paths to VCF dataset directories, as
            returned by `resolve_vcf_dataset_dirs`.

    Returns:
        A `ValidationReport` with every full-tier check's outcome recorded.
    """
    report = ValidationReport("CRISPRme input validator (full scan)")
    chrom_length_cache: Dict[str, Optional[int]] = {}

    for vcf_dir in vcf_dataset_dirs:
        if not os.path.isdir(vcf_dir):
            continue  # already reported by run_lightweight
        vcf_files = [
            os.path.join(vcf_dir, f)
            for f in sorted(os.listdir(vcf_dir))
            if f.endswith(".vcf.gz")
        ]
        chrom_issues, chrom_by_file = check_vcf_filename_chrom(vcf_files)
        for vf in vcf_files:
            fname = os.path.basename(vf)
            expected_chrom = chrom_by_file.get(fname, "")
            if expected_chrom not in chrom_length_cache:
                chrom_length_cache[expected_chrom] = _genome_chrom_length(
                    genomedir, expected_chrom
                )
            report.add(
                check_vcf_full_scan(
                    vf, expected_chrom, chrom_length_cache[expected_chrom]
                ),
                ok_message=f"{fname}: full scan passed (chromosome coverage, "
                "AF consistency, FILTER ratio, POS bounds, multiallelic/"
                "breakend/duplicate/phasing survey)",
            )
    return report


def _read_vcf_dataset_names(vcf_config_path: str) -> List[str]:
    """Reads dataset names from a `--vcf` config file.

    Mirrors the parsing already done inline in `complete_search()`
    (`crisprme.py:1191-1208`) so the validator resolves the same dataset
    directories the actual pipeline will use.

    Args:
        vcf_config_path: Path to the `--vcf` config file (one dataset name
            per line).

    Returns:
        A list of dataset names, one per non-empty line.
    """
    names = []
    with open(vcf_config_path, "r") as fin:
        for line in fin:
            if line.strip():
                if line[-2] == "/":
                    line = line[:-2]
                names.append(os.path.basename(line.strip()))
    return names


def resolve_vcf_dataset_dirs(
    vcf_config_path: str, current_working_directory: str
) -> List[str]:
    """Resolves a `--vcf` config file into actual VCF dataset directory paths.

    Args:
        vcf_config_path: Path to the `--vcf` config file.
        current_working_directory: CRISPRme's current working directory (the
            directory containing the `VCFs/` folder).

    Returns:
        Absolute paths to each dataset's directory under `VCFs/`.
    """
    return [
        os.path.join(current_working_directory, "VCFs", name)
        for name in _read_vcf_dataset_names(vcf_config_path)
    ]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Standalone runner for CRISPRme's input validator "
        "(for manual testing outside of complete-search)."
    )
    parser.add_argument("--genome", required=True)
    parser.add_argument("--vcf", help="--vcf config file (dataset names, one per line)")
    parser.add_argument("--cwd", default=os.getcwd(), help="current working directory (for resolving VCFs/)")
    parser.add_argument("--guide", required=True)
    parser.add_argument("--pam", required=True)
    parser.add_argument("--gene_annotation")
    parser.add_argument("--samplesID")
    parser.add_argument(
        "--full_input_validate",
        action="store_true",
        help="also run the slower, full-file-scan checks",
    )
    ns = parser.parse_args()

    vcf_dirs = resolve_vcf_dataset_dirs(ns.vcf, ns.cwd) if ns.vcf else []
    result = run_lightweight(
        ns.genome,
        vcf_dirs,
        ns.guide,
        ns.pam,
        ns.gene_annotation,
        ns.samplesID,
    )
    result.write()
    has_errors = result.has_errors()
    if ns.full_input_validate:
        full_result = run_full(ns.genome, vcf_dirs)
        full_result.write()
        has_errors = has_errors or full_result.has_errors()
    raise SystemExit(1 if has_errors else 0)
