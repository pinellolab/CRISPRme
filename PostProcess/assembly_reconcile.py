"""
This module reconciles CRISPRme off-target predictions across two personal
genome haplotype assemblies (e.g. paternal/maternal), producing one combined
report in hg38 coordinates.

Each haplotype's predicted off-targets are lifted independently to hg38 via
`liftOver` (using that haplotype's own chain file), then joined directly on
hg38 coordinates: found on both haplotypes is homozygous-equivalent, found on
only one is heterozygous-equivalent, and predictions that don't lift at all
are haplotype-non-mappable -- invisible to any reference-genome-based search.

This logic was developed and validated interactively against real HG01255
(HPRC) data before being ported here -- see
`assembly_search/assembly_search_generalized_071326.ipynb` in the project
root for the exploratory version and the real numbers it produced. Two real
bugs were caught during that validation and are fixed here from the start:

1. The alternative-alignments file contains many rows per genomic locus with
   the *exact* same coordinate (different mismatch/bulge interpretations of
   one site -- measured up to 233 rows at one locus).
2. A subtler case: ~5.7% of loci have a near-duplicate 1-3bp away on the same
   chrom+strand, every one involving a bulge -- the *same* physical site
   reported at a shifted anchor position because a bulge changes the
   alignment's registration.

Both are handled by `cluster_collapse()`, which reuses the exact greedy
chained-gap clustering algorithm CRISPRme's own `--merge` step uses
(`merge_contiguous_targets.py:531-596`), rather than exact-coordinate
matching. The `merge_bp` passed to every function in this module should match
the `--merge` value used for the underlying `complete-search` runs, since
this is reusing CRISPRme's own definition of "one site."
"""

from typing import Dict, List, Optional, Tuple

import glob
import os
import shutil
import subprocess

import pandas as pd

PRED_COLS = [
    "Spacer+PAM", "Chromosome", "Start_coordinate_(fewest_mm+b)",
    "Strand_(fewest_mm+b)", "Aligned_spacer+PAM_(fewest_mm+b)",
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)", "Aligned_protospacer+PAM_ALT_(fewest_mm+b)",
    "Mismatches_(fewest_mm+b)", "Bulges_(fewest_mm+b)", "CFD_score_(fewest_mm+b)",
]
CHUNKSIZE = 10**6
# if more than this fraction of a haplotype's predictions have a chromosome
# missing from the chromAlias file, that's almost certainly a genome/alias
# naming mismatch, not a handful of legitimately-unmappable decoy contigs --
# raise instead of silently folding all of them into "non-mappable"
CHROM_ALIAS_MISMATCH_ERROR_RATIO = 0.5


def cluster_collapse(
    df: pd.DataFrame, chrom_col: str, strand_col: str, pos_col: str,
    score_col: str, merge_bp: int,
) -> pd.DataFrame:
    """Collapses rows to one per proximity cluster.

    Uses the same greedy chained-gap algorithm CRISPRme's own `--merge` step
    uses: sort by position, start a new cluster whenever the gap from the
    previous point (same chrom+strand) exceeds `merge_bp` -- clusters chain,
    this is not a fixed window -- then keep the best-scoring row per cluster.

    Args:
        df: Rows to collapse.
        chrom_col: Column name holding the chromosome.
        strand_col: Column name holding the strand.
        pos_col: Column name holding the position to cluster on.
        score_col: Column name to pick the best row per cluster (higher wins;
            NaN sorts last).
        merge_bp: Maximum gap (bp) between consecutive points to stay in the
            same cluster. Should match the `--merge` value used to generate
            the underlying CRISPRme results.

    Returns:
        One row per cluster, the highest-`score_col` row in each.
    """
    df = df.sort_values([chrom_col, strand_col, pos_col]).reset_index(drop=True)
    cluster_ids: List[int] = []
    cluster_id = 0
    prev_chrom = prev_strand = prev_pos = None
    for chrom, strand, pos in zip(df[chrom_col], df[strand_col], df[pos_col]):
        if chrom != prev_chrom or strand != prev_strand or (pos - prev_pos) > merge_bp:
            cluster_id += 1
        cluster_ids.append(cluster_id)
        prev_chrom, prev_strand, prev_pos = chrom, strand, pos
    df = df.assign(_cluster_id=cluster_ids)
    df = df.sort_values(score_col, ascending=False, na_position="last")
    df = df.drop_duplicates(subset=["_cluster_id"], keep="first")
    return df.drop(columns=["_cluster_id"]).reset_index(drop=True)


def find_results_prefix(results_dir: str) -> str:
    """Derives a haplotype's CRISPRme results filename prefix by locating its
    `*_integrated_results.tsv` file, rather than reconstructing CRISPRme's
    internal naming convention (guide+PAM+genome+mm+bMax) by hand.

    Args:
        results_dir: A `complete-search` output folder.

    Returns:
        The filename prefix shared by `<prefix>_integrated_results.tsv` and
        `<prefix>_all_results_with_alternative_alignments.tsv`.

    Raises:
        FileNotFoundError: If zero or more than one `*_integrated_results.tsv`
            is found (e.g. a multi-guide run, not yet supported here).
    """
    suffix = "_integrated_results.tsv"
    matches = glob.glob(os.path.join(results_dir, f"*{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one *{suffix} file in {results_dir}, found "
            f"{len(matches)}: {matches}"
        )
    return os.path.basename(matches[0])[: -len(suffix)]


# submit_job_automated_new_multiple_vcfs.sh renames log_error.txt to this,
# unconditionally, as its very last meaningful step (right after echoing
# "JOB END") -- so its presence is the one signal that a run reached its
# true end, not just that some result files happened to get written along
# the way.
LOG_ERROR_NO_CHECK_FILENAME = "log_error_no_check.txt"


def haplotype_search_complete(results_dir: str) -> bool:
    """Checks whether a haplotype's `complete-search` output represents a
    genuinely finished run, not just one that produced some result files
    before failing partway through.

    `*_integrated_results.tsv` (what `find_results_prefix` looks for) is
    written well before the pipeline's actual end -- a crash during the
    zip/db_creation/mail-send/cleanup steps that follow it would still leave
    that file behind. `log_error.txt` only gets renamed to
    `log_error_no_check.txt` as the pipeline's unconditional last step, so
    requiring both together is the accurate "did this really finish" check,
    used to decide whether `assembly_search` can skip re-running an
    already-completed haplotype search.

    Args:
        results_dir: A `complete-search` output folder.

    Returns:
        True if both the integrated results file and the post-completion
        log rename are present.
    """
    if not os.path.isdir(results_dir):
        return False
    try:
        find_results_prefix(results_dir)
    except FileNotFoundError:
        return False
    return os.path.isfile(os.path.join(results_dir, LOG_ERROR_NO_CHECK_FILENAME))


def load_chrom_alias(chrom_alias_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Loads a chromAlias file into assembly->ucsc and ucsc->genbank mappings.

    Args:
        chrom_alias_file: Path to an HPRC-style `*.chromAlias.txt` file
            (tab-separated, columns `# assembly`, `ucsc`, `genbank`).

    Returns:
        `(assembly_to_ucsc, ucsc_to_genbank)`.
    """
    chrom_alias_df = pd.read_csv(chrom_alias_file, sep="\t")
    assembly_to_ucsc = dict(zip(chrom_alias_df["# assembly"], chrom_alias_df["ucsc"]))
    ucsc_to_genbank = dict(zip(chrom_alias_df["ucsc"], chrom_alias_df["genbank"]))
    return assembly_to_ucsc, ucsc_to_genbank


def load_crisprme_predictions(
    results_dir: str, prefix: str, merge_bp: int, cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """Loads and combines one haplotype's CRISPRme output into one
    deduplicated, one-row-per-genomic-locus table.

    Args:
        results_dir: A `complete-search` output folder.
        prefix: Filename prefix (from `find_results_prefix`).
        merge_bp: Passed to `cluster_collapse` -- should match the `--merge`
            value used for this haplotype's `complete-search` run.
        cols: Result columns to keep. Defaults to `PRED_COLS`.

    Returns:
        One row per genomic locus, with a stable `off_target_id` column.
    """
    if cols is None:
        cols = PRED_COLS

    alt_file = os.path.join(results_dir, f"{prefix}_all_results_with_alternative_alignments.tsv")
    integrated_file = os.path.join(results_dir, f"{prefix}_integrated_results.tsv")

    with open(alt_file) as f:
        num_lines = sum(1 for _ in f)
    num_chunks = (num_lines - 1) // CHUNKSIZE + 1
    chunks = [
        chunk
        for chunk in pd.read_csv(alt_file, sep="\t", usecols=cols, chunksize=CHUNKSIZE)
    ]
    alt_df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=cols)

    integrated_df = pd.read_csv(integrated_file, sep="\t", usecols=cols)
    combined = alt_df.merge(integrated_df, on=cols, how="outer").drop_duplicates(ignore_index=True)

    combined = cluster_collapse(
        combined, "Chromosome", "Strand_(fewest_mm+b)",
        "Start_coordinate_(fewest_mm+b)", "CFD_score_(fewest_mm+b)", merge_bp,
    )
    combined["off_target_id"] = combined.index.astype(str)
    return combined


def check_chrom_alias_coverage(
    preds: pd.DataFrame, ucsc_to_genbank: Dict[str, str], name: str
) -> None:
    """Sanity-checks that the chromAlias file actually names this haplotype's
    chromosomes, before treating any mismatch as ordinary haplotype-private
    biology.

    A handful of predictions on contigs the chromAlias file doesn't cover
    (e.g. small unplaced/decoy scaffolds) is expected and gets folded into
    the non-mappable count downstream. But if a *large* fraction of
    predictions have an unrecognized chromosome, that's a strong signal the
    genome FASTA's chromosome naming doesn't actually match the supplied
    chromAlias file (wrong file, wrong assembly, or a naming convention the
    file doesn't cover) -- worth failing clearly and immediately rather than
    silently reporting a mostly-meaningless "almost everything is
    haplotype-private" result.

    Args:
        preds: This haplotype's predictions (must have a `Chromosome` column).
        ucsc_to_genbank: This haplotype's chromAlias `ucsc`->`genbank` mapping.
        name: Haplotype name, for the error message.

    Raises:
        ValueError: If the unmatched-row fraction exceeds
            `CHROM_ALIAS_MISMATCH_ERROR_RATIO`.
    """
    total = len(preds)
    if total == 0:
        return
    known = set(ucsc_to_genbank)
    unmatched = preds.loc[~preds["Chromosome"].isin(known)]
    ratio = len(unmatched) / total
    if ratio > CHROM_ALIAS_MISMATCH_ERROR_RATIO:
        examples = sorted(unmatched["Chromosome"].unique())[:5]
        raise ValueError(
            f"{name}: {len(unmatched)}/{total} predictions ({ratio:.0%}) have a "
            f"chromosome not found in the {name} chromAlias file's 'ucsc' column "
            f"(e.g. {examples}) -- this usually means the genome FASTA's "
            f"chromosome naming doesn't match the --chrom-alias-{name} file "
            "provided, not that these are genuinely unmappable contigs. Both "
            "haplotype searches already completed successfully; only "
            "reconciliation is affected by this. Fix the chromAlias file or "
            "genome naming, then re-run assembly-search -- already-completed "
            "haplotype searches will be reused, not re-run."
        )


def build_offtarget_bed(
    preds: pd.DataFrame, ucsc_to_genbank: Dict[str, str], bed_path: str
) -> Tuple[str, set]:
    """Writes a BED file of one haplotype's predicted off-target coordinates,
    keyed by each row's `off_target_id`, in the chain file's chromosome
    naming (GenBank accessions, per the HPRC chromAlias convention).

    Returns:
        `(bed_path, dropped_ids)` -- `dropped_ids` are the `off_target_id`s
        excluded because their chromosome has no chromAlias mapping. These
        are just as unmappable-to-hg38 as a genuine liftOver failure (there's
        no way to look them up in the chain file at all), so callers should
        fold them into the same non-mappable accounting as
        `load_unlifted_ids`, not just discard them.
    """
    bed = preds[["Chromosome", "Start_coordinate_(fewest_mm+b)", "off_target_id"]].copy()
    bed = bed.rename(columns={"Start_coordinate_(fewest_mm+b)": "chromEnd"})
    bed["chromStart"] = bed["chromEnd"] - 1
    bed["chrom"] = bed["Chromosome"].map(ucsc_to_genbank)
    dropped_ids = set(bed.loc[bed["chrom"].isna(), "off_target_id"])
    bed = bed.dropna(subset=["chrom"])
    bed = bed[["chrom", "chromStart", "chromEnd", "off_target_id"]]
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    return bed_path, dropped_ids


def check_liftover_available() -> None:
    """Verifies the `liftOver` binary is on PATH before any expensive work
    starts.

    `liftOver` (UCSC tool, bioconda package `ucsc-liftover`) is a required
    dependency for `assembly-search` -- bundled into the Docker image
    alongside crispritz/crisprme, but not yet part of CRISPRme's own bioconda
    recipe for the plain `mamba install crisprme` path (a separate,
    not-yet-done follow-up). Checking this before launching two multi-hour
    haplotype searches means a missing dependency fails immediately and
    clearly instead of only surfacing at the very last reconciliation step.

    Raises:
        RuntimeError: If `liftOver` isn't found on PATH.
    """
    if shutil.which("liftOver") is None:
        raise RuntimeError(
            "liftOver (UCSC tool) not found on PATH -- required for "
            "assembly-search. Install it with `mamba install -c bioconda "
            "ucsc-liftover` (or `conda install`), or use a CRISPRme Docker "
            "image built after this dependency was added."
        )


def run_liftover(bed_path: str, chain_file: str, mapped_path: str, unmapped_path: str) -> Tuple[str, str]:
    """Runs UCSC `liftOver` on a BED file. Assumes `liftOver` is on PATH --
    see `check_liftover_available`.

    Raises:
        RuntimeError: If `liftOver` exits non-zero.
    """
    cmd = f"liftOver {bed_path} {chain_file} {mapped_path} {unmapped_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"liftOver failed:\n{result.stderr}")
    return mapped_path, unmapped_path


def load_lifted_bed(mapped_path: str) -> pd.DataFrame:
    return pd.read_csv(
        mapped_path, sep="\t", header=None,
        names=["hg38_chr", "hg38_start", "hg38_end", "off_target_id"],
        dtype={"off_target_id": str},
    )


def load_unlifted_ids(unmapped_path: str) -> set:
    """liftOver writes failed features back out in BED format, with
    '#'-prefixed comment lines explaining why each one failed to map."""
    ids = set()
    with open(unmapped_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.strip().split("\t")[3])
    return ids


LOG_VERBOSE_FILENAME = "log_verbose.txt"
LOG_ERROR_FILENAME = "log_error.txt"


def reconcile_haplotypes(
    haplotypes: Dict[str, dict], workdir: str, merge_bp: int = 3
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Runs the full reconciliation pipeline for exactly two haplotypes.

    Writes its own progress to `<workdir>/log_verbose.txt` and, on failure,
    the error to `<workdir>/log_error.txt` -- the same filenames every other
    CRISPRme run (`complete-search`, `complete-test`) already writes, opened
    the same way `complete_search()` opens them: truncated fresh at the
    start of every call, not accumulated (`crisprme.py`'s own
    `open(..., "w")` on both files, right before launching the pipeline).
    `complete_search()` never actually re-enters an existing output folder
    (an earlier check blocks that), but `reconcile_haplotypes` deliberately
    does support being retried into the same `workdir` (see
    `haplotype_search_complete` in this module and its use in
    `crisprme.py`'s `assembly_search`), so each retry's log reflects only
    that attempt.

    Args:
        haplotypes: Exactly two entries, keyed by haplotype name (e.g.
            "paternal", "maternal"). Each value must have:
            `chrom_alias_file`, `chain_file`, `results_dir`. `results_prefix`
            is derived automatically via `find_results_prefix` if not given.
        workdir: Directory for intermediate BED files and the log files
            described above (created if missing).
        merge_bp: Locus-clustering threshold -- should match the `--merge`
            value used for both haplotypes' `complete-search` runs.

    Returns:
        `(combined, summary)` where `combined` is one row per reconciled
        locus with an `origin` column (`"<a>_only"`, `"<b>_only"`, or
        `"both"`), and `summary` has counts per category plus each
        haplotype's non-mappable (un-liftable) count.

    Raises:
        ValueError: If `haplotypes` doesn't have exactly two entries.
    """
    if len(haplotypes) != 2:
        raise ValueError(
            f"reconcile_haplotypes requires exactly 2 haplotypes, got {len(haplotypes)} "
            "(N-way ploidy is not yet supported)"
        )
    os.makedirs(workdir, exist_ok=True)
    log_verbose_path = os.path.join(workdir, LOG_VERBOSE_FILENAME)
    log_error_path = os.path.join(workdir, LOG_ERROR_FILENAME)
    open(log_verbose_path, "w").close()
    open(log_error_path, "w").close()

    def log(msg: str) -> None:
        with open(log_verbose_path, "a") as f:
            f.write(msg + "\n")

    try:
        predictions: Dict[str, pd.DataFrame] = {}
        lifted: Dict[str, pd.DataFrame] = {}
        unlifted_ids: Dict[str, set] = {}

        for name, cfg in haplotypes.items():
            log(f"Loading {name} predictions...")
            prefix = cfg.get("results_prefix") or find_results_prefix(cfg["results_dir"])
            _, ucsc_to_genbank = load_chrom_alias(cfg["chrom_alias_file"])

            predictions[name] = load_crisprme_predictions(cfg["results_dir"], prefix, merge_bp)
            check_chrom_alias_coverage(predictions[name], ucsc_to_genbank, name)
            log(f"{name}: {len(predictions[name])} unique predicted loci")

            bed_path = os.path.join(workdir, f"{name}_offtargets.bed")
            _, no_chrom_alias_ids = build_offtarget_bed(predictions[name], ucsc_to_genbank, bed_path)

            log(f"Lifting {name} predictions to hg38...")
            mapped_path = os.path.join(workdir, f"{name}_offtargets_lifted.bed")
            unmapped_path = os.path.join(workdir, f"{name}_offtargets_not_lifted.bed")
            run_liftover(bed_path, cfg["chain_file"], mapped_path, unmapped_path)

            lifted[name] = load_lifted_bed(mapped_path)
            # a missing chromAlias entry is just as unmappable-to-hg38 as a
            # genuine liftOver failure -- there's no chain-file name to look
            # it up under at all -- so it belongs in the same non-mappable
            # count, not silently dropped (see build_offtarget_bed's docstring)
            unlifted_ids[name] = load_unlifted_ids(unmapped_path) | no_chrom_alias_ids
            log(f"{name}: {len(unlifted_ids[name])} non-mappable (no hg38 equivalent)")

        hg38_predictions: Dict[str, pd.DataFrame] = {}
        for name in haplotypes:
            preds = predictions[name].copy()
            preds["off_target_id"] = preds["off_target_id"].astype(str)
            merged = preds.merge(lifted[name], on="off_target_id", how="inner")
            merged = cluster_collapse(
                merged, "hg38_chr", "Strand_(fewest_mm+b)", "hg38_start",
                "CFD_score_(fewest_mm+b)", merge_bp,
            )
            hg38_predictions[name] = merged

        names = list(haplotypes.keys())
        a, b = names[0], names[1]

        log("Combining lifted predictions across haplotypes...")
        combined = hg38_predictions[a].merge(
            hg38_predictions[b],
            on=["hg38_chr", "hg38_start", "hg38_end"],
            how="outer",
            suffixes=(f"_{a}", f"_{b}"),
            indicator=True,
        )
        origin_map = {"left_only": f"{a}_only", "right_only": f"{b}_only", "both": "both"}
        combined["origin"] = combined["_merge"].map(origin_map)
        combined = combined.drop(columns=["_merge"])

        origin_counts = combined["origin"].value_counts()
        summary = {
            "both": int(origin_counts.get("both", 0)),
            f"{a}_only": int(origin_counts.get(f"{a}_only", 0)),
            f"{b}_only": int(origin_counts.get(f"{b}_only", 0)),
            f"{a}_non_mappable": len(unlifted_ids[a]),
            f"{b}_non_mappable": len(unlifted_ids[b]),
        }
        log("Reconciliation complete:")
        for category, count in summary.items():
            log(f"  {category}: {count}")
    except Exception as e:
        with open(log_error_path, "a") as f:
            f.write(f"{e}\n")
        raise
    return combined, summary
