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
import subprocess

import pandas as pd

PRED_COLS = [
    "Spacer+PAM", "Chromosome", "Start_coordinate_(fewest_mm+b)",
    "Strand_(fewest_mm+b)", "Aligned_spacer+PAM_(fewest_mm+b)",
    "Aligned_protospacer+PAM_REF_(fewest_mm+b)", "Aligned_protospacer+PAM_ALT_(fewest_mm+b)",
    "Mismatches_(fewest_mm+b)", "Bulges_(fewest_mm+b)", "CFD_score_(fewest_mm+b)",
]
CHUNKSIZE = 10**6


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


def build_offtarget_bed(preds: pd.DataFrame, ucsc_to_genbank: Dict[str, str], bed_path: str) -> str:
    """Writes a BED file of one haplotype's predicted off-target coordinates,
    keyed by each row's `off_target_id`, in the chain file's chromosome
    naming (GenBank accessions, per the HPRC chromAlias convention).
    """
    bed = preds[["Chromosome", "Start_coordinate_(fewest_mm+b)", "off_target_id"]].copy()
    bed = bed.rename(columns={"Start_coordinate_(fewest_mm+b)": "chromEnd"})
    bed["chromStart"] = bed["chromEnd"] - 1
    bed["chrom"] = bed["Chromosome"].map(ucsc_to_genbank)
    bed = bed.dropna(subset=["chrom"])
    bed = bed[["chrom", "chromStart", "chromEnd", "off_target_id"]]
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    return bed_path


def run_liftover(
    bed_path: str, chain_file: str, mapped_path: str, unmapped_path: str, env: str = "liftover_env"
) -> Tuple[str, str]:
    """Runs UCSC `liftOver` (via the dedicated `liftover_env` conda env --
    not a `crisprme-*` env dependency) on a BED file.

    Raises:
        RuntimeError: If `liftOver` exits non-zero.
    """
    cmd = f"liftOver {bed_path} {chain_file} {mapped_path} {unmapped_path}"
    result = subprocess.run(f"mamba run -n {env} {cmd}", shell=True, capture_output=True, text=True)
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


def reconcile_haplotypes(
    haplotypes: Dict[str, dict], workdir: str, merge_bp: int = 3
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Runs the full reconciliation pipeline for exactly two haplotypes.

    Args:
        haplotypes: Exactly two entries, keyed by haplotype name (e.g.
            "paternal", "maternal"). Each value must have:
            `chrom_alias_file`, `chain_file`, `results_dir`. `results_prefix`
            is derived automatically via `find_results_prefix` if not given.
        workdir: Directory for intermediate BED files (created if missing).
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

    predictions: Dict[str, pd.DataFrame] = {}
    lifted: Dict[str, pd.DataFrame] = {}
    unlifted_ids: Dict[str, set] = {}

    for name, cfg in haplotypes.items():
        prefix = cfg.get("results_prefix") or find_results_prefix(cfg["results_dir"])
        _, ucsc_to_genbank = load_chrom_alias(cfg["chrom_alias_file"])

        predictions[name] = load_crisprme_predictions(cfg["results_dir"], prefix, merge_bp)

        bed_path = os.path.join(workdir, f"{name}_offtargets.bed")
        build_offtarget_bed(predictions[name], ucsc_to_genbank, bed_path)

        mapped_path = os.path.join(workdir, f"{name}_offtargets_lifted.bed")
        unmapped_path = os.path.join(workdir, f"{name}_offtargets_not_lifted.bed")
        run_liftover(bed_path, cfg["chain_file"], mapped_path, unmapped_path)

        lifted[name] = load_lifted_bed(mapped_path)
        unlifted_ids[name] = load_unlifted_ids(unmapped_path)

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
    return combined, summary
