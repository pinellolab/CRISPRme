#!/usr/bin/env python
"""Estimate the CRISPRme post-analysis resource budget for a search.

CRISPRme's post-analysis cost is dominated by the number of *in-budget*
alignments retained per site. That count grows combinatorially with the
allowed mismatches and DNA/RNA bulges, and again with local variant density
(each IUPAC/variant position can fan a single site out into a powerset of
alternative targets, see ``iupac_decomposition`` in
``new_simple_analysis.py``). Past a threshold the intermediate ``best*``
tables and the in-memory clusters explode: a single chromosome at
mm6 + 2 DNA + 2 RNA bulges produced ~207 GB of ``best*`` intermediates and
OOMed the analysis at ~15-58 GB RAM (issues #106/#107/#108).

This module is a dependency-light, standalone heuristic. It does NOT parse a
real run; it turns the search parameters into a coarse, order-of-magnitude
budget and a risk tier so a user (or a wrapper script) can be warned *before*
launching a run that will blow up.

IMPORTANT: every number produced here is a HEURISTIC estimate, not a
measurement or a guarantee. Real usage depends heavily on the genome, the
VCF/variant density, the PAM, and the number of guides. Treat the output as a
"go / be careful / probably don't" signal, not a capacity plan.

CLI usage:
    python PostProcess/search_budget.py --guide-length 23 --mm 6 --bdna 2 --brna 2
    python PostProcess/search_budget.py --mm 4 --bdna 1 --brna 1
    python PostProcess/search_budget.py --mm 6 --bdna 2 --brna 2 --fail-over EXTREME

The CLI always exits 0 (it is advisory) unless ``--fail-over TIER`` is given,
in which case it exits non-zero when the estimated tier is at or above TIER.
That opt-in hard stop is the #107 guard.
"""

from dataclasses import dataclass, asdict
from math import comb
import argparse
import sys


# --------------------------------------------------------------------------- #
# Heuristic constants. These are deliberately coarse and are documented as
# such. They were chosen so the tiers line up with the observed field data:
#   - mm<=4 & bulges<=1                 -> OK       (routine runs)
#   - moderate mm/bulges                -> HEAVY    (watch disk/RAM)
#   - mm6 + 2 DNA + 2 RNA bulges        -> EXTREME  (~200 GB disk / OOM RAM)
# --------------------------------------------------------------------------- #

# Risk tiers, ordered from cheapest to most expensive.
TIERS = ("OK", "HEAVY", "EXTREME")
_TIER_RANK = {tier: i for i, tier in enumerate(TIERS)}

# Multiplicity-score thresholds separating the tiers. The score is a
# dimensionless proxy for "how many in-budget alignments a typical site can
# fan out into" (see ``_multiplicity_score``). Tuned against the mm6+2+2
# blow-up (score ~6e10, EXTREME) and routine mm<=4 & bulges<=1+1 runs
# (score <=~3e7, OK). Escalation across the grid:
#   mm4+1+1 ~2.7e7 -> OK ; mm5+1+1 ~1.8e8 -> HEAVY ; mm6+2+1 ~8.8e9 -> HEAVY ;
#   mm6+2+2 ~6.3e10 -> EXTREME.
_HEAVY_THRESHOLD = 1.0e8
_EXTREME_THRESHOLD = 1.0e10

# Order-of-magnitude anchors for the disk/RAM heuristic, pinned to the
# reported field numbers at mm6 + 2 DNA + 2 RNA on a single large chromosome:
# ~207 GB of best* intermediates and a peak RAM footprint that OOMed in the
# ~15-58 GB range. We scale disk/RAM off the multiplicity score relative to
# that anchor point.
_ANCHOR_SCORE = None  # filled in lazily below via _multiplicity_score
_ANCHOR_DISK_GB = 207.0
_ANCHOR_RAM_GB = 50.0

# Floor estimates for a trivial run (a few hundred MB of intermediates, a
# couple GB of working memory), so tiny searches don't report ~0.
_FLOOR_DISK_GB = 0.3
_FLOOR_RAM_GB = 2.0


@dataclass
class BudgetEstimate:
    """Coarse, heuristic post-analysis budget for a CRISPRme search.

    All fields except ``risk`` are order-of-magnitude ESTIMATES.
    """

    guide_length: int
    mismatches: int
    bulges_dna: int
    bulges_rna: int
    variants_per_kb: float
    # Dimensionless proxy for in-budget alignment fan-out per site.
    multiplicity_score: float
    # Order-of-magnitude peak intermediate disk footprint (best* tables etc.).
    peak_disk_gb: float
    # Order-of-magnitude peak working-set RAM.
    peak_ram_gb: float
    # One of TIERS.
    risk: str

    def to_dict(self):
        return asdict(self)


def _edit_configurations(guide_length, mismatches, bulges_dna, bulges_rna):
    """Count, roughly, the distinct edited-alignment configurations.

    This is the combinatorial core of the blow-up. For a guide of length L we
    approximate the number of ways to place up to ``mismatches`` substitutions
    and up to ``bulges_dna``/``bulges_rna`` bulges along the alignment.

    - Mismatches: choose positions C(L, k) and (implicitly) a substituted base.
      We fold the per-position base choices into a small multiplier rather than
      the full 3**k so the score stays in a tractable range; the point is the
      relative growth, not an exact count.
    - Bulges: each additional DNA/RNA bulge adds roughly L placement choices
      and lengthens the effective search window, so we treat total bulges as a
      further combinatorial factor over an extended length.
    """
    L = max(1, int(guide_length))
    mm = max(0, int(mismatches))
    bdna = max(0, int(bulges_dna))
    brna = max(0, int(bulges_rna))

    # Mismatch placements (cumulative up to mm), with a modest per-mismatch
    # base-choice multiplier to reflect A/C/G/T substitutions without the full
    # 3**k explosion.
    mm_conf = 0.0
    for k in range(0, mm + 1):
        mm_conf += comb(L, min(k, L)) * (1.7 ** k)

    # Bulge placements: an extended alignment window of length L+bulges, with
    # cumulative placement of up to (bdna+brna) bulges over it. DNA and RNA
    # bulges compound because they can co-occur.
    total_bulges = bdna + brna
    window = L + total_bulges
    bulge_conf = 0.0
    for b in range(0, total_bulges + 1):
        bulge_conf += comb(window, min(b, window))

    return mm_conf * bulge_conf


def _multiplicity_score(guide_length, mismatches, bulges_dna, bulges_rna,
                        variants_per_kb):
    """Dimensionless proxy for in-budget alignments retained per site.

    Combines the edited-alignment configuration count with a variant-density
    fan-out factor. Variant density matters because each variant/IUPAC
    position inside a target can split it into a powerset of alternative
    targets during ``iupac_decomposition`` -- so density enters exponentially
    over the (small) expected number of variants overlapping a target window.
    """
    conf = _edit_configurations(guide_length, mismatches, bulges_dna,
                                bulges_rna)

    # Expected variants overlapping a target window of ~guide_length bases.
    L = max(1, int(guide_length))
    exp_vars = (float(variants_per_kb) / 1000.0) * L
    # Powerset-style fan-out, capped so a dense region doesn't produce an
    # astronomically large (and meaningless) score.
    variant_factor = min(2.0 ** exp_vars, 64.0)

    return conf * variant_factor


# Compute the anchor score once, for the mm6 + 2 DNA + 2 RNA / L=23 case with
# no variant fan-out, so the disk/RAM scaling is pinned to field data.
_ANCHOR_SCORE = _multiplicity_score(23, 6, 2, 2, 0.0)


def _risk_tier(score):
    if score >= _EXTREME_THRESHOLD:
        return "EXTREME"
    if score >= _HEAVY_THRESHOLD:
        return "HEAVY"
    return "OK"


def estimate_budget(guide_length, mismatches, bulges_dna, bulges_rna,
                    variants_per_kb=None):
    """Estimate the post-analysis resource budget for a CRISPRme search.

    Parameters
    ----------
    guide_length : int
        Length of the guide/spacer (e.g. 20 or 23).
    mismatches : int
        Maximum mismatches allowed in the search (mm).
    bulges_dna : int
        Maximum DNA bulges allowed (bDNA).
    bulges_rna : int
        Maximum RNA bulges allowed (bRNA).
    variants_per_kb : float, optional
        Local variant density (variants per kilobase) in the searched region.
        Drives the IUPAC/variant fan-out. Defaults to 0 (reference-only search)
        when None.

    Returns
    -------
    BudgetEstimate
        Dataclass with the multiplicity score, order-of-magnitude peak disk and
        peak RAM estimates (both clearly heuristic), and a risk tier in
        {"OK", "HEAVY", "EXTREME"}.

    Notes
    -----
    Thresholding intent (documented, not a guarantee):
      - mm <= 4 and (bDNA + bRNA) <= 1  -> OK
      - escalating with mm/bulges       -> HEAVY
      - mm 6 + 2 DNA + 2 RNA            -> EXTREME (~200 GB disk / OOM RAM)
    """
    if variants_per_kb is None:
        variants_per_kb = 0.0

    score = _multiplicity_score(guide_length, mismatches, bulges_dna,
                                bulges_rna, variants_per_kb)
    tier = _risk_tier(score)

    # Order-of-magnitude disk/RAM: scale linearly off the field-anchored point.
    # Sub-linear anchoring would understate the tail, so we keep it linear in
    # the score ratio, then apply floors for trivial runs.
    ratio = score / _ANCHOR_SCORE if _ANCHOR_SCORE else 0.0
    peak_disk_gb = max(_FLOOR_DISK_GB, ratio * _ANCHOR_DISK_GB)
    peak_ram_gb = max(_FLOOR_RAM_GB, ratio * _ANCHOR_RAM_GB)

    return BudgetEstimate(
        guide_length=int(guide_length),
        mismatches=int(mismatches),
        bulges_dna=int(bulges_dna),
        bulges_rna=int(bulges_rna),
        variants_per_kb=float(variants_per_kb),
        multiplicity_score=score,
        peak_disk_gb=peak_disk_gb,
        peak_ram_gb=peak_ram_gb,
        risk=tier,
    )


def _format_gb(value):
    """Human-friendly, order-of-magnitude GB formatting."""
    if value >= 1000.0:
        return f"~{value / 1000.0:.1f} TB"
    if value >= 10.0:
        return f"~{value:.0f} GB"
    return f"~{value:.1f} GB"


def format_report(est):
    """Render a human-readable, multi-line report for a BudgetEstimate."""
    lines = []
    lines.append("CRISPRme post-analysis budget estimate (HEURISTIC)")
    lines.append("=" * 52)
    lines.append(f"  guide length     : {est.guide_length}")
    lines.append(f"  mismatches (mm)  : {est.mismatches}")
    lines.append(f"  DNA bulges       : {est.bulges_dna}")
    lines.append(f"  RNA bulges       : {est.bulges_rna}")
    lines.append(f"  variants / kb    : {est.variants_per_kb:g}")
    lines.append("-" * 52)
    lines.append(f"  alignment-multiplicity score : {est.multiplicity_score:.3g}")
    lines.append(f"  est. peak disk (best* etc.)  : {_format_gb(est.peak_disk_gb)}  [heuristic]")
    lines.append(f"  est. peak RAM                : {_format_gb(est.peak_ram_gb)}  [heuristic]")
    lines.append(f"  risk tier                    : {est.risk}")
    lines.append("-" * 52)

    if est.risk == "OK":
        lines.append("  OK: routine search space. No special action expected.")
    else:
        lines.append(f"  WARNING: {est.risk} search space.")
        lines.append("  This run may consume a very large amount of disk and RAM")
        lines.append("  in the post-analysis stage and could OOM. Guidance:")
        lines.append("    - Reduce mismatches and/or DNA/RNA bulges if possible")
        lines.append("      (cost grows combinatorially with each one).")
        lines.append(f"    - Expect on the order of {_format_gb(est.peak_disk_gb)} of")
        lines.append(f"      best* intermediates and up to {_format_gb(est.peak_ram_gb)} RAM.")
        lines.append("    - Consider running per-chromosome to cap peak footprint.")
        lines.append("    - Ensure ample scratch disk and swap before launching.")
        if est.risk == "EXTREME":
            lines.append("    - EXTREME: at this setting a single large chromosome has")
            lines.append("      been observed to produce ~207 GB of intermediates and")
            lines.append("      OOM. Strongly reconsider the parameters.")
    lines.append("")
    lines.append("  NOTE: all disk/RAM figures are order-of-magnitude estimates,")
    lines.append("  not guarantees. Actual usage depends on genome, VCF density,")
    lines.append("  PAM, and guide count.")
    return "\n".join(lines)


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="search_budget.py",
        description=(
            "Estimate the CRISPRme post-analysis resource budget for a search "
            "and warn before a run that will blow up (issues #106/#107/#108). "
            "All disk/RAM figures are heuristic order-of-magnitude estimates."
        ),
    )
    parser.add_argument("--guide-length", type=int, default=23,
                        help="Guide/spacer length (default: 23).")
    parser.add_argument("--mm", type=int, required=True,
                        help="Maximum mismatches allowed.")
    parser.add_argument("--bdna", type=int, default=0,
                        help="Maximum DNA bulges allowed (default: 0).")
    parser.add_argument("--brna", type=int, default=0,
                        help="Maximum RNA bulges allowed (default: 0).")
    parser.add_argument("--variants-per-kb", type=float, default=None,
                        help="Local variant density (variants/kb) driving "
                             "IUPAC fan-out. Default: reference-only (0).")
    parser.add_argument("--fail-over", choices=TIERS, default=None,
                        help="Opt-in hard stop (#107): exit non-zero when the "
                             "estimated risk tier is at or above this tier. "
                             "Without this flag the tool is purely advisory "
                             "and always exits 0.")
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    est = estimate_budget(
        guide_length=args.guide_length,
        mismatches=args.mm,
        bulges_dna=args.bdna,
        bulges_rna=args.brna,
        variants_per_kb=args.variants_per_kb,
    )

    print(format_report(est))

    if args.fail_over is not None:
        if _TIER_RANK[est.risk] >= _TIER_RANK[args.fail_over]:
            print(
                f"\nERROR: risk tier {est.risk} >= --fail-over {args.fail_over}; "
                "exiting non-zero.",
                file=sys.stderr,
            )
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
