# CRISPRme+ — changes emerging from the 2.1.x lineage

This document tracks everything that diverges from the last stable **2.1.x** release
as the **2.2.0 / CRISPRme+** lineage takes shape. It is the single place to see
"what is new," "what is faster," and — critically — **what breaks backward
compatibility**. Version stays `2.2.0` and the package/CLI stays `crisprme` (so
bioconda/Docker users are not disrupted); the "+" denotes the diverging,
optimization-focused lineage that **may break retro-compatibility** where noted.

Status legend: ✅ done & validated · 🔬 validated in isolation, integration pending ·
🚧 in progress · ⏸ pinned (waiting on an external dependency).

---

## 1. Correctness — variant / PAM handling

### 1.1 IUPAC-aware PAM filter on raw targets ✅
`PostProcess/pam_filter.py` filters raw CRISPRitz targets to those whose PAM region
can satisfy the requested PAM, **IUPAC-aware** so variant-*created* PAMs are kept
(e.g. an allele that turns `AGA`→`AGG` for NGG). This makes results correct
regardless of which index served the search, and is a no-op for PAM-specific and
direct (`-r`) searches.

*Why it matters:* the full pipeline never re-filtered reference targets by PAM, so a
degenerate/pamless index leaked non-matching PAMs into the results. Now fixed.

### 1.2 Pamless (NNN) index fallback ✅
`submit_job` now falls back to a degenerate **`NNN_<b>_<ref>[+<vcf>]`** index (for
both the reference and the variant/enriched genome) when a PAM-specific index is
absent. Combined with 1.1, **one pamless index can serve any PAM correctly**.
Validated end-to-end: pamless + filter reproduces the per-PAM (NGG) result exactly,
including the `PAM_creation` column and variant-sample attribution.

### 1.3 High-variant-density cap + greedy representative + report ✅
`new_simple_analysis.py` guards the IUPAC decomposition (the historical
`# the time of the universe` loop). A protospacer window overlapping more than
`CRISPRME_IUPAC_CAP` (default **10**) ambiguity codes would otherwise enumerate up
to 2^k haplotype sequences (dense/hypervariable regions like MHC/HLA; unbounded for
unphased VCFs). Above the cap it reports a **single greedy representative** per
haplotype instead of the full 2^k: at each ambiguity position it picks the
min-mismatch allele (preferring an alt on ties, which keeps PAM-region variant
alleles so **PAM-creation is preserved**). Because mismatch is additive per position,
this greedy target is provably the exact **min-mismatch / max-CFD** haplotype — the
single most concerning off-target — without enumerating anything. Validated:
4000/4000 random cases (both strands, with/without bulges) greedy == brute-force
argmin; and on a PAM-creating variant the greedy reproduces the full-enumeration
result exactly (creation + sample attribution identical). The region is logged — with
the **union of alt-allele-carrying samples** — to a new results file:

```
<results>/...high_variant_density_regions.bed
#chrom  start  end  guide  n_variants  samples_with_alt
```

*Naming note:* deliberately **not** "low-complexity" (that means repetitive/homopolymer
in genomics — the opposite of this). These are **high variant density** regions.

*Phased vs unphased:* phased VCFs (1000G, HGDP) prune impossible allele combinations
via sample-set intersection; unphased VCFs (e.g. gnomAD) cannot and are the worst
case. The cap protects both and is unavoidable for unphased data — so it is a
region-based mitigation, not an input rejection.

---

## 2. Performance — `--max-total-edits` enforced inside the search ✅

CRISPRme issue **#107**. Independent `--mm/--bDNA/--bRNA` limits combine additively
(6 + 2 + 2 = **10** edits), which explodes the search space in variant-aware /
pamless runs. `--max-total-edits` caps the **total** edits per alignment.

**What changed vs 2.1.x:** the cap is now **pruned inside the CRISPRitz TST search**
(`searchOnTST.cpp`, new optional `--max-edits` arg) so over-budget alignments are
**never generated** — not merely dropped afterward. A post-search `awk` drop remains
as a backstop for the `-r`/brute-force path and older binaries.

**Measured (chr22, one guide, pamless NNN index, 16 threads):**

| total-edits cap | reference search | rows |
|---|---|---|
| unbounded (2.1.x behavior) | **intractable** (>50 min, never completes) | — |
| 6 | **too heavy** (>10 min on chr22 alone, killed) → hours genome-wide | — |
| **5 (new default)** | ~88 s | ~24k |
| 4 | ~19 s | ~1.1k |

**Verdict:** `5` is the sweet spot — tractable and biologically ample; `6` is
already ~7× heavier on one chromosome and does not scale genome-wide; `4` is very
fast if a run needs to be lighter. Hence the default of **5**.

The in-search prune is **exact**: output equals the unbounded search filtered to
`Total ≤ cap` (verified byte-for-byte, zero false adds/drops).

**End-to-end (chr22 + 1000G, pamless NNN, 6/2/2, default max-edits 5 + density cap):**
the search that in 2024 never completed (2.3 TB, 392 GB/chr intermediate, 0 results)
now finishes in **16.5 min → 81 MB, 4,477 integrated results**, clean (no error log).

*Nuance — the cap bounds the SEARCH, not the post-decomposition final count.* On the
enriched genome an ambiguity code counts as a match if any allele matches the guide,
so the search enforces `≤ cap` on that basis. Variant **decomposition** then resolves
each IUPAC code to a specific allele; for the alt allele a previously-matched position
can become a mismatch, so a small tail of final variant targets legitimately shows
`>cap` total edits (observed: ~6% of rows at 6–9 edits). This is correct — those are
real variant off-targets — and is why the summary histogram had to be resized (see
the `process_summaries.py` fix below).

Plumbing: `crisprme.py --max-total-edits` → `submit_job` (26th positional) →
`crispritz.py search --max-edits` → `searchTST` argv[11] (also threaded to the indels
search via `pool_search_indels.py`).

### 2.1 Summary-histogram overflow fix (`process_summaries.py`) ✅
Pre-existing latent bug: the per-guide edit-distribution histogram was hard-sized
`[10][bulge+1]` but indexed by `Bulge_Size + Total` (up to `mms + 2*bulge`, e.g. 14
for 6/2/2) and `Total` (up to `mms + bulge`). It never surfaced on PAM-constrained
searches (few high-edit targets) but pamless / dense-variant runs raise
`IndexError: list index out of range` and abort summary processing. The array is now
sized to fit the index expression (`_new_distributions()`), so high-edit variant
targets (see the decomposition nuance above) no longer crash the run.

---

## 3. Web app — graphical Settings / Data Manager ✅

A browser-only, non-expert path to add data (local mode only; disabled on the public
`--website` server):
- **Genomes** — server-side fetch from UCSC (by assembly, e.g. `susScr11`) or
  HuggingFace, register a path, or chunked/resumable browser upload.
- **Indexes** — HuggingFace download or local build from an installed genome+PAM+bulges.
- **VCF datasets** — HF/URL/path/chunked upload, with reference-genome tagging so a
  VCF is only offered for its matching reference.
- **Annotations / PAMs** — upload/compose.
- Background jobs with progress, staging + atomic promote, disk pre-flight, input
  validation; publish-to-HF is maintainer-gated.

Plus web fixes: HOST/PORT binding bug, dead graphical-report wiring, `/load`
empty-state, stuck download banner, friendly PAM/Cas labels, genome-aware variant
defaults, `4/1/1` default search parameters.

---

## 4. Data / distribution

- **Compressed VCFs** used and stored as bgzip (`.vcf.gz` + `.tbi`), per-chromosome. ✅
- **HuggingFace footprint reduced**: reference genomes and annotations re-uploaded
  compressed (genomes 3.27→0.98 GB, annotations 1.65→0.11 GB), HGDP compressed. ✅
- **Default shipped index** (pamless NNN, combined 1000G+HGDP) — 🚧 pending the
  feasibility conclusion (see §2 + the max-edits verdict) before build + host.

---

## 5. Branding ✅
New CRISPRme+ logo (`assets/crisprme-logo.svg`, embedded Outfit font, transparent),
white navbar, favicons, vendored font, README + website typography.

---

## 6. Coordinated with CRISPRitz

- **PR #36 "enricher AF/FILTER robustness"** (Ann): `enricher.py`/`enricher.cpp`
  accept `FILTER='.'` (not only `PASS`, e.g. HPRC vcfwave output) and match the `AF`
  INFO key exactly (not a 2-char prefix that also hit `AF_afr=` etc.), with a
  multiallelic AF-count guard. Touches **only** the enricher — **provably disjoint
  from** our `searchOnTST.cpp` `--max-edits` change (`CRISPR-Cas-Tree/`), so the two
  compose with **zero conflict**. ✅ Integrated into our combined crispritz-src
  working tree (both changes coexist, verified); merge cleanliness confirmed.
  ⏸ Ann's PR still open upstream — bump the crispritz pin once it merges.
- **CRISPRme validator mirror** (`validate_inputs.py`): `ENRICHER_PASS_VALUES` now
  `("PASS", ".")` to match #36. Forward-compatible: on an older crispritz that only
  accepts `PASS`, a `.` record is simply enriched-out — the validator's superset is a
  warning, never a false pass.

---

## BREAKING CHANGES (2.1.x → 2.2.0 / CRISPRme+)

> The 2.2.0 lineage intentionally diverges and **may not be retro-compatible**. 2.1.x
> remains the last stable release. Breaking changes are listed here explicitly.

1. **`--max-total-edits` now defaults to `5`** (was effectively unlimited =
   `mm+bDNA+bRNA`). A run of `--mm 6 --bDNA 2 --bRNA 2` that previously reported
   alignments with up to 10 total edits will now **omit** alignments with >5 total
   edits. These are biologically implausible off-targets (many mismatches *and*
   several bulges), but the output set changes.
   **Escape hatch (restore 2.1.x behavior):** pass `--max-total-edits <mm+bDNA+bRNA>`
   (e.g. `--max-total-edits 10`) to effectively disable the cap.

2. **New results file** `..._high_variant_density_regions.bed` is written per
   variant search. Targets in windows above `CRISPRME_IUPAC_CAP` (default 10)
   ambiguity codes are represented by a bounded/greedy row rather than the full 2^k
   enumeration, so counts in those extreme regions differ from an uncapped 2.1.x run.
   **Escape hatch:** set `CRISPRME_IUPAC_CAP=-1` to disable the cap (full enumeration).

3. **Pamless (NNN) index fallback** changes which index may serve a search when a
   PAM-specific index is absent (previously an error / silent rebuild). Results are
   equivalent (PAM enforced by the post-search filter), but the index-selection
   behavior differs.

### NOT broken (verified retro-compatible)
- Existing CLI invocations run unchanged (no args removed; new args optional).
- Cached indexes still work (`searchTST`'s `--max-edits` is an optional 11th arg;
  omitted ⇒ disabled ⇒ 2.1.x behavior; direct `crispritz.py search` without
  `--max-edits` is unchanged).
- Package/CLI name and version (`crisprme`, `2.2.0`) unchanged for bioconda/Docker.

---

---

## 7. 2.1.x end-of-life policy

2.1.x is the **last stable line** and is being retired. Decision:

- **One final `2.1.14` patch**, then hard-freeze/archive. It backports **only the
  non-feature correctness fixes** that could affect users mid-analysis on 2.1.x:
  1. **stderr-is-fatal** — `[ -s $logerror ]` checks fail any stage whose subprocess
     merely writes to stderr (hit bit post-analysis in 2.1.13).
  2. **input-validator hardening**.
- **No feature backports.** `--max-edits`, the density cap, pamless fallback,
  Settings/Data Manager, compressed HF, branding, and the CRISPRitz changes are
  **2.2.0-only** — they depend on the diverging pipeline and are not worth re-basing
  onto a retiring line.
- After 2.1.14: all development on 2.2.0 / CRISPRme+; 2.1.x archived on the old repo.

---

*Last updated during 2.2.0 development. Flip 🔬/🚧/⏸ items as they are validated and
integrated.*
