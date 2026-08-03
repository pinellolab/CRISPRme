# CRISPRme 2.1.12 — Per-step timing (from-scratch full-genome run on ml007)

Purpose: measure where wall-clock time goes in a genuine new-user full-genome
`complete-test --vcf_dataset 1000G` (4/1/1), to target CRISPRme+ optimizations.

**Machine:** ml007 — 256 cores, 1 TB RAM, `/srv/local` SSD. Run launched 19:16:34.
**Method:** download phases timed via file mtimes; compute phases will be read from
CRISPRme's `Results/*/log.txt` (authoritative `Stage<TAB>Start/End` markers).

## Measured so far

| Stage | Start | End | Duration | Size | Notes |
|---|---|---|---|---|---|
| Setup + **genome download+extract** + annotations/samples | 19:16:34 | ~19:40:13 | **~24 min** | 3.1 GB genome | hg38.chromFa.tar.gz + extract all contigs |
| **1000G VCF download** (phase3 GRCh38) | 19:40:13 | *ongoing* | **>97 min** (incomplete) | 8.4 GB / 10 of 23 chroms | **~1.4 MB/s — the dominant bottleneck** |
| VCF preprocessing (bcftools) | — | — | TBD | — | not started |
| Genome index build (CRISPRitz TST) | — | — | TBD | — | not started |
| Variant/enriched search (CRISPRitz 2.7.1) | — | — | TBD | — | 256 cores available |
| Post-processing (merge/annotate/score CFD+CRISTA) | — | — | TBD | — | candidate for Rust rewrite |
| Index-genome Reference (buildTST, C++) | 10:30:32 | 10:39:51 | **~9 min** | — | fast, C++ |
| Index-genome Variant (enriched) | 10:39:51 | 10:49:49 | **~10 min** | — | fast, C++ |
| Indexing Indels | 10:49:49 | 10:49:54 | **~5 s** | — | trivial |
| Off-targets search (searchTST, C++) | 10:49:54 | *running* | TBD | — | multi-threaded (~16 cores); **peak RAM 55.6 GB** |
| Post-processing (merge/annotate/score) | — | — | TBD | — | Rust candidate |
| Report + image generation | — | — | TBD | — | — |

**🔴 Peak RAM finding:** the search stage peaked at **55.6 GB** (genome-wide 1000G) —
**above the README's stated 32 GB minimum.** The README's "64 GB recommended" is the
real practical floor for whole-genome variant searches; the 32 GB "minimum" should be
reworded (doc fix + scalability note). Contrast: enrichment/index stages stayed ≤13 GB.

## Findings / optimization targets (CRISPRme+)

1. **Download dominates — by far.** VCF download alone is projected at **~3+ hours**
   (~16 GB 1000G phase3 at ~1.4 MB/s from the current source). This is the #1
   new-user pain point and pure network cost, not compute.
   → **Mirror reference data on Hugging Face / a fast CDN**, and/or
   → **Ship precomputed indexes** for the default references (1000G+HGDP+TOPMed+AllofUs
      union; Pangenome 2.0) so users skip download+preprocess+index entirely.
2. **Genome extract preserves 2014 file mtimes** — cosmetic, but means extraction
   time must be measured from process logs, not mtimes.
3. Compute-stage timings (preprocess / index / search / **post-process** / plots)
   to be filled from `log.txt` once the search runs — post-processing is the
   expected Rust-rewrite target; this run will quantify its share.

## Source benchmark (measured on ml007, 15s sustained, concurrent with the live run)

| Source | Throughput | Notes |
|---|---|---|
| EBI 1000G (`ftp.1000genomes.ebi.ac.uk`) — current VCF source | **1.4 MB/s** | trans-Atlantic (EMBL-EBI), rate-limited |
| **Hugging Face** (Cloudfront CDN) | **9.2 MB/s** | **~6.5× faster**; conservative (test ran under load) |
| UCSC (`hgdownload.soe.ucsc.edu`) — current genome source | 1.4 MB/s | also throttled — mirror this too |

- ml007 sustained ~10 MB/s aggregate during the test → EBI/UCSC are the bottleneck, not ml007's link.
- 16 GB 1000G: **~3.2 h (EBI) → ~29 min (HF)**. HF likely faster on an idle link.
- Caveat: HF fair-use bandwidth policy; fine for dataset hosting at this scale (LFS).
- Bigger win: **precomputed indexes on HF** skip download+preprocess+index-build entirely.

## Live finding: the "Add-variants" stage is the compute bottleneck

- **What it is:** `crispritz.py add-variants` → **`enricher.py`** (pure Python,
  `#!/usr/bin/env python`, in **CRISPRitz** `Python_Scripts/Enrichment/`). Builds the
  variant-enriched genome from the VCFs, **per chromosome**.
- **Observed (measured):** genome-wide enrichment ran **~12h30m** (Add-variants
  Start Aug 2 22:00:45 → End Aug 3 10:30:32). The active
  `enricher.py` sits at **~101% CPU = single-threaded** (1 of 256 cores), ~5 GB RSS,
  processing **one chromosome at a time — it ignores `--thread`**. This single stage
  dwarfs everything else and makes a genuine from-scratch genome-wide 1000G run
  impractical (~half a day before the search even starts).
- **Optimization (CRISPRitz work item, not CRISPRme):**
  1. **Parallelize across chromosomes** (embarrassingly parallel — 23 chroms → ~10–20×
     with minimal change). Biggest, cheapest win.
  2. **Rust rewrite of `enricher.py`** for per-variant speed on top of that.
  → Pair with the CRISPRitz 2.8.0 / CRISPRme+ track.

## FINAL — genome-wide run completed (Job Done, EXIT=0)

**Total wall-clock: ~13h25m** (Aug 2 22:00:44 → Aug 3 11:25:35).

| Stage | Duration | Share |
|---|---|---|
| **Add-variants (enrichment, single-thread Python `enricher.py`)** | **12h30m** | **~93%** |
| Index-genome Reference (C++ buildTST) | 9m19s | |
| Index-genome Variant (C++) | 9m58s | |
| Indexing Indels | 5s | |
| **Off-targets search (C++ searchTST, ~16 cores)** | 28m56s | ~4% |
| Post-analysis SNPs | 2m51s | |
| Post-analysis INDELs | 25s | |
| Merging Targets | 18s | |
| Annotating results | 1m3s | |
| Creating images | 58s | |
| Integrating results | 1m5s | |
| Creating database | 3s | |

**Peak RAM: 97.7 GB** — and it's the **`Post-analysis SNPs`** stage, NOT the search.
Memory curve: enrichment/index ≤13 GB; search ~30 GB steady (56 GB peak at the
reference→variant index handoff, 11:05); then `Post-analysis SNPs` fanned out to
**~161 parallel processes** and their aggregate RSS spiked to **97.7 GB at 11:20:06**
(stage lasted only ~3 min), then collapsed to ~0.
- **Memory bottleneck = unbounded parallelism in Post-analysis SNPs** (≈160 workers,
  each loading result data). Distinct from the *time* bottleneck (enrichment).
- **Fix (cheap, high-impact):** bound the post-analysis worker pool to a memory
  budget / core count. Stage is only ~3 min, so throttling costs little time but could
  drop peak RAM ~98 GB → ~32–64 GB, making genome-wide fit on normal hardware.
- Caveat: 97.7 GB is summed RSS across 161 procs (over-counts shared pages); true
  unique memory is lower, but the fan-out spike is real.

**Outputs (correct + complete):** 55,719 integrated off-targets; 171,130 alt-alignments;
6 summary files (guide+samples × CFD/CRISTA/fewest); 23 plots in imgs/; best/altMerge; db.

### Bottom line for CRISPRme+ optimization
- **93% of wall-clock is one single-threaded Python stage (`enricher.py`, CRISPRitz).**
  Parallelize across chromosomes (~10–20×, near-free) + Rust → ~13.5 h collapses toward ~1 h.
- **Download (~3 h EBI) + enrichment (~12.5 h)** are the whole cost; the actual search +
  scoring + plots total ~40 min. → HF mirror + **precomputed indexes** removes almost all of it.
- Whole-genome 1000G needs **~64–100 GB RAM**, not 32 GB (doc fix S1).
