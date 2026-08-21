# CRISPRme / CRISPRitz — project status & open threads

Living tracker so nothing gets lost. Forward plan lives in
[`ROADMAP.md`](ROADMAP.md); measured findings in
[`PERF_TIMING_2112.md`](PERF_TIMING_2112.md) and
[`DOC_AUDIT_2112.md`](DOC_AUDIT_2112.md).

_Last updated: 2026-08-03._

## Release lines (branch isolation)

| Line | Branch | Purpose |
|---|---|---|
| CRISPRme released | `main` | 2.1.12 (shipped: GitHub tag + Bioconda + multi-arch Docker) |
| CRISPRme 2.1.13 | `2.1.13-dev` | patch: robustness + memory cap |
| CRISPRme+ v2.4.0 | `main` | current stable: Ann's search, Python 3.11 + Dash 2.x; dict-less 1000G+HGDP variant index is the **default** distribution; shareable report, COSMIC annotation, high-variant-density flagging |
| CRISPRitz released | `master` | 2.7.1 |
| CRISPRitz 2.8.0 | `2.8.0-dev` | major: Python 3.11, C++ enricher, parallelism |

## Open PRs

**CRISPRme**
- #137 → `main` — CRISPRme+ roadmap (this plan)
- #138 → `main` — companion-guide doc fixes
- #136 → `2.1.13-dev` — download retry/resume + post-analysis 64 GB memory cap
- #113 → `2.2.0-dev` — Ann's assembly-search (diploid) [@anncir1]
- #131 → `2.2.0-dev` — Python 3.11 + Dash 1.x→2.x web migration
- _(merged: #135 — README + website)_

**CRISPRitz** (all → `2.8.0-dev`)
- #26 — byte-identical C++ enricher (replaces single-threaded Python enricher.py)
- #27 — memory-bounded cross-chromosome parallelism + use compiled enricher
- #24 — Python 3.11 integrated modernization + Docker + CI
- #21 — modernize: Python 3.11 + drop stale bytecode
- #23 → `modernize-python` — azimuth on-target scoring via CRISPR-HAWK (stacked on #21)

## Key measured findings (genome-wide 1000G e2e, 4/1/1)

- Total ~13.5 h; **~93 % (12.5 h) is variant enrichment** (was single-threaded Python).
  → fixed by CRISPRitz #26 (C++) + #27 (parallel).
- **Peak RAM ~98 GB in post-analysis** (not the search). → capped to 64 GB (CRISPRme #136).
- Download slow/brittle: EBI ~1.4 MB/s, no retry. **HF ~6.5× faster.** → #136 retry/resume + HF plan.

## Open threads / TODO (CRISPRme+)

- [ ] **HF data hosting** — repo `lucapinello/crisprme-data` (private, dataset card up).
      Upload the complete-test reference data (genome + 1000G VCFs + annotations +
      PAMs + samplesIDs) from ml007. _(in progress)_
- [ ] **Precomputed indexes** on HF for the default references (skip download+enrich+index).
- [ ] **Default references**: 1000G+HGDP+TOPMed+All of Us union; Pangenome 2.0; NNN+NGG PAMs.
- [ ] **Model organisms**: add pig (Sus scrofa) + mouse (Mus musculus).
- [ ] Transfer the HF dataset repo to a `pinellolab` org (currently under `lucapinello`).
- [ ] Wire CRISPRme `setup` to pull from HF instead of EBI/UCSC.
- [ ] Manuel/Ann: review the 2.8.0 / 2.2.0 stacks; divide work; manuscript.

## Deferred (from 2.1.13 roadmap §C/D — see ROADMAP.md)

- #112 input validator (fix samplesID false-positive), #107 max-edits filter,
  #108 deep scalability streaming, separate DNA/RNA bulge columns, pickled-model
  modernization (CRISTA/azimuth sklearn shims).
