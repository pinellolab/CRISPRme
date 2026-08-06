# Changelog

All notable changes to CRISPRme are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crisprme` skill.

## [Unreleased]

## [2.2.0] - 2026-08-06

### Added
- Graphical **Settings / Data Manager** page in the web interface: add reference
  genomes (UCSC by assembly name — e.g. the pig `susScr11` — HuggingFace, or a
  direct URL), precomputed indexes (download from HuggingFace or build locally
  from an installed genome + PAM), VCF datasets (HuggingFace or register an
  existing server folder), annotations (BED upload), and nucleases/PAMs (a small
  form). New data lands in the local data folder and is auto-discovered by the
  search form. Long operations run as detached jobs on a dedicated executor with
  live progress, so they never starve the search slots. Mutations are local-mode
  only; publishing an index to the shared HuggingFace repo is maintainer-only.
  `download --what genome` gained `--source {hf,ucsc,url}` (+ `--url`) so the CLI
  and web share one non-human-genome download path.
- Python 3.11 modernization: pipeline fixes for pandas 2.x / matplotlib 3.x and
  a Dash 1.x → 2.x web-app migration, plus a Python-3.11 Docker image built from
  source (CRISPRitz 2.8.0) (#131).
- `assembly-search` subcommand: off-target search on a personal diploid genome
  assembly (two haplotypes, no VCF), reconciled to hg38 via liftOver (#113).
- Reference-index UX: `build-index-only` pre-builds the reusable CRISPRitz
  reference index without running a search, and `complete-search --index-path`
  reuses a prebuilt/staged index library (a missing index is a hard error rather
  than a silent rebuild).
- HuggingFace data distribution: `download` fetches reference data (genome,
  annotations, PAMs, sample IDs, VCFs, precomputed indexes) from a HuggingFace
  dataset repository over its CDN, and `publish-index` uploads a locally built
  index for reuse. Default repo `lucapinello/crisprme-data`, overridable via
  `--hf-repo` / `CRISPRME_HF_REPO`. `setup`/`complete-test` also try HuggingFace
  first and fall back transparently to the original UCSC/EBI/Sanger sources
  (#140, #141).
- `complete-search --max-total-edits N`: cap the total edits (mismatches +
  bulges) per reported alignment; over-cap targets are dropped right after the
  search, shrinking intermediate files and post-analysis time (#107).

### CI
- New `unit tests` workflow: fast, hermetic byte-compile + network-free HF/index
  unit tests on every code PR.
- New `web e2e (playwright)` workflow: builds the py3.11 image, serves the web
  app, and drives Chromium to assert every Dash 2.x page renders (no blank pages
  / JS errors).
- `validate-benchmarks` gained a new-subcommand dispatch + unit-test smoke step.

### Changed
- Clearer failure reporting: when a search fails, CRISPRme now prints *which
  stage* failed (from the per-stage log) and the last lines of the error log,
  instead of only "run failed — see log_error.txt". Makes failures actionable
  for non-expert users.

### Fixed
- Web interface (Dash 2.x) hardening, from a full Playwright stress test of the
  running app:
  - The web server no longer crashes on a from-source install. Dash 2.x's
    `app.run()` lets the `HOST` environment variable override the host argument,
    and the from-source conda env sets `HOST` to a non-bindable compiler build
    triple (`x86_64-conda-linux-gnu`); the server now forces the intended host/
    port so it binds correctly.
  - The **Query Genomic Region** and **Personal Risk Cards** result tabs no
    longer return HTTP 500. Both callbacks type-checked their inputs before the
    "no click yet" guard, so Dash's initial (empty) render raised a `TypeError`;
    the guard now runs first, and Filter/Generate with nothing selected is a
    graceful no-op.
  - Removed the dead cross-origin "skeleton" stylesheet (blocked by browsers on
    every page); the layout already uses the Bootstrap grid.
  - The nuclease dropdown collapses case-variant duplicate PAM files so each
    nuclease is listed once.
- Zero-hit searches now complete cleanly with an empty result instead of
  aborting. A search that finds no off-targets (e.g. a very stringent
  guide/parameter combination) previously failed part-way through post-analysis
  ("off-targets post-analysis (reference) failed", then a cascade through the
  rsID / summary / integration steps, all of which assumed at least one target).
  Added a zero-target guard to the reference SNP post-analysis (mirroring the
  INDELs one), made `remove_n_and_dots.py` tolerate a header-only report, and
  added a high-level short-circuit that emits an empty-but-valid result and exits
  0 when no off-targets are found. Verified end-to-end on ml007 (empty result
  exits 0; a normal with-hits search is unaffected).

### Documentation
- New `docs/DOCKER_QUICKSTART.md` + a README quickstart callout: a few-command,
  no-conda / no-410 GB path to the web interface for non-experts — fast HF data
  download, a prebuilt index, then `web-interface` in the browser, with
  troubleshooting and "install more indexes as needed".
- Data-setup guide: documented what `setup` produces (including the combined
  1000G+HGDP config files), the HuggingFace fast-download path, and a new
  "Prebuild, reuse, and share the reference index" section.
- README: added a from-source install path and a reference-index /
  data-distribution commands section.

## [2.1.13] - 2026-08-04

### Added
- Robust `complete-test` downloads: HTTP retry with linear backoff plus
  checksum-verified resume, so a dropped connection or truncated transfer from a
  slow/flaky host retries (and skips already-verified files) instead of aborting
  a multi-hour run (#136).
- Memory-bounded post-analysis: per-chromosome SNP/indel workers are capped to a
  RAM budget (default 64 GB, overridable via `CRISPRME_MAX_MEM_GB` and
  `CRISPRME_POSTPROC_WORKER_GB`), preventing the peak-RAM spike observed on
  genome-wide 1000G runs (#136).

### Fixed
- Post-analysis worker-count diagnostic is now written to stdout (`log_verbose`)
  instead of stderr. The caller treats a non-empty stderr log as a fatal
  post-analysis failure, so the informational message previously aborted every
  run with a false "post-analysis (snps) failed" even when the analysis
  succeeded (#136).

## [2.1.12] - 2026-08-01

### Added
- Cas9 and Cas12a validation benchmark examples with precomputed brute-force
  ground-truth references and an extensible registry
  (`test/benchmark/benchmarks.json`); `validate-test` compares CRISPRme
  off-targets against them (#116).
- A fast, dependency-free Rust port of the brute-force ground-truth generator
  (`test/benchmark/rust/`) that produces output identical to the Python
  generator (#116).
- `validate-benchmarks` GitHub Actions CI that runs the full
  `complete-test` → `validate-test` round-trip on Linux for every registered
  benchmark (#116, #120).
- Native Apple Silicon (`linux/arm64`) Docker image via a multi-arch `buildx`
  workflow; the multi-arch image is published to Docker Hub on release tags
  (#121).
- Low-memory startup warning (`PostProcess/memory_check.py`) advising Docker
  Desktop users to raise the memory limit when total RAM is below ~32 GB (#121).
- Search-space budget estimator and guard (`PostProcess/search_budget.py`) that
  warns before resource-explosive runs, with a companion analysis in
  `docs/SCALABILITY_ANALYSIS.md` (#118).
- Release tooling: a `release-crisprme` Claude Code skill, `docs/RELEASING.md`,
  this `CHANGELOG.md`, and `scripts/prepare_release.py` (#117); plus docs for the
  fast Rust generator and using the release skill (#122).
- Input-format hardening: `--gene_annotation` is now auto-sorted like
  `--annotation`; clear errors for multi-line PAM files and malformed PAM
  filenames; and a non-fatal warning when a degenerate PAM motif is combined with
  bulges (the CRISPRitz #105 crash, fixed in CRISPRitz 2.7.1) (#125).
- Documentation: `docs/INPUT_FORMATS.md` (PAM/Cas12a/VCF/chromosome/annotation/
  bulge guidance, #124), `docs/ISSUE_AUDIT.md` (a review of all open + closed
  issues, #126), and `docs/ROADMAP.md` (the 2.1.12 → 2.1.13 plan, #128).

### Fixed
- Multiple bugs affecting new VCF dataset processing, concurrent runs, and
  gnomAD handling (#96).
- Guard the CFD PAM-score lookup against non-ATCG PAM bases present in real hg38,
  which previously crashed the run with a `KeyError` (#94, #125).
- Guard the post-analysis intermediate reads against non-UTF-8 bytes
  (`UnicodeDecodeError`) so a stray byte no longer aborts a run (#52, #125).
- `validate-test` now exits non-zero when a benchmark mismatches or the
  `complete-test` output is missing, instead of passing silently (#116).
- `complete-test` downloads 1000 Genomes VCFs over HTTPS instead of FTP,
  unblocking FTP-restricted and CI networks (#116).

### Changed
- Synced the `Dockerfile` `crisprme_version` pin to the released version, which
  had lagged behind `crisprme.py` and the Bioconda recipe (#117).

## [2.1.11] - 2026-07-06

### Fixed
- Consistent `.gz` suffix for compressed annotation files, removing a mismatch
  between file content and file name in the standard and personal-annotation
  pipelines (`crisprme.py`).
- `_check_personal_annotation` now strips a `.gz` suffix before processing so
  downstream steps see the expected filename.
- Removed double-compression of annotation files, preventing corrupted or
  unreadable annotation outputs.

### Changed
- Refactored per-chromosome test genome-directory handling: replaced
  `ensure_hg38_directory` with `_assign_genome_directory_name` (whole-genome
  builds in `<dest>/hg38`, single-chromosome builds in `<dest>/hg38_<chrom>`).
- `download_genome_data` now returns the resolved genome directory path directly.
- Renamed downloaded sample-ID files to include genome/dataset context
  (`hg38_{dataset}.samplesID.txt`) via a centralized mapping.
- Reduced default web-interface example-run parameters (max mismatches 6 → 4,
  DNA bulges 2 → 1, RNA bulges 2 → 1) for a faster example run.

### Added
- New `_write_sg1617_file` helper to generate an `sg1617` guide fixture during
  legacy database setup.
- New guides `docs/crisprme_data_setup_051826.md` (CLI setup/usage) and
  `docs/crisprme_web_interface_user_guide.md` (web interface walkthrough).
- "Setup Legacy Database" section in `README.md`.

### Removed
- Unused `_bgzip_ann_data` helper and stale typing imports/docstrings.

## [2.1.10] - 2026-05-29

### Added
- New `setup` command for automated download/installation of CRISPRme reference
  resources directly from the command line, removing the need for manual
  downloads from the CRISPRme website.

### Changed
- Updated documentation and container configuration to support the new
  installation workflow and improve reproducibility across deployments.

## [2.1.9] - 2026-01-16

### Added
- `validate-test` functionality: validates `complete-test` off-target
  predictions against brute-force ground-truth alignments derived from
  1000 Genomes variant data (PR #92).
- Chromosome-level validation support (single chromosome or genome-wide).

### Changed
- Strengthened the `complete-test` pipeline and refined benchmarking utilities
  for correctness and reproducibility with population-scale variant data.

## [2.1.8] - 2025-12-10

### Fixed
- Targeted fixes and stability improvements across annotation handling,
  off-target search execution, test coverage, and report generation.
- Fixed and updated the `complete-test` workflow to align with the latest
  internal logic and outputs (PR #86).

### Changed
- Updated Docker distribution and testing environment for consistent runtime
  behavior (PR #85).

## [2.1.7] - 2025-03-10

### Added
- Support for exome-based gnomAD VCF files.

### Fixed
- Corrected handling of indels in off-target reporting (previously missing or
  misreported).
- Fixed the visualization of the Personal Card tab in the GUI and website.
- Closed file streams after creating graphical reports to prevent memory issues.

### Changed
- Improved contiguous-target merge logic (always reports the leftmost off-target
  in a cluster).
- Improved error-tracing logic.

## [2.1.6] - 2024-11-27

### Added
- Support for gnomAD v4.1 VCFs, including joint variant files.
- New `complete-test` function for quick single-chromosome and full-genome
  off-target testing.

### Fixed
- Corrected indel handling in off-target reporting to eliminate false positives.
- Fixed the alternative-alignments report to list all possible alignments.
- Resolved Personal Card tab visualization issues in the GUI and website.
- Addressed a memory-overflow issue when merging contiguous targets.

### Changed
- Upgraded the DockerHub image with the latest fixes.

[Unreleased]: https://github.com/pinellolab/CRISPRme/compare/v2.1.13...HEAD
[2.1.13]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.13
[2.1.12]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.12
[2.1.11]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.11
[2.1.10]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.10
[2.1.9]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.9
[2.1.8]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.8
[2.1.7]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.7
[2.1.6]: https://github.com/pinellolab/CRISPRme/releases/tag/v2.1.6
