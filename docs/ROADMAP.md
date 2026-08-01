# CRISPRme / CRISPRitz roadmap

This document records the release strategy: a **stable 2.1.12** on the current
Python (so the Dash web app keeps working) with every known bug fixed, followed
by a **2.1.13** modernization + new-features release (Python 3.11, Dash upgrade,
dependency bumps, deferred features).

---

## 2.1.12 — stability release (current, Python 3.8)

**Goal:** a stable, fully-tested release with all known bugs fixed. Deliberately
**stays on Python 3.8** so the Dash web interface and its dependency stack keep
working unchanged — **no modernization risk**.

**Included (merged to `main`):**
- VCF / gnomAD / concurrent-run fixes (#96); annotation `.gz` handling (#109);
  malformed-cluster-row guard (#114).
- Cas9 + Cas12a validation benchmarks + brute-force ground truth + Rust generator,
  and the `validate-benchmarks` CI (#116, #120).
- Scalability guard `search_budget.py` + analysis (#118).
- Native Apple-Silicon (arm64) multi-arch Docker image + low-memory startup
  warning (#121).
- Input hardening (#125): non-ATCG PAM guard (**#94**), `--gene_annotation`
  sorting, multi-line-PAM + PAM-filename validation, a `UnicodeDecode` guard
  (**#52**), and a **#105** guard (clear warning on degenerate-PAM + bulge).
- Docs: `INPUT_FORMATS.md`, `ISSUE_AUDIT.md`, `RELEASING.md`, release skill.
- Full audit of all 8 open + 58 closed issues — nothing closed by mistake.

**CRISPRitz side (2.7.1):** the #105 heap-overflow fix (PR #19), a CI safety net
(PR #20), and the stale-`.pyc` cleanup (#17). CRISPRme 2.1.12 runs on the existing
**crispritz 2.7.0** with the CRISPRme-side #105 guard, so it is not blocked on the
crispritz/Bioconda chain. Repinning to crispritz 2.7.1 can fold in at tag time if
Bioconda is ready, else it moves to 2.1.13.

**Explicitly deferred out of 2.1.12 (→ 2.1.13):** the input validator (#112, has a
`samplesID` false-positive that aborts valid 1000G runs), Python 3.11, the Dash
upgrade, and the assembly-search feature (#113).

---

## 2.1.13 — modernization + new features (next)

### A. Python 3.11 migration
Both tools are currently pinned to Python 3.8 (EOL). Assessed difficulty: **MEDIUM**.
- **CRISPRitz → 2.8.0:** drop stale `.pyc` (#17, done in PR #21); bump the pins;
  **adopt CRISPR-HAWK's `scores/azimuth`** (a modern, sklearn-1.1-loadable model) to
  restore on-target Doench scoring — the only hard part, since CRISPRitz's own
  azimuth model is a sklearn-0.21 pickle; **fix the broken standalone Dockerfile**
  (fails today at `conda install python=3.8`). Groundwork: PRs #21, `azimuth-crisprhawk`.
  The off-target (CFD) pipeline already runs on 3.11 (tested).
- **CRISPRme → 2.2.0:** pipeline fixes — `DataFrame.append`→`pd.concat`
  (`process_summaries.py`), `seaborn-poster`→`seaborn-v0_8-poster` (4 sites),
  bare `python`→`sys.executable` (`generate_sample_card.py`, `crisprme.py`); a
  **CRISTA pickle forward-compat shim** (`sklearn.ensemble.forest`→`_forest`,
  `sklearn.tree.tree`→`_classes`) + pin `numpy<2`; repin to the crispritz-py3.11
  build. Groundwork: PR `py311-core`. **Prerequisite:** a crispritz **Python-3.11
  Bioconda build** (compiled; linux-64 + linux-aarch64) must exist first.

### B. Dash web-app modernization (the long pole)
The web app is written for **Dash 1.x** (`dash 1.10`, `flask 1.1`, `werkzeug 1.0`,
`itsdangerous<2.0`, and the retired `dash-core-components`/`dash-html-components`).
Migrate to **Dash 2.x**: rewrite 11 files' imports to `from dash import dcc, html,
Input, Output, State`, `app.run_server`→`app.run`, and bump Flask/Werkzeug/
itsdangerous/plotly. Not covered by the benchmark CI — **validate with Playwright**
click-through (as in the original audit).

### C. Deferred features / bugs
- **#112 input validator** — re-enable after fixing the `samplesID` sample-check
  false-positive (it currently `ERROR`s on the valid 1000G workflow). It's the
  FDA-requested pre-flight check.
- **#113 assembly-search** subcommand (re-test on current `main` + review).
- **#107** true `--max-edits` (total mm+bulges) filter.
- **#108** deep scalability fix (stream `best*` via `sort -m`; flush at cluster
  boundaries) — see `docs/SCALABILITY_ANALYSIS.md`.
- Separate **DNA/RNA bulge count** columns in the output.
- CRISPRitz **#13** (non-human/scaffold genome crash — check if the #105 fix
  resolves it) and the multi-motif / additional-PAM improvements.

### D. Pickled-model modernization (cross-cutting)
Both tools ship models pickled under sklearn ~0.21 (paths removed in 0.22): CRISTA
(CRISPRme) and azimuth (CRISPRitz). Fix via a `sys.modules` alias shim and/or
re-export, and **validate scores numerically** against the benchmark ground truth
(a clean unpickle is not sufficient).

---

## Order of operations (priority: stability first)

1. **Ship CRISPRme 2.1.12** — finalize docs `#124` + changelog `#123`, tag `v2.1.12`
   (triggers the arm64 Docker publish + Bioconda autobump). Old Python, working Dash,
   all bugs fixed.
2. **CRISPRitz 2.7.1** (parallel) — Manuel reviews/merges #19 (#105) + #20 (CI) +
   the #17 cleanup → tag `v2.7.1` → Bioconda → optionally repin CRISPRme.
3. **2.1.13 modernization** — Python 3.11 (CRISPRitz 2.8.0 → CRISPRme 2.2.0) + Dash
   2.x + deferred features, each gated by the benchmark CI (and Playwright for the web).

## Validation gates
- `validate-benchmarks` (byte-identical Cas9 + Cas12a off-targets, incl. CFD **and**
  CRISTA) gates every CRISPRme change.
- CRISPRitz CI (build + AddressSanitizer #105-regression + native arm64).
- Playwright click-through for the Dash web app.
- Multi-arch Docker build (linux/amd64 + linux/arm64) for both images.
