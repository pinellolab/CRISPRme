# CRISPRme issue audit (pre-2.1.12)

A full review of every open and closed GitHub issue, done to make the 2.1.12
release genuinely stable: confirm that closed issues were not closed by mistake
or since regressed, and that the release covers the open problem surface.

**Scope:** 8 open + 58 closed issues, audited against `main` (post #96/#109/#114/#116/#118/#120/#121).
**Headline:** No issue was closed by mistake. A handful are "closed but not code-provable
from this repo" (external recipe / server-side / user-closed) — none are wrong. One real
defensive gap (#52) and one watch item (indel false-positives) are folded into 2.1.12.

Legend for action: **MUST-FIX** = in 2.1.12 · **2.1.13** = deferred · **doc** = documentation ·
**close** = resolvable at/after release · **WAI** = works-as-intended.

## Open issues (8)

| Issue | Summary | Status on `main` | Action |
|---|---|---|---|
| #115 | Coordination meta-issue (ours) | Its listed work is done (#96/#109/#114 merged; #106/#107/#108 addressed) | **close** when release is cut |
| #111 | Does `--vcf` support SVs? | SVs unsupported — CRISPRme resolves SNV/indel to explicit REF/ALT only; no symbolic/BND handling | **doc** (state SNV/indel-only) |
| #110 | Request for benchmark sgRNAs + reference results | Satisfied by #116 (Cas9 `sg1617` + Cas12a HBG benchmarks, brute-force refs, `validate-test` CI) | **close** post-release |
| #106 | v2.1.10 failed with `sg1617.txt` | Actual crash was `cluster.dict.py:137` IndexError on a truncated intermediate row; now guarded by #114 (`cluster.dict.py:143`). Not the annotation/OOM bugs | **close** (robustness fix; root-cause scalability tracked in #108) |
| #94 | Fails on 1000G 30xNYGC dataset | Core failures fixed by #96 (chr-from-filename, haploid GT, pandas dtype). **Residual:** non-ATCG PAM `KeyError` at `new_simple_analysis.py:123` | **MUST-FIX** (guarded in #125), then **close** |
| #108 | complete-test chr1 OOM in `new_simple_analysis.py` | Mitigated by the `search_budget.py` advisory guard (#118) + `docs/SCALABILITY_ANALYSIS.md`; the in-pipeline OOM is unchanged | **2.1.13** (deep fix D2) |
| #107 | Filter by max total number of edits | No true in-search `--max-edits` filter; `search_budget.py --fail-over` only warns | **2.1.13** (feature) |
| #105 | Crash with bulge=1 on WTN PAM | CRISPRitz C++ heap overflow (`searchOnTST.cpp:1073`), triggered by odd-length degenerate PAM + bulge; not Python-fixable | CRISPRme **guard** in 2.1.12 (#125); C++ fix → **CRISPRitz 2.7.1** |

## Closed issues — verification (highlights)

All 58 verified. Correctly-closed-and-fixed (code confirmed in `main`) include: **#103** (tabix flock/mtime
reindex, `annotation.py:63-75`), **#98** (`bMax` auto-computed `crisprme.py:1179`; annotation optional), **#81/#54**
(rsID-in-AF cleaned in the converter + GT guards `resultIntegrator.py:144-146`), **#64/#73** (gnomAD v4.1 joint
+ exome handling in `convert_gnomAD_vcfs.py`), **#74/#68** (contiguous-merge leftmost/sample-priority), **#80**
(dropped `python` prefix in `complete_test.py`), **#72** (stderr→log_error tracing), **#77/#67/#59** (indel
false-positive filter `remove_bad_indel_targets.py` + #96 chr-assignment fix). Correctly-closed WAI: **#93**
(CFD is NGG-calibrated — a perfect NG on-target legitimately scores ~0.016), **#100** (Azimuth is dead code),
**#99** (index already reused), **#82** (`--debug` is source-only), **#76** (guide-N convention), **#95** (→ design tool).

### Flags — closed but not code-provable here (none are wrong)
- **#52** `UnicodeDecodeError` in `cluster.dict.py` — the only real *defensive gap*: intermediate reads used strict UTF-8 (#114 only guarded `IndexError`). **Fixed in #125** (`errors="replace"`).
- **#63** pysam missing from the conda recipe — fix lives in `bioconda-recipes/meta.yaml` (external); verify on the recipe side.
- **#45** Cas12a web "no results" — server-side outage, not repo code.
- **#84** missing InDel-file error — closed by the reporter; nominal cleanup remedy, no pinpointed guard. Low risk.

### Watch item
- **#67/#77/#59** indel false-positives — the coordinate-overlap filter (`remove_bad_indel_targets.py:62-79`) is
  present and active, and #96 fixed the indel chromosome-assignment bug (a major source), but the exact
  *sequence-level* false positive #77 described isn't provably eliminated by an overlap test alone. Recommend a
  targeted regression test on the #77 dataset.
- **#39** radar-chart `guide[count]` IndexError — legitimately closed (user error), but the code still lacks a
  guide-length guard; candidate for a future friendly-error check.

## 2.1.12 disposition summary
- **MUST-FIX (in flight):** #94 non-ATCG PAM guard + #52 decode guard + FDA input-hardening (PR #125); #105 CRISPRme guard (PR #125); blend the #112 input validator.
- **Deferred to 2.1.13:** #105 CRISPRitz C++ fix (→ CRISPRitz 2.7.1), #107 (max-edits filter), #108 (deep OOM fix), #113 (assembly-search).
- **Doc-only:** #111 (SV), #94 lowercase/non-ATCG FASTA, #93 CFD NGG-calibration note (in `docs/INPUT_FORMATS.md`).
- **Close at/after release:** #106, #110, #115, #94.
- **External/verify elsewhere:** #63 (bioconda recipe).
