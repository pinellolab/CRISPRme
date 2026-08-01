# CRISPRme Post-Analysis Scalability Analysis

Cost map and ranked deep-fix design for the post-analysis blow-up tracked in
issues **#106** (disk explosion), **#107** (missing guard / hard stop), and
**#108** (OOM in the analysis stage).

## Summary

Post-analysis cost is dominated by the number of **in-budget alignments
retained per site**. That count grows combinatorially with allowed mismatches
and DNA/RNA bulges, and again with local variant density (each variant/IUPAC
position inside a target window can fan a single site out into a powerset of
alternative targets). Empirically, at **mm6 + 2 DNA + 2 RNA bulges** a single
large chromosome produced **~207 GB** of `best*` intermediate tables and OOMed
the analysis at a **~15-58 GB** RAM working set.

This document is the precise cost map (with `file:line` citations verified
against the tree at branch `scale-mitigation`, commit off `origin/main`
`d2e5ab2`) plus a ranked deep-fix design. The safe, collision-free mitigation
shipped alongside it is `PostProcess/search_budget.py` — an estimate/guard that
warns *before* a run blows up (advisory, with an opt-in `--fail-over` hard stop
implementing #107).

> **Line-number note.** The task brief cited approximate lines against an
> earlier revision. The citations below were re-verified against the current
> `PostProcess/new_simple_analysis.py` (903 lines) and
> `PostProcess/remove_contiguous_samples.py` (695 lines). Corrections vs. the
> brief are called out inline as **[corrected: ...]**.

---

## 1. Cost map — RAM

File: `PostProcess/new_simple_analysis.py`

### 1.1 Whole-chromosome `genomeStr` load — line 686
```python
686  genomeStr = inFasta.readlines()  # lettura fasta del chr
687  genomeStr = "".join(genomeStr).upper()
688  ...
689  genomeStr = genomeStr.replace("\n", "")
```
The entire chromosome FASTA is read into a single Python string, uppercased,
and `\n`-stripped — transiently holding **~3x** the chromosome size in memory
(the `readlines` list, the joined string, and the replaced copy all coexist).
For chr1 (~250 Mbp) this is a multi-hundred-MB fixed baseline before any
targets are processed. It is a **fixed** cost (does not scale with the search
budget) but sets the RAM floor.
**[matches brief: ~line 686]**

### 1.2 `iupac_decomposition()` SNP-combination powerset — lines 198-390
```python
198  def iupac_decomposition(split, guide_no_bulge, guide_no_pam, cluster_to_save):
...
258      for size in range(countIUPAC):  # the time of the universe
...
390                          cluster_to_save.append(final_line)
```
For each target overlapping `countIUPAC` variant/IUPAC positions, this builds
`totalDict[count][size]` layer by layer, intersecting sample sets to enumerate
**every co-occurring multi-SNP combination** (the author's own comment on
line 258 reads `# the time of the universe`). This is the powerset fan-out:
targets in variant-dense regions expand into many alternative rows, each
appended to `cluster_to_save` (lines 386-390). This is the primary driver of
both the RAM working set and the downstream `best*` disk volume, and it
compounds multiplicatively with the mismatch/bulge budget.
**[matches brief: ~lines 198-390]**

### 1.3 `cluster_to_save` flushes at a fixed 100k-row count — line 808
```python
808  if len(cluster_to_save) >= 100000:
809      # after reading 100k lines from file and creating the cluster, start processing it
810      clusters_with_scores = calculate_scores(cluster_to_save)
```
The in-memory batch is flushed on a fixed **row count** (`>= 100000`), **not on
cluster boundaries**. In a dense region a single logical cluster can exceed
100k expanded rows, so the flush splits clusters — and, worse, between flushes
the list plus its per-target score copies (see 1.4) are all resident. The
100k threshold is a magic constant; it bounds neither cluster size nor peak
memory in variant-dense regions.
**[corrected: brief said ~line 808 "cluster_to_save flushes at a 100k-row
count" — the flush condition is at line 808, `>= 100000`. Confirmed exact.]**

### 1.4 List copies in `calculate_scores()` — lines 617-629
```python
617  def calculate_scores(cluster_to_save):
...
624      for target in cluster_to_save:  # calculate CFD score for each target
625          target_CFD = target.copy()
626          cluster_with_CFD_score.append(preprocess_CFD_score(target_CFD))
...
629      cluster_with_CRISTA_score = preprocess_CRISTA_score(cluster_to_save)
...
679      return [cluster_with_CFD_score, cluster_with_CRISTA_score]
```
`calculate_scores` builds a **full CFD copy** of the batch (`target.copy()` per
row, line 625, accumulated into `cluster_with_CFD_score`) **and** a full CRISTA
list (line 629), then returns both. So at the flush point the peak is roughly
**3x the batch**: the original `cluster_to_save` plus two scored copies, all
live simultaneously. Combined with 1.2 and 1.3, this is the mechanism behind
the #108 OOM.
**[corrected: brief said ~line 617 "list copies in calculate_scores()" — the
copy is at line 625; the function spans 617-679. Confirmed.]**

---

## 2. Cost map — Disk

File: `PostProcess/new_simple_analysis.py`

### 2.1 `best*` fan-out written per-row, before consolidation — lines 818/820/825 and 870/872/877
Two write blocks emit one line per (expanded) target into three parallel
tables:
```python
# flush block (>= 100k rows), lines 812-826:
818      cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
820      mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
825      crista_best.write("\t".join(target) + "\t" + str(0) + "\n")
# tail block (final < 1M rows), lines 863-877:
870      cfd_best.write("\t".join(target) + "\t" + str(0) + "\n")
872      mmblg_best.write("\t".join(target) + "\t" + str(0) + "\n")
877      crista_best.write("\t".join(target) + "\t" + str(0) + "\n")
```
Every expanded alternative target from `iupac_decomposition` (1.2) is written
to **three** `best*` files (`.bestCFD.txt`, `.bestmmblg.txt`, `.bestCRISTA.txt`)
**before any cluster consolidation happens**. These files hold the full,
un-deduplicated fan-out — on the order of **~48 GB each** at mm6+2+2 on a large
chromosome, i.e. the observed **~207 GB** aggregate (three tables + the
CFDGraph and transient sort copies).
**[corrected: brief cited ~lines 818/870. Confirmed: first block writes at
818/820/825, second block at 870/872/877. The `>= 100000` guard is line 808;
there is no code at "line 870" that differs in intent from the flush block.]**

### 2.2 Consolidation happens later, in `remove_contiguous_samples.py`
File: `PostProcess/remove_contiguous_samples.py`
```python
7    def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
...
37       final_list_best_ref = list()
44       final_list_best_var = list()
80       if sort_order == "score":  # sort each cluster, pick best target(s)
```
`get_best_targets` (line 7) is where each cluster is finally reduced to its
best target(s) — but it runs **downstream**, consuming the already-materialized
`best*` tables. Before it runs, the driver script sorts each full table on
disk. From `PostProcess/submit_job_automated_new_multiple_vcfs.sh`:
```
608  tail -n +2 $final_res.bestCFD.txt    | LC_ALL=C sort ... >> $final_res.tmp && mv $final_res.tmp $final_res.bestCFD.txt
611  tail -n +2 $final_res.bestCRISTA.txt | LC_ALL=C sort ... >> $final_res.tmp && mv $final_res.tmp $final_res.bestCRISTA.txt
614  tail -n +2 $final_res.bestmmblg.txt  | LC_ALL=C sort ... >> $final_res.tmp && mv $final_res.tmp $final_res.bestmmblg.txt
```
Each `sort` writes a **full second copy** (`$final_res.tmp`) of a ~48 GB table
before `mv` replaces the original — the peak-disk doubling that D1 targets.

---

## 3. Ranked deep-fix design

Ranked by value-to-risk. **D1** and **D2** preserve results by design; **D3**
changes results by design and is opt-in.

### D1 — Stream per-chrom `best*` through `sort -m` (drop a ~48 GB copy)

**Change.** Instead of writing each full `best*` table, then re-reading and
re-sorting it with `sort` + `mv` (a full temp copy — §2.2), write each
chromosome's `best*` output **already sorted per chromosome** and use
`sort -m` (merge of pre-sorted streams) at the aggregation step. `sort -m` does
not re-sort; it merges pre-sorted inputs in a single streaming pass, so the
transient `$final_res.tmp` full copy is eliminated. Per-chromosome outputs are
already produced separately, so each is small enough to sort within its own
partition.

**Files touched.** `PostProcess/submit_job_automated_new_multiple_vcfs.sh`
(lines ~608-614) and the per-chrom emit ordering in
`PostProcess/new_simple_analysis.py` (the write blocks, §2.1) — and/or the
merge glue in `PostProcess/merge_close_targets_cfd.sh`.

**Correctness risk.** Low-to-moderate. Requires that per-chrom inputs are
sorted on **exactly** the same key `sort -m` merges on (`-k16,16 -k5,5 -k7,7n
-k21,21rg -k11,11n` for CFD; note the different key for `bestmmblg`). A key
mismatch silently produces a mis-ordered merge, so the sort keys must be
asserted identical at emit and merge time. `LC_ALL=C` must be set on both
sides. No rows are added or dropped.

**Validation recipe.** On a mid-size run (single medium chromosome, mm4+1+1,
a real VCF): capture `integrated_results.tsv` before the change; apply D1;
re-run; `sort` both final TSVs identically and `diff` — expect **byte-identical**
results. Also assert peak disk (`du -s` sampled during the run, or `/usr/bin/time
-v` max RSS-agnostic disk high-water via `df`) drops by roughly one `best*`
table.

### D2 — Flush at cluster boundaries + remove `calculate_scores` copies (fix the #108 OOM)

**Change.** Two coupled edits in `PostProcess/new_simple_analysis.py`:
1. Replace the fixed `len(cluster_to_save) >= 100000` flush (line 808) with a
   **cluster-boundary flush**: accumulate a cluster's expanded rows, and flush
   when the next input row crosses the cluster key
   (`Real_Guide, Chromosome, Cluster_Position`). This bounds the resident batch
   by one cluster instead of an arbitrary 100k rows, and stops splitting
   clusters across flushes (which currently also corrupts the per-cluster
   "best" selection at the boundary).
2. In `calculate_scores` (lines 617-679), drop the per-target `target.copy()`
   (line 625) and the parallel full CFD/CRISTA lists; score **in place** /
   stream each scored row straight to the writer, so peak is ~1x the batch
   instead of ~3x.

**Correctness risk.** Moderate. The boundary detection must match the exact key
the downstream consolidation (`remove_contiguous_samples.get_best_targets`)
assumes, or "best" selection changes. Removing `.copy()` (line 625) requires
confirming `preprocess_CFD_score` / `preprocess_CRISTA_score` do not mutate the
row in a way the *other* scorer then reads — currently the copy hides exactly
that coupling, so it must be untangled (score CRISTA from the original before
CFD mutates, or have each scorer take/return its own row). This is the fix that
directly resolves the #108 OOM.

**Validation recipe.** Same mid-size run as D1. Diff final
`integrated_results.tsv` before/after — expect identical (D2 is a memory-layout
change, not a results change). Additionally, capture peak RAM with
`/usr/bin/time -v` (Max RSS) before/after on a variant-dense chromosome and
confirm the working set no longer scales with region density — it should track
one cluster, not the whole 100k batch.

### D3 — Optional total-edit cap (changes results by design; opt-in)

**Change.** Add an optional `--max-total-edits N` (or reuse an existing total
budget) that hard-caps `mismatches + bulges_dna + bulges_rna` retained per
alignment during `iupac_decomposition` (§1.2), pruning the deepest, lowest-value
tail of the powerset before it is ever written to `best*`.

**Correctness risk.** High **by design** — this **removes** high-edit-distance
targets from the output, so results differ. It must be **off by default**,
clearly documented, and surfaced in the run metadata so a user knows the search
was truncated. It is the last-resort lever when D1+D2 are not enough and the
parameters cannot be reduced.

**Validation recipe.** Not a byte-diff (results intentionally differ). Instead:
run with and without the cap on a mid-size run; confirm every row present in the
**capped** output is present and identical in the **uncapped** output (the cap
only ever removes rows, never alters retained ones); confirm all removed rows
have `Total > N`; quantify the disk/RAM reduction.

---

## 4. Coordination

**D1, D2, and D3 all edit files that are actively being changed by open PRs
#96, #112, and #113** (the PostProcess pipeline rework):
`PostProcess/new_simple_analysis.py`,
`PostProcess/remove_contiguous_samples.py`,
`PostProcess/submit_job_automated_new_multiple_vcfs.sh`, and
`PostProcess/merge_close_targets_cfd.sh`. Landing D1/D2/D3 from a separate
branch would collide with that in-flight work and risk silent result changes in
the exact code paths those PRs are refactoring.

**Therefore D1/D2/D3 must be authored by, or closely coordinated with, the
maintainers** (cc @ManuelTgn @anncir1) on top of / after PRs #96/#112/#113
land, using the validation recipes above as the acceptance gate.

**The only collision-free deliverables** — and the ones shipped in this PR —
are:
- `PostProcess/search_budget.py` — the pre-run estimate/guard (advisory, with
  opt-in `--fail-over` hard stop for #107), and
- this document (`docs/SCALABILITY_ANALYSIS.md`).

Neither touches any file the open PRs are editing.
