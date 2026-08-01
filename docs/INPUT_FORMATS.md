# CRISPRme Input Formats and Behavior Guide

This document explains how to prepare CRISPRme inputs (PAM, guide, VCF,
annotation, and genome files) and how a few key parameters change the results.
It is example-driven and written to help you avoid the most common mistakes that
silently produce wrong or empty output.

Every rule below was verified against the code in this repository. Where a
behavior is enforced by a specific line, the file and line are cited so you can
check it yourself.

## Table of contents

1. [PAM file format](#1-pam-file-format)
2. [Guide file format](#2-guide-file-format)
3. [PAM filename convention](#3-pam-filename-convention)
4. [One PAM motif per file](#4-one-pam-motif-per-file)
5. [How the PAM affects results and cluster merging](#5-how-the-pam-affects-results-and-cluster-merging)
6. [Chromosome naming](#6-chromosome-naming)
7. [Annotation BEDs](#7-annotation-beds)
8. [VCF requirements](#8-vcf-requirements)
9. [Pre-flight validation](#9-pre-flight-validation)
10. [Bulge reporting](#10-bulge-reporting)
11. [Guides with no hits](#11-guides-with-no-hits)
12. [Large / expensive searches](#12-large--expensive-searches)

---

## 1. PAM file format

A PAM file is a **single line** with two whitespace-separated fields:

```
<motif> <length>
```

- `<motif>` is the full spacer+PAM window written as a string, with the spacer
  region as `N`s.
- `<length>` is the PAM length **and its position**: a **positive** value means
  the PAM is on the **3′** side of the protospacer (like SpCas9); a **negative**
  value means the PAM is on the **5′** side (like Cas12a/Cpf1).

The sign of the length is what tells CRISPRme where the PAM sits, and its
absolute value is the PAM length. This is parsed in `crisprme.py` (~lines
1198–1211): if the length is `< 0`, the PAM is taken from the **start** of the
motif (`pam_begin = True`); otherwise it is taken from the **end**.

### 3′ PAM — SpCas9 (NGG)

The committed example `test/data/PAMs/20bp-NGG-SpCas9.txt` contains exactly:

```
NNNNNNNNNNNNNNNNNNNNNGG 3
```

That is 20 `N`s for the 20 nt spacer, followed by `NGG` (3 nt PAM), with length
`3` (positive → PAM is 3′). The total motif is 23 characters.

### 5′ PAM — Cas12a (TTTV)

The committed example `test/data/PAMs/TTTV-20bp-Cas12a.txt` contains exactly:

```
TTTVNNNNNNNNNNNNNNNNNNNN -4
```

That is `TTTV` (4 nt PAM) followed by 20 `N`s for the 20 nt spacer, with length
`-4` (negative → PAM is 5′). The total motif is 24 characters.

### Cas12a: WRONG vs RIGHT

This is the single most common Cas12a mistake.

**WRONG** — positive length and the literal PAM bases inside the guide:

```
# PAM file (WRONG for Cas12a)
TTTVNNNNNNNNNNNNNNNNNNNN 4
```

A positive length tells CRISPRme the PAM is 3′, so it will slice the PAM off the
**wrong end** of the motif. Likewise, putting the literal `TTT` at the front of
your **guide** (see §2) is wrong.

**RIGHT** — negative length, and the PAM region represented as part of the
motif:

```
# PAM file (RIGHT for Cas12a)
TTTVNNNNNNNNNNNNNNNNNNNN -4
```

The `-4` places the 4 nt PAM on the 5′ side; the guide file then uses `N`s for
the PAM region (see §2).

---

## 2. Guide file format

The guide file lists one spacer per line. Two rules matter:

1. **The guide length must equal the spacer length of the PAM motif.** For a
   Cas12a `TTTV` + 20 nt spacer motif, each guide is 24 characters
   (4 for the PAM region + 20 for the spacer).
2. **The PAM region of the guide is written as `N`s** — you do not put the
   literal PAM bases (e.g. `TTTV` or `TTT`) into the guide.

The committed example `test/data/Guides/hbg_cas12a_test_guide.txt` contains
exactly:

```
NNNNCCTTGTCAAGGCTATTGGTC
```

Here the leading `NNNN` is the 5′ PAM region (matching `TTTV`, 4 nt) and
`CCTTGTCAAGGCTATTGGTC` is the 20 nt spacer. The guide is 24 characters, matching
the 24-character motif from the Cas12a PAM file in §1.

For a 3′ SpCas9 guide, the `N`s go at the **end** (e.g. a 20 nt spacer followed
by `NNN` for the `NGG` PAM region), matching the 23-character SpCas9 motif.

---

## 3. PAM filename convention

CRISPRme parses a **nuclease label** out of the PAM *filename* with:

```python
nuclease = os.path.basename(pamfile).split(".")[0].split("-")[2]
```

(`crisprme.py`, line 1215.) It strips the extension, splits the basename on `-`,
and takes the **third** field (index 2). So the filename needs at least three
`-`-separated fields:

```
<length>-<motif>-<Cas>.txt
```

for example `20bp-NGG-SpCas9.txt` or `TTTV-20bp-Cas12a.txt`.

**Important:** the three fields are cosmetic labels used for naming/reporting.
The **file contents** (the motif and the signed length on line 1) are what
actually drive the search. If your filename has fewer than three `-`-separated
fields, the label parse will fail (`IndexError`); if the label is "wrong" but
the contents are right, the search still runs correctly — only the reported
nuclease label is affected.

---

## 4. One PAM motif per file

CRISPRme reads **only the first line** of the PAM file:

```python
pam_char = pam_file.readline()
```

(`crisprme.py`, line 1199.) A PAM file with multiple lines will **silently** use
only line 1 — the rest is ignored, with no warning.

To search **more than one PAM motif at once**, do not stack multiple lines.
Instead, encode the alternatives as a **single IUPAC motif**. For example, to
match both `NAG` and `NGG` PAMs, use `NRG` (`R` = A or G):

```
NNNNNNNNNNNNNNNNNNNNNRG 3
```

CRISPRme expands IUPAC codes during the search, so one degenerate motif covers
the set.

---

## 5. How the PAM affects results and cluster merging

### The PAM is part of the score

The PAM is not just a filter — it contributes to the **CFD score**. In
`PostProcess/new_simple_analysis.py` the score is multiplied by a
PAM-dependent factor:

```python
score *= pam_scores[pam]
```

(`PostProcess/new_simple_analysis.py`, line 123, inside `calc_cfd`.) So an
otherwise identical off-target can score differently depending on the PAM it
was found with (e.g. an `NGG` protospacer scores higher than the same site with
an `NAG` PAM).

### Cluster merging: two output files

Off-target hits that fall at the **same genomic site** (contiguous
positions/alignments) are grouped into a cluster, and CRISPRme keeps only the
**single best-scoring representative** of each cluster in the main results. See
`retrieve_best_target` in `PostProcess/merge_contiguous_targets.py` (~lines
440–528): the best target is written to the primary output, and every remaining
member of the cluster is written to the alternative-alignments output.

This produces two result tables:

- **`*_integrated_results.tsv`** — one row per cluster: the best-scoring
  representative, annotated and used for downstream figures/summaries.
- **`*_all_results_with_alternative_alignments.tsv`** — the alternative
  alignments that were merged away. For example, if the same locus is reachable
  with both an `NGG` and an `NAG` PAM and the `NGG` variant wins, the `NGG` row
  goes to `*_integrated_results.tsv` and the `NAG` row goes to
  `*_all_results_with_alternative_alignments.tsv`.

If you are hunting for a specific alignment (a particular PAM variant, bulge
configuration, or SNP combination) and don't see it in the integrated results,
check the alternative-alignments file — it was very likely merged into a
higher-scoring representative rather than dropped.

---

## 6. Chromosome naming

CRISPRme expects **UCSC-style, `chr`-prefixed** chromosome names throughout your
FASTA, VCF, and BED inputs:

- Supported: `chr1` … `chr22`, `chrX`, `chrY`, `chrM`.
- **Not** supported: Ensembl-style bare IDs (`1`, `2`, `X`), and
  punctuated/unplaced contigs (e.g. `KI270728.1`).

The genome indices are `chr`-prefixed, and the pipeline normalizes/expects the
prefix (see `PostProcess/pool_post_analisi_indel.py`, `_normalize_chrom`, lines
21–26: *"the genome indices use 'chr'-prefixed names"*). Mixing naming schemes
between your FASTA, VCF, and BED files will cause chromosomes to silently not
match.

**VCF filenames must be dot-segmented.** CRISPRme identifies the chromosome of
each VCF by splitting the *filename* on `.` and finding a segment that starts
with `chr`. So the chromosome must be its own dot-separated segment, e.g.
`MyCohort.chr1.vcf.gz` or `ALL.chr1.filtered.vcf.gz`. A name like
`MyCohort_chr1.vcf.gz` (chromosome buried in an underscore-delimited prefix)
will **not** be recognized. See `docs/crisprme_data_setup_051826.md`
("File naming requirement", ~line 234).

**Trailing blank lines are fine.** Config and sample `.txt` files may have
trailing blank lines; empty lines are ignored during parsing.

---

## 7. Annotation BEDs

CRISPRme accepts two distinct annotation inputs, and they route to **different
result columns** — they are not interchangeable:

- **`--annotation`** — a functional/regulatory BED (e.g. ENCODE cCREs). Its
  fourth column is used as the annotation label.
- **`--gene_annotation`** — a gene-model BED (e.g. GENCODE) used to locate the
  closest gene / genic context.

Both accept either a plain `.bed` or a bgzipped `.bed.gz` file (the sort helper
transparently handles both; see `_sort_annotation` in `crisprme.py`, lines
694–713).

**Sorting.** BEDs must be sorted with BEDOPS/GRCh ordering. As of 2.1.12,
**both `--annotation` and `--gene_annotation` are auto-sorted for you** — each is
routed through `_sort_annotation` (via `_check_annotation` / `_check_gene_annotation`
in `crisprme.py`), which decompresses if needed, sorts with `sort-bed`, and
re-compresses. You can pass either a plain `.bed` or a `.bed.gz`, sorted or not, and
CRISPRme normalizes it.

---

## 8. VCF requirements

For a variant-aware search, each VCF record needs:

- a **`FILTER`** value — use `PASS` (records not marked `PASS` may be dropped by
  filtering);
- **`INFO/AF`** — a **numeric** allele frequency;
- **`FORMAT/GT`** — genotype fields for the samples.

### Variant types

- **SNPs and indels** are both supported. A VCF containing **only SNPs** or
  **only indels** works fine — the missing category simply yields an empty
  branch of results.
- **`complex` and `mnp` records are treated as indels.**
- **Structural variants are not supported.** Symbolic ALT alleles such as
  `<DEL>`, `<INS>`, `<DUP>` are skipped, and breakend (BND) notation can corrupt
  parsing — **remove SV records before running.**

### Recommended prep

Normalize and pre-filter your VCF before a run, for example:

```bash
# keep PASS records, drop symbolic/SV alleles, ensure chr-prefixed names
bcftools view -f PASS input.vcf.gz \
  | bcftools view -e 'ALT[0]~"<"' \
  -Oz -o clean.vcf.gz
```

and confirm every record carries a numeric `INFO/AF` and a `FORMAT/GT` field.
See §9 for catching these problems before a long run.

---

## 9. Pre-flight validation

The most expensive failures are the ones you discover **after** a multi-hour
search — a missing `INFO/AF`, a non-`PASS` `FILTER`, an absent `FORMAT/GT`, a
sample-ID mismatch, or a chromosome-naming mismatch.

> **Note (current `main`):** a dedicated input validator
> (`PostProcess/validate_inputs.py`, and a `--full_input_validate` option) is
> **not yet present** in this repository — it is referenced as a forthcoming,
> PR-added module in `.github/workflows/validate-benchmarks.yml` (line 47) but
> has not been merged. Until it lands, validate inputs manually before a long
> run:
>
> - `bcftools view -h your.vcf.gz` — confirm `FILTER`, `INFO/AF`, and
>   `FORMAT/GT` are declared in the header.
> - `bcftools query -f '%CHROM\n' your.vcf.gz | sort -u` — confirm `chr`-prefixed
>   chromosome names (see §6).
> - Cross-check that your sample IDs in the VCF match your samplesIDs file.
>
> (Do **not** confuse this with `PostProcess/validate.py` / the `validate-test`
> subcommand, which checks *off-target predictions* against brute-force
> benchmarks — it is a correctness regression test, not an input validator.)

---

## 10. Bulge reporting

CRISPRme reports bulges with a **single combined `Bulge_Size` column**
(DNA + RNA together) plus a **`Bulge_type` column** that indicates which kind of
bulge it is. There are currently **no separate per-type (DNA vs RNA) bulge
counts** in the output. See the result header in
`PostProcess/analisi_indels_NNN.py` (line 404), which lists `Bulge_type` and a
single `Bulge_Size` field.

---

## 11. Guides with no hits

A guide with **zero in-budget off-targets** (given your mismatch/bulge limits)
produces **empty results**, not an error. An empty results table for a guide is
a valid outcome, not a failure — it simply means nothing matched within the
allowed budget.

---

## 12. Large / expensive searches

Search cost grows **combinatorially** with mismatches and DNA/RNA bulges, and
again with local variant density. Allowing **more than one bulge in both DNA and
RNA simultaneously** is especially expensive: a single chromosome at mm6 + 2 DNA
+ 2 RNA bulges has produced **~207 GB of intermediate files and OOMed** the
analysis (see the header of `PostProcess/search_budget.py`).

**Estimate before you launch.** A dependency-light heuristic ships in the repo:

```bash
python PostProcess/search_budget.py --mm 6 --bdna 2 --brna 2
```

It returns a coarse resource estimate and a **risk tier** (`OK`, `HEAVY`,
`EXTREME`). An `EXTREME` tier corresponds to roughly **~200 GB of disk / OOM
RAM** — treat it as "probably don't run this as-is." You can also make the
estimate a hard gate with `--fail-over EXTREME` (exits non-zero at/above that
tier). Every number is a heuristic, not a guarantee — use it as a
go / be-careful / don't signal.

---

## See also

- `docs/crisprme_data_setup_051826.md` — setting up genomes, VCF datasets, and
  annotations.
- `README.md` — installation and full command reference.
