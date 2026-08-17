# Combined multi-source VCF panels for CRISPRme (single-scan enrichment)

**Status:** validated **genome-wide** (1000G + HGDP) — the full AF-provenance panel
was built and both the **NGG** and **pamless NNN** combined indexes were produced
from a single enrichment, with pooled AF confirmed correct. Opening this to compare
with Ann's merge/provenance approach and converge on one.

## Motivation
CRISPRme enriches a reference with a VCF's variants and scans once. Running N
resources (1000G, HGDP, AllOfUs, TOPMed, gnomAD, …) as N separate datasets means N
enrichments + N scans + a post-hoc result merge. **Merging the resources into one
panel first → enrich once, scan once** — a large, growing time saving. The catch is
doing the merge without losing information we care about (per-source allele
frequency and which dataset a variant came from).

## What the merge produces (provenance, per variant)
| INFO field | Meaning |
|---|---|
| `AF` | **Pooled** allele frequency, recomputed as AC/AN across ALL merged samples (`bcftools +fill-tags`). This is what CRISPRme's enricher reports. |
| `AF_<SRC>` | Each source's **original** AF, verbatim; `.` where that source lacked the variant. |
| `SRC` | Comma-list of sources containing the variant: `1000G,HGDP` (shared), `HGDP` (HGDP-unique), … Derived from which `AF_<SRC>` are present. |

So provenance is recoverable two ways: implicitly (`AF_1000G` present ⇔ in 1000G) and
explicitly (`SRC`). A variant **unique** to one source is immediately identifiable.

## How it works (see `merge_vcf_panels.sh`)
Per chromosome, for each source: normalize contigs to `chr`-prefixed (to match the
hg38 reference), rename `INFO/AF → INFO/AF_<SRC>`, index. Then
`bcftools merge -m none` the sources, `bcftools +fill-tags -- -t AF` to compute the
pooled `AF`, and derive `INFO/SRC` from a fast sites-only pass over the `AF_<SRC>`
fields. `build_combined_panel.sh` runs this over all chromosomes (parallel,
resumable), assembles the panel `samplesID` (union), and calls
`crisprme.py build-index-only` to enrich + index.

## Validation

### Genome-wide (1000G + HGDP) — the full build
- **All 23 chromosomes merged** (chr1–22 + chrX), 3,477 samples each, AF-provenance
  intact (`AF` pooled + `AF_1000G`/`AF_HGDP` + `SRC`).
- **Both combined indexes built from ONE shared enrichment:**
  `NGG_3_hg38+hg38_1000G_HGDP` (8.7 GB + 730 MB INDELS) and the pamless
  `NNN_3_hg38+hg38_1000G_HGDP` (144 GB + 15 GB INDELS).
- **Pooled AF correct at scale:** `chr22:10516173 → AF=0.024992` (`AF_1000G=0.02`,
  `AF_HGDP=0.030541`, `SRC=1000G,HGDP`).
- **Why this matters — the earlier plain merge was wrong:** a naive
  `bcftools merge` (no `+fill-tags`) leaves shared 1000G∩HGDP variants carrying only
  ONE source's AF (e.g. the HGDP 0.0305 instead of the pooled 0.0250 above) — 15/1000
  sampled sites mismatched. This tooling recomputes the pooled AF, so the shipped
  combined index reports correct frequencies.
- **max-edits tolerance:** the pamless index is searched with `--max-total-edits 5`
  (CRISPRitz ≥2.8.1); this bounds the pamless search-space explosion in
  variant-dense regions (validated separately).

### chr22 detail (provenance cross-check)
- **Samples:** 3,477 = 2,548 (1000G) + 929 (HGDP); cohorts disjoint (0 shared).
- **Record provenance cross-checks exactly:** SRC = 572,013 `1000G` + 697,942 `HGDP`
  + 487,066 `1000G,HGDP` = 1,757,021 total; and 572,013+487,066 = **1,059,079**
  (= 1000G source records) and 697,942+487,066 = **1,185,008** (= HGDP source records).
- **Pooled AF correct:** e.g. chr22:10516173 `AC=156 AN=6242 → AF=0.024992`, with
  `AF_1000G=0.02`, `AF_HGDP=0.030541` preserved.
- **Enricher reads the pooled AF:** CRISPRitz #36's exact-key `AF` match selects
  `AF=0.024992`, correctly skipping the earlier `AF_1000G`/`AFR_AF` fields (the old
  2-char-prefix match would have grabbed `AF_1000G`).
- **End-to-end:** the panel enriches hg38_chr22 and builds the pamless `NNN_3` index
  (+ INDELS) without error. [`e2e_afmerge.sh`]

## Trade-offs (intentional, documented)
- **Phasing is not preserved across sources.** Cohorts are disjoint, so cross-source
  multi-variant haplotypes have no real carrier; single-variant off-targets (the vast
  majority) are fully correct. The single-scan speed win is the point.
- **Sites-only sources** (gnomAD, TOPMed: no per-sample genotypes) contribute variant
  positions + `AF_<SRC>` but no samples; the pooled `AF` then reflects only the
  genotyped sources. Add them as extra `SOURCES` entries.
- `bcftools merge -m none` can emit same-POS records for multiallelic sites; #36's
  enricher guards the AF-count mismatch this can create.

## Dependencies
- `bcftools >= 1.18` with the `+fill-tags` plugin; `htslib` (`bgzip`/`tabix`).
- **CRISPRitz #36** (enricher AF/FILTER robustness) for enrichment (exact-key pooled
  AF read + multiallelic AF-count guard) and **#37/#38** (`--max-edits`) for the
  bounded pamless search — both shipped in **CRISPRitz ≥ 2.8.1**.
- **Scratch space:** the per-chromosome intermediates are multi-GB and a genome-wide
  N-way-parallel run writes tens–hundreds of GB of temp. `merge_vcf_panels.sh` now
  puts its work dir on `$TMPDIR` (fall back: the output volume) rather than `/tmp` —
  point `TMPDIR` at a volume with ample free space (a genome-wide 1000G+HGDP run
  needs well over 100 GB of temp + ~160 GB for the final pamless index).

## Reproduce
```bash
# one chromosome (configure SOURCES in the script):
BCFTOOLS=bcftools ./merge_vcf_panels.sh /path/to/cwd/VCFs/hg38_1000G_HGDP chr22
./validate_merge.sh /path/to/cwd/VCFs/hg38_1000G_HGDP/merged.chr22.vcf.gz 1000G HGDP
# genome-wide + enrich + index:
./build_combined_panel.sh /path/to/crisprme_working_dir hg38_1000G_HGDP 20bp-NNN-NO-PAM 2 2 4
```

## To compare with Ann
- Merge mechanics: `-m none` vs `-m both`/normalization; how you handle multiallelic
  records and same-POS duplicates.
- Provenance: is `SRC` (+ `AF_<SRC>`) the same shape you use, or a different tag?
- Pooled vs per-source AF: do you recompute pooled AF, keep per-source, or both?
- Sites-only sources (gnomAD/TOPMed) and phasing policy.
