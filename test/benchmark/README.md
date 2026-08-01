# CRISPRme validation benchmark (brute-force ground truth)

`crisprme.py validate-test` checks a CRISPRme installation by comparing the
off-targets that `complete-test` predicts against an **independent brute-force
ground truth** generated with a different method (exhaustive dynamic-programming
alignment instead of CRISPRme's TST index search).

## Files
- `brute-force-1000G/brute_force_1000G.tsv` — Cas9 reference (SpCas9, NGG PAM,
  guide `sg1617`, thresholds mm=4 / DNA-bulge=1 / RNA-bulge=1) on the hg38 +
  1000G variant-enriched genome.
- `brute-force-cas12-1000G/brute_force_cas12_1000G.tsv` — Cas12a reference
  (TODO: generate with the pipeline below).
- `generate_brute_force.py` — the ground-truth generator (see credit below).

## Reproducible generation pipeline
```bash
# 1. Build the variant-enriched genome the same way CRISPRme does internally
#    (per-chromosome FASTA with 1000G SNVs encoded as IUPAC bases).
#    This is the `add-variants` step of a normal complete-search/complete-test.

# 2. Run the brute-force generator per chromosome, e.g. for Cas9 (sg1617):
python test/benchmark/generate_brute_force.py \
    --fasta <chrN.enriched.fa> \
    --rna CTAACAGTTGCTTTTATCACNGG \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --chrom chrN --output chrN.cas9.tsv

#    ...and for Cas12a (TTTV 5' PAM, 23 nt spacer): put the PAM at the 5' end
python test/benchmark/generate_brute_force.py \
    --fasta <chrN.enriched.fa> \
    --rna TTTV<CAS12_SPACER> \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --chrom chrN --output chrN.cas12.tsv

# 3. Concatenate per-chromosome TSVs, sort, and record the MD5 (BFTARGETSMD5 in
#    PostProcess/utils.py) so validate.py can integrity-check the download.
```

Thresholds are kept at **mm=4 / bulge 1+1** to keep both the search and the
brute-force generation cheap (matches `complete_test.py`). Runtime note: the
generator is a pure-Python exhaustive scanner (~minutes per Mb), so a full
chromosome is a batch job.

## Important: genome provenance
The brute force must be generated from **the exact same variant-enriched genome
CRISPRme searches** (its `add-variants` output). Using a differently-enriched
FASTA yields spurious/missing hits in variant-dense regions (an all-IUPAC window
matches almost anything). Pin the enriched genome build alongside the reference.

## Credit / license
The dynamic-programming brute-force checker is by **Benjamin Vyshedskiy**
(https://github.com/benjaminvyshedskiy/Dynamic_checker). `generate_brute_force.py`
vendors it with attribution; only output formatting (tab-delimited, underscore
column names) was adapted. **TODO (maintainers):** confirm licensing/permission
with the author and add an appropriate license header, as the source repository
does not currently include a license.
