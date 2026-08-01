# brute_force_gen (Rust)

A fast, dependency-free (std-only) Rust port of
[`generate_brute_force.py`](../generate_brute_force.py), the exhaustive
brute-force CRISPR off-target ground-truth generator used by CRISPRme's
`validate-test`.

It is a **drop-in replacement** for `generate_brute_force.py`: same CLI, same
semantics, and it produces **byte-for-byte identical output** (as a set of TSV
rows). It is validated against the committed reference
(`../brute-force-1000G/brute_force_1000G.tsv`, chr22 rows: 3495 rows, 0 missing,
0 extra — case-insensitive on the RNA/DNA columns).

## Build

```sh
cargo build --release
```

This produces `./target/release/brute_force_gen`. No external crates are
required (Rust std only, edition 2021).

## Usage

Identical flags to the Python generator:

```
--fasta        Per-chromosome (enriched) single-record FASTA
--rna          Guide spacer+PAM (IUPAC allowed)
--max-mismatches   (default 4)
--max-dna-gaps     (default 1)
--max-rna-gaps     (default 1)
--chrom        Chromosome name written to the CHR column (e.g. chr22)
--output       Output TSV path
--pam          Concrete PAM bases with wildcard N removed
               (e.g. GG for NGG, TTTV for Cas12a); empty disables PAM filtering
--pam-5prime   PAM is at the 5' end of the guide (Cas12a); default 3' (Cas9)
```

### Cas9 example (matches the committed reference)

```sh
./target/release/brute_force_gen \
    --fasta chr22.enriched.fa \
    --rna CTAACAGTTGCTTTTATCACNGG \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --chrom chr22 --pam GG \
    --output cas9_chr22.tsv
```

### Cas12a example (5' PAM)

```sh
./target/release/brute_force_gen \
    --fasta chr22.enriched.fa \
    --rna TTTVCCTTGTCAAGGCTATTGGTC \
    --max-mismatches 4 --max-dna-gaps 1 --max-rna-gaps 1 \
    --chrom chr22 --pam TTTV --pam-5prime \
    --output cas12a_chr22.tsv
```

## Output

Tab-separated with the same header and columns as the Python generator:

```
CHR  RNA  DNA  Strand  Start  END  mismatches  gaps_in_RNA  gaps_in_DNA
```

## Notes

- The aligner enumerates **all** in-budget alignments (mismatch / DNA-gap /
  RNA-gap budget), not just the optimal-scoring one, then applies the PAM +
  edge-gap filters — exactly matching `generate_brute_force.py` and the
  vendored reference algorithm it is derived from.
- The FASTA reader expects a single-record (per-chromosome) FASTA; it skips the
  `>` header and concatenates the sequence lines, preserving case (variant
  enrichment encodes SNVs as IUPAC bases).
