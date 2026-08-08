# Precomputed CRISPRme indexes on HuggingFace

Bulge-enabled CRISPRme searches need a CRISPRitz **index** of the reference
genome. Building it is the single most expensive one-time step of a search. That
index depends only on the **genome + PAM + bulge count** (not on the guides or
the variant dataset), so it can be built once and reused — or built once by a
maintainer, published to HuggingFace, and downloaded by everyone else.

This document covers the download → build → publish workflow. It complements
Section 3.5 of the data-setup guide (`docs/crisprme_data_setup_051826.md`).

## Layout on HuggingFace

Indexes live under `indexes/` in the CRISPRme dataset repo (default
`lucapinello/crisprme-data`, override with `--hf-repo` or `CRISPRME_HF_REPO`):

```
indexes/
  NGG_2_hg38.tar.gz          # SpCas9 (NGG), up-to-1-bulge index of hg38
  TTTV_2_hg38.tar.gz         # Cas12a (TTTV), up-to-1-bulge index of hg38
  ...
```

Each tarball unpacks to a single `genome_library/<name>/` directory plus a
`manifest.json` describing its provenance (PAM, bulge level, genome, build
timestamp). The folder name is `<PAM>_<bMax+1>_<ref>` — exactly what
`complete-search` looks for — so a downloaded index is used with no extra steps.

## Download (end users)

```bash
# fetch a prebuilt index straight into genome_library/
crisprme.py download --what index --index-name NGG_2_hg38 --path "$CRISPRME_DIR"

# then search, pointing at that library (or just run from $CRISPRME_DIR)
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NGG-SpCas9.txt \
  --guide my_guide.txt --mm 4 --bDNA 1 --bRNA 1 \
  --index-path "$CRISPRME_DIR/genome_library" \
  --output my_search
```

The search finds the prebuilt index and skips the build entirely.

## Build (maintainers)

`build-index-only` builds the index without running a search. Pass the **same**
`--genome`/`--pam`/`--bDNA`/`--bRNA` end users will search with, so the folder
name matches:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NGG-SpCas9.txt \
  --bDNA 1 --bRNA 1 --thread 16 --path "$CRISPRME_DIR"
# -> genome_library/NGG_2_hg38/
```

## Publish (maintainers)

Upload the built index to the dataset repo (needs an HF **write** token — via
`--token` or the `HF_TOKEN` env var; never commit it):

```bash
export HF_TOKEN=hf_...        # your write token, in the shell only
crisprme.py publish-index --index genome_library/NGG_2_hg38
# -> uploaded to indexes/NGG_2_hg38.tar.gz (with a manifest.json inside)
```

## manifest.json

Every published index carries a small provenance manifest inside its tarball:

```json
{
  "name": "NGG_2_hg38",
  "created_at": "2026-08-05T12:00:00+00:00",
  "pam": "NGG",
  "index_bmax": "2",
  "genome": "hg38"
}
```

It is surfaced (build timestamp) when the index is downloaded, and otherwise
ignored — the index directory itself is what `complete-search` consumes.

## Notes

- An index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`;
  a different PAM or a higher bulge count needs its own index.
- Variant-enriched (genome + VCF) indexes are dataset-specific and are built
  locally per run; only the reference index is shared this way.
- If `--index-path` is given but no matching index is found there,
  `complete-search` fails fast with a clear message rather than silently
  rebuilding — so a missing/wrong download is caught immediately.
