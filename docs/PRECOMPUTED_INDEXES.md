# Precomputed CRISPRme indexes on HuggingFace

For a given `(genome, PAM, VCF-set)`, CRISPRme's two slowest one-time stages are:

1. **Genome enrichment** — `crispritz.py add-variants` builds the variant-enriched
   FASTAs, the per-chromosome SNP dictionaries (`my_dict_<chr>.json`) and the
   fake-indel FASTAs.
2. **TST index build** — `crispritz.py index-genome` builds the
   `genome_library/<true_pam>_<bMax>_<name>` trie indexes.

These artifacts are fully reusable — once built for the default reference they
never change. `PostProcess/crisprme_index.py` packages them and ships them
through the CRISPRme HuggingFace dataset so end users can download a precomputed
index and skip both stages.

The tool reuses the HF plumbing from `PostProcess/utils.py`: `HF_DATA_REPO`
(env `CRISPRME_HF_DATA_REPO`, default `lucapinello/crisprme-data`) and the same
`repo_type="dataset"` convention as `hf_fetch`.

## HuggingFace layout

Every index lives under the `indexes/<name>/` prefix of the HF dataset:

```
indexes/<name>/
    manifest.json                                       # provenance
    genome_library/<true_pam>_<bMax>_<name>/...         # reference TST index
    genome_library/<true_pam>_<bMax>_<name>+<vcf>/...   # variant TST index (if --vcf)
    Genomes/<name>+<vcf>/...                             # enriched genome (if --vcf)
    Genomes/<name>+<vcf>_INDELS/...                      # fake-indel FASTAs (if --vcf)
    Dictionaries/...                                     # variant dictionaries (if --vcf)
```

`manifest.json` records provenance so a downloaded index is self-describing:

```json
{
  "name": "hg38",
  "genome": "hg38",
  "pam": "20bp-NGG-SpCas9.txt",
  "pam_name": "NGG",
  "bMax": 2,
  "vcf": "1000G",
  "thread": 8,
  "crispritz_version": "...",
  "created_by": "...",
  "created_at": "2026-...T...Z",
  "note": "default reference index for CRISPRme+"
}
```

## Build an index (maintainers, on the big machine)

crispritz must be installed (it ships in the CRISPRme conda/Docker environment).
The reference index is built for **bMax, bMax+1 and 1** — matching the pipeline's
`idx_folder1/2/3` scheme — so any search up to the requested bulge count can reuse it.

```bash
# reference-only index
python PostProcess/crisprme_index.py build \
    --genome Genomes/hg38 \
    --pam PAMs/20bp-NGG-SpCas9.txt \
    --name hg38 \
    --bmax 2 --thread 8 \
    --work-dir /scratch

# reference + variant (enriched) index
python PostProcess/crisprme_index.py build \
    --genome Genomes/hg38 \
    --pam PAMs/20bp-NGG-SpCas9.txt \
    --name hg38 \
    --vcf VCFs/1000G \
    --bmax 2 --thread 8 \
    --work-dir /scratch \
    --note "default 1000G reference index for CRISPRme+"
```

Artifacts land in `<work-dir>/indexes-build/<name>/` (`genome_library/`,
`Genomes/`, `Dictionaries/`, `manifest.json`).

## Upload an index (maintainers)

Requires an HF token with write access to the dataset repo (`hf auth login`).

```bash
python PostProcess/crisprme_index.py upload --name hg38
# custom source folder / repo:
python PostProcess/crisprme_index.py upload \
    --name hg38 --path /scratch/indexes-build/hg38 --repo lucapinello/crisprme-data
```

Uploads via `huggingface_hub.upload_large_folder` (multi-worker, resumable) to
`indexes/hg38/`, falling back to `upload_folder`, then prints the tree URL.

## List available indexes (anyone)

```bash
python PostProcess/crisprme_index.py list
# ->
# Precomputed indexes on lucapinello/crisprme-data (indexes/<name>/):
#   hg38    (12.4 GB)
```

## Download an index and skip enrichment/indexing (end users)

```bash
cd /path/to/crisprme-working-dir
python PostProcess/crisprme_index.py download --name hg38 --dest .
```

The tool snapshot-downloads `indexes/hg38/*` and **promotes** `genome_library/`,
`Genomes/` and `Dictionaries/` to the root of `--dest`. Run your CRISPRme search
from that directory (as the current working directory) and it will find the
precomputed index and enriched genome instead of rebuilding them, saving the two
slowest stages. The manifest is copied to `manifest.<name>.json` for reference.
