# Settings / Data Manager — add data from the browser

CRISPRme's web interface has a **Settings / Data Manager** page where you add
everything a search needs — reference genomes, indexes, variant (VCF) datasets,
annotations and PAMs — without touching the command line. Whatever you add lands
in your local data folder and shows up in the search form automatically.

Open it from the **"Settings / Data Manager"** button on the home page, or go to
`/settings` directly (e.g. `http://localhost:8080/settings`).

> Data management runs only in **local mode** (how you run CRISPRme on your own
> machine / Docker). On the public web server the page is read-only, because it
> would change shared data.

## What you can do

### Add a reference genome
Type an assembly name and pick a source:
- **UCSC** — downloads the assembly by name from UCSC (e.g. the pig genome
  `susScr11`). CRISPRme handles both UCSC layouts (per-chromosome archive or a
  single multi-FASTA) automatically.
- **HuggingFace** — pulls a genome we host on the CRISPRme data repo.
- **Direct URL** — any `.fa.gz` / `.tar.gz` link.
- **Upload** — pick a local genome file; it is sent in chunks, so there is no
  browser size limit.

The same thing from the CLI:
```bash
crisprme.py download --what genome --ref susScr11 --source ucsc
```

### Add an index (or build one)
Bulge-enabled searches need a genome **index**.
- **Download** a ready-made index from HuggingFace by name (e.g. `NGG_2_hg38`).
- **Build locally** from an installed genome + PAM + bulge counts.
- **Pre-build a variant-aware index**: pick a VCF dataset in the build panel and
  CRISPRme enriches the genome with the variants and indexes it ahead of time, so
  your first variant-aware search is fast instead of slow. (CLI equivalent:
  `crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NGG-SpCas9.txt --bDNA 1 --bRNA 1 --vcf VCFs/1000G`.)

### Add a VCF dataset
Variant datasets are large, so they are **fetched server-side** rather than
uploaded from the browser: from HuggingFace, a URL, an existing folder already on
the machine (registered in place), or a chunked upload for a local file.

### Add an annotation
Upload a `.bed` file (annotations are small).

### Define a nuclease / PAM
Fill in a short form — name, guide length, motif (e.g. `NGG`, `TTTV`), and whether
the PAM is at the 3′ (Cas9-like) or 5′ (Cas12a-like) end — and CRISPRme writes the
PAM file for you.

## Storage and cleanup
The right-hand **Storage** panel shows how much disk each category
(genomes / VCFs / indexes / annotations) is using, the total, and how much free
disk you have. Each download card also states the space it typically needs, so
you can check before starting.

The **Remove installed data** card lets you delete any downloaded genome, index,
VCF dataset or annotation to free space — with a confirmation prompt. Anything you
delete can be downloaded again later.

## How long operations behave
Downloads and index builds run in the background on a dedicated slot (so they
never slow down searches) and report **live progress** in the Activity panel.
When one finishes, the installed-data and storage panels update automatically and
the new data is immediately available in the search form.
