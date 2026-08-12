# CRISPRme: Setup and Usage Guide

Installation (Section 1) and setup of the default data (Section 2) need to be completed only
once. After that, running a new search on this data simply involves editing a single guide file
and re-running one search command (Section 3). Information on how to include a new VCF dataset (Section 4) and a new PAM (Section 5) are also detailed.

---

## Section 1. Install CRISPRme

*Estimated time: ~5 minutes on a standard internet connection.*

### Configure Bioconda channels (one-time)

```bash
mamba config --add channels bioconda
mamba config --add channels defaults
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```

### Create and activate the environment

```bash
mamba create -n crisprme python=3.8 crisprme -y
mamba activate crisprme
```

### Verify the installation

```bash
crisprme.py --version
```

---

## Section 2. Set up the default reference data (one-time)

*Estimated disk usage: ~410 GB total (~16 GB 1000G + ~390 GB HGDP + ~4 GB genome). Download time depends heavily on FTP server speed and your network; expect several hours.*

This section downloads the standard legacy bundle — hg38 genome, 1000 Genomes Phase 3 (~16 GB), and HGDP (~390 GB) — and builds all config and PAM files automatically via the built-in `setup` command. You run this once per machine, then skip straight to Section 3 for every search.

### 2a. Choose your working directory

```bash
# EDIT THIS: set the path where CRISPRme will store all data and results
CRISPRME_DIR="$HOME/my_crisprme_run"
```

### 2b. Run crisprme.py setup

With your `crisprme` environment active, run the setup command, passing your chosen directory with `--path`:

```bash
crisprme.py setup --path "$CRISPRME_DIR"
```

The command creates the directory structure if it does not exist and runs all data download and configuration steps automatically. It is resumable: if a connection drops, re-running the same command skips any files that already pass an MD5 integrity check. See Appendix A2 for step-by-step details on exactly what the setup command does.

**Testing setup with a single chromosome first (optional)**

To verify the full pipeline on a small download before committing to the full ~410 GB, pass `--chrom` with any canonical chromosome label (chr1–chr22 or chrX):

```bash
crisprme.py setup --chrom chr22 --path "$CRISPRME_DIR"
```

This downloads only the chr22 genome FASTA and VCF files for both 1000G and HGDP. Once you have confirmed the setup, re-run with `--chrom all` to fetch the full genome.

### If setup fails partway through

By default, each asset is checked with an MD5 digest before downloading. Simply re-running the same command will resume from where it left off, skipping anything that is already valid. To force a full re-download regardless of what is already present, add `--force`. 

### Common issues

- **FTP timeout**: re-run the setup command; incomplete files will be retried automatically.
- **`bgzip` not found**: ensure the `crisprme` environment is active (`mamba activate crisprme`).
- **Disk full**: check available space with `df -h .` before starting; the full bundle needs ~410 GB.

### 2c. What the setup command produces

When `setup` finishes, your working directory holds everything a search needs — `Genomes/hg38/`, `Annotations/`, `PAMs/`, the per-dataset variant folders `VCFs/hg38_1000G/` and `VCFs/hg38_HGDP/`, `samplesIDs/`, and two **combined config files that are generated for you and already span both default datasets**:

- `vcf.config.txt` — lists the variant folders to search; after the default setup it contains **both** `hg38_1000G` and `hg38_HGDP` (one per line), so a single `complete-search --vcf vcf.config.txt` searches 1000G **and** HGDP together and combines their results.
- `samplesIDs.config.txt` — the matching per-dataset sample-ID files, in the same order.

You never have to merge these by hand — that is what "combine everything" means here. To search only one dataset, point `--vcf`/`--samplesID` at a config listing just that folder (see Section 4, Step 3).

### 2d. Faster downloads via HuggingFace (optional)

`setup` pulls the genome and variants from their original UCSC/FTP sources, which can be slow. When the data has been published to a HuggingFace dataset repository, you can instead fetch it over HF's CDN with the `download` command — typically much faster, and resumable:

```bash
# reference bundle in one shot: genome + annotations + PAMs + sample-ID files
crisprme.py download --what all --path "$CRISPRME_DIR"

# add a variant dataset (delivered bgzipped, ready to search)
crisprme.py download --what vcf --dataset 1000G --path "$CRISPRME_DIR"
crisprme.py download --what vcf --dataset HGDP  --path "$CRISPRME_DIR"
```

The default repository is `lucapinello/crisprme-data`; override it with `--hf-repo <org/name>` or the `CRISPRME_HF_REPO` environment variable. FASTA is decompressed automatically after download; VCFs stay bgzipped. You can freely mix sources — e.g. `download` the genome but `setup` the variants, or vice versa.

---

## Section 3. Run a search

*This is the everyday workflow. Everything from here is fast and lightweight.*

> **Before running any search**, ensure the CRISPRme environment is active: `mamba activate crisprme`

### 3a. Choose a PAM file

The setup command pre-creates PAM files for 17 common nucleases in `PAMs/`. The examples in this guide use the standard SpCas9 file:

```
PAMs/20bp-NGG-SpCas9.txt
```

Contents:

```
NNNNNNNNNNNNNNNNNNNNNGG 3
```

The 20 leading `N`s define the spacer length; `NGG` is the PAM sequence; `3` is the PAM length in bp. To see all available files, run `ls PAMs/`. For a nuclease not covered, see Section 5 to create a custom PAM file.

### 3b. Write a guide file

Create a plain text file with one guide RNA sequence per line. Each sequence is the spacer followed by `N`s to pad out the PAM (one `N` per PAM nucleotide).

```bash
# EDIT THIS: replace the sequence with your guide
GUIDE_SEQ="CTAACAGTTGCTTTTATCACNNN"
GUIDE_FILE="sg1617_guide.txt"

printf '%s\n' "$GUIDE_SEQ" > "$GUIDE_FILE"
```

For multiple guides, add one sequence per line.

### 3c. Run the search

The setup command creates a combined config file (`vcf.config.txt`) that includes both 1000G and HGDP. Use it directly in the search command:

```bash
# EDIT THESE: set your guide file and output folder name
GUIDE_FILE="sg1617_guide.txt"
OUTPUT_NAME="sg1617_search_1"

crisprme.py complete-search \
  --genome Genomes/hg38 \
  --pam PAMs/20bp-NGG-SpCas9.txt \
  --guide "$GUIDE_FILE" \
  --vcf vcf.config.txt \
  --samplesID samplesIDs.config.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 \
  --bDNA 1 \
  --bRNA 1 \
  --merge 3 \
  --output "$OUTPUT_NAME" \
  --thread 16
```

Results will be written to `Results/$OUTPUT_NAME/`.

> Whole-genome searches take several hours. Run inside `tmux` or `screen` so the job continues if your connection drops.

To search against only one dataset (1000G or HGDP), create a single-line config file and pass it instead:

```bash
printf 'hg38_1000G\n'              > vcf.config.1000G.txt
printf 'hg38_1000G.samplesID.txt\n' > samplesIDs.config.1000G.txt
```

Then substitute `vcf.config.1000G.txt` and `samplesIDs.config.1000G.txt` for the `--vcf` and `--samplesID` flags.

### 3d. Run another search

The first search builds subfolders in `Genomes/`, `genome_library/`, and `Dictionaries/`. Some parameter changes (PAM, VCF, bulge count) will trigger rebuilds; see Appendix A1 for the full breakdown. To search a new guide, only `--guide` and `--output` need to change:

```bash
# EDIT THESE: set your guide file and output folder name
GUIDE_SEQ="GAGTCCGAGCAGAAGAAGAANNN"
GUIDE_FILE="sg_other_guide.txt"
OUTPUT_NAME="sg_other_search_1"

printf '%s\n' "$GUIDE_SEQ" > "$GUIDE_FILE"

# Rerun — no re-download, no re-index
crisprme.py complete-search \
  --genome Genomes/hg38 \
  --pam PAMs/20bp-NGG-SpCas9.txt \
  --guide "$GUIDE_FILE" \
  --vcf vcf.config.txt \
  --samplesID samplesIDs.config.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 \
  --bDNA 1 \
  --bRNA 1 \
  --merge 3 \
  --output "$OUTPUT_NAME" \
  --thread 16
```

### 3e. Directory structure after completing Section 3

This is exactly what you should see on disk after running both example searches above:

```
my_crisprme_run/
├── Annotations/
│   ├── dhs+encode+gencode.hg38.bed.gz
│   └── gencode.protein_coding.bed.gz
├── Dictionaries/                            # populated after first search, reused after
├── Genomes/                                 # populated further on first search
│   └── hg38/
│       ├── chr1.fa
│       ├── chr2.fa
│       └── ...                              # chr3–chrY + unplaced scaffolds
├── genome_library/                          # populated after first search, reused after
├── PAMs/
│   ├── 20bp-NGG-SpCas9.txt
│   ├── 20bp-NGC-SpCas9.txt
│   └── ...                                  # 17 PAM files total (see Section 3a)
├── Results/
│   ├── sg1617_search_1/
│   └── sg_other_search_1/
├── samplesIDs/
│   ├── hg38_1000G.samplesID.txt
│   └── hg38_HGDP.samplesID.txt
├── VCFs/
│   ├── hg38_1000G/
│   │   ├── ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
│   │   ├── ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
│   │   └── ...                              # chr3–chrX (23 files total)
│   └── hg38_HGDP/
│       ├── hgdp_wgs.20190516.full.chr1.vcf.gz
│       ├── hgdp_wgs.20190516.full.chr2.vcf.gz
│       └── ...                              # chr3–chrX (23 files total)
├── samplesIDs.config.txt
├── sg1617_guide.txt
├── sg_other_guide.txt
└── vcf.config.txt
```

---

## Section 3.5. Prebuild, reuse, and share the reference index

Bulge-enabled searches need a CRISPRitz **index** of the reference genome. By default `complete-search` builds it on the first run — into `genome_library/<PAM>_<bulges+1>_<genome>` — and silently reuses it on later runs from the same working directory (that is why "Run another search" in Section 3d is fast: *no re-index*). Three commands let you control this explicitly, which matters for batch jobs, shared clusters, and skipping the build entirely.

### Build the index once, ahead of time

`build-index-only` builds exactly that index **without** running a search — handy on a large machine, or before launching a batch of jobs:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 \
  --pam PAMs/20bp-NGG-SpCas9.txt \
  --bDNA 1 --bRNA 1 \
  --thread 16 \
  --path "$CRISPRME_DIR"
```

Pass the **same** `--genome`, `--pam`, `--bDNA`, `--bRNA` you will use in `complete-search`, so the folder name matches and the search reuses it. (Zero-bulge searches need no index and this command will say so.)

### Point a search at a prebuilt (or shared) index

`--index-path` tells `complete-search` to look for the reference index in a specific library instead of building one under the working directory. If a matching index isn't there, the run **fails fast with a clear message** (rather than silently rebuilding) — which is what you want on a read-only shared mount:

```bash
crisprme.py complete-search \
  ... \
  --index-path /shared/crisprme/genome_library
```

### Share indexes via HuggingFace

Publish a locally built index so other machines (or collaborators) can skip the build entirely:

```bash
# upload — needs an HF write token: export HF_TOKEN=...  (never commit it)
crisprme.py publish-index --index genome_library/NGG_2_hg38

# elsewhere: download it straight into genome_library/
crisprme.py download --what index --index-name NGG_2_hg38 --path "$CRISPRME_DIR"
```

Then run the search with `--index-path "$CRISPRME_DIR/genome_library"` (or simply from `$CRISPRME_DIR`) and it reuses the downloaded index. The index folder name encodes the PAM, bulge count and genome, so an index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`.

---

## Section 4. Add a new VCF dataset

*Use this section whenever you want to run searches against a cohort not included in the default setup. Adding any new dataset follows the same four steps every time. The 1000 Genomes 2021 high-coverage dataset (3,202 samples, phased, GRCh38/hg38) is used as a worked example throughout.*

Adding a new dataset requires four things: the per-chromosome VCF files, a samplesIDs metadata file, and two short config files that tell CRISPRme where to find them.

**Prerequisite.** This section assumes you have already set up the reference genome and annotations (Section 2, or `crisprme.py download --what all`). A new VCF dataset is searched *against* that existing `Genomes/hg38/` reference — you do not re-download the genome to add a cohort.

**File naming requirement.** CRISPRme identifies the chromosome of each VCF by splitting the filename on `.` and looking for a segment that starts with `chr`. Your filenames must follow this pattern — the chromosome label must be its own dot-separated segment, e.g. `MyCohort.chr1.vcf.gz` or `ALL.chr1.filtered.vcf.gz`. A filename like `MyCohort_chr1.vcf.gz` (chromosome embedded inside an underscore-delimited prefix) will not work.

---

### Step 1. Place your VCF files

Create a subdirectory under `VCFs/` named after your dataset. Place one bgzipped `.vcf.gz` file per chromosome inside it. Do not include `.tbi` index files.

```bash
mkdir -p VCFs/MyCohort
# copy or symlink your per-chromosome bgzipped VCF files here
# filenames must have the chromosome as a standalone dot-separated segment:
#   MyCohort.chr1.vcf.gz  ✓
#   MyCohort_chr1.vcf.gz  ✗  (chromosome embedded in prefix — will not work)
```

The subdirectory name (`MyCohort` here) becomes the dataset identifier used in config files and results output.

### Step 2. Place a samplesIDs file

Create a tab-separated file in `samplesIDs/` with a header line and four columns: `#SAMPLE_ID`, `POPULATION_ID`, `SUPERPOPULATION_ID`, `SEX` (values: `male` or `female`).

```
#SAMPLE_ID	POPULATION_ID	SUPERPOPULATION_ID	SEX
NA12878	CEU	EUR	female
HG00096	GBR	EUR	male
```

```bash
cp /path/to/my_samples_metadata.txt samplesIDs/samplesIDs.MyCohort.txt
```

### Step 3. Create config files

Two one-line config files tell CRISPRme which VCF directory and samplesIDs file to use:

```bash
printf 'MyCohort\n'                  > vcf.config.MyCohort.txt
printf 'samplesIDs.MyCohort.txt\n'   > samplesIDs.config.MyCohort.txt
```

The value in `vcf.config` must exactly match the subdirectory name under `VCFs/`.
The value in `samplesIDs.config` must exactly match the filename under `samplesIDs/`.

### Step 4. Run the search

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 \
  --pam PAMs/20bp-NGG-SpCas9.txt \
  --guide sg1617_guide.txt \
  --vcf vcf.config.MyCohort.txt \
  --samplesID samplesIDs.config.MyCohort.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 \
  --bDNA 1 \
  --bRNA 1 \
  --merge 3 \
  --output sg1617_MyCohort_search_1 \
  --thread 16
```

To search multiple datasets simultaneously, list each on its own line in a combined config file:

```bash
printf 'hg38_1000G\nMyCohort\n'                                      > vcf.config.combined.txt
printf 'hg38_1000G.samplesID.txt\nsamplesIDs.MyCohort.txt\n'         > samplesIDs.config.combined.txt
```

---

### Worked example: 1000 Genomes 2021 high-coverage dataset

The updated 1KG 2021 release (3,202 samples, eagle2/shapeit2-duohmm phased, GRCh38) is a direct upgrade over the legacy 1KG dataset included in the default setup.

> **Known issues in the current CRISPRme release — read before downloading.**
> The two issues below affect this specific dataset and are fixed in a forthcoming update. Workarounds are provided for each.
>
> **Issue 1 — chrX haploid genotypes.**
> The 2021 release uses eagle2 phasing for chrX, which stores hemizygous male genotypes as bare `0` or `1` (no `|` separator) rather than the diploid `0|0`/`0|1`/`1|1` encoding used for autosomes.
> **Workaround: download and use chr1–22 only. Do not include the chrX VCF.**
>
> **Issue 2 — raw FTP filenames do not follow CRISPRme naming conventions.**
> The FTP files are named `CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered...vcf.gz`. The chromosome is embedded inside the first underscore-delimited segment rather than appearing as its own dot-separated element.
> **Workaround: rename the files after downloading (Step 2 below).**

**Step 1 — Download the phased VCFs (chr1–22 only)**

*Approximate download size: ~20 GB. Run inside `tmux` or `screen`.*

```bash
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

mkdir -p VCFs/hg38_1000G_2021
for chrom in $CHROMS; do
    wget -c "${BASE_URL}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz" \
         -O "VCFs/hg38_1000G_2021/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
done
```

**Step 2 — Rename files for CRISPRme compatibility**

```bash
cd VCFs/hg38_1000G_2021
for f in CCDG_*.vcf.gz; do
    mv "$f" "${f/CCDG_14151_B01_GRM_WGS_2020-08-05_/1KG_2021.}"
done
cd ../..
```

Files are now named `1KG_2021.chr1.filtered.shapeit2-duohmm-phased.vcf.gz`, etc., with the chromosome as a standalone segment.

**Step 3 — Prepare the samplesIDs file**

```bash
wget -q -O /tmp/1kg_3202_ped.txt \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

# File is space-delimited: FamilyID SampleID FatherID MotherID Sex Population Superpopulation
# Sex: 1 = male, 2 = female
awk 'BEGIN{OFS="\t"} NR==1{print "#SAMPLE_ID","POPULATION_ID","SUPERPOPULATION_ID","SEX"; next}
     {sex=($5==1?"male":"female"); print $2,$6,$7,sex}' /tmp/1kg_3202_ped.txt \
    > samplesIDs/hg38_1000G_2021.samplesID.txt
```

Spot-check the result:

```bash
head -3 samplesIDs/hg38_1000G_2021.samplesID.txt
wc -l samplesIDs/hg38_1000G_2021.samplesID.txt   # expect 3203 lines (header + 3202 samples)
```

**Step 4 — Create config files and run the search**

```bash
printf 'hg38_1000G_2021\n'                    > vcf.config.1KG2021.txt
printf 'hg38_1000G_2021.samplesID.txt\n'      > samplesIDs.config.1KG2021.txt
```

```bash
# EDIT THESE: set your guide file and output folder name
GUIDE_FILE="sg1617_guide.txt"
OUTPUT_NAME="sg1617_1KG2021_search_1"

crisprme.py complete-search \
  --genome Genomes/hg38 \
  --pam PAMs/20bp-NGG-SpCas9.txt \
  --guide "$GUIDE_FILE" \
  --vcf vcf.config.1KG2021.txt \
  --samplesID samplesIDs.config.1KG2021.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 \
  --bDNA 1 \
  --bRNA 1 \
  --merge 3 \
  --output "$OUTPUT_NAME" \
  --thread 16
```

To search the legacy 1KG dataset (from Section 2) and the 2021 dataset together:

```bash
printf 'hg38_1000G\nhg38_1000G_2021\n'                                      > vcf.config.1KG_both.txt
printf 'hg38_1000G.samplesID.txt\nhg38_1000G_2021.samplesID.txt\n'          > samplesIDs.config.1KG_both.txt
```

---

## Section 5. Add a new PAM / nuclease

*Use this section when you want to search with a nuclease not covered by the pre-created PAM files in Section 3a.*

### 5a. Write the PAM file

PAM files define the spacer length and PAM sequence. The format is:

```
<spacer_Ns><PAM_sequence> <PAM_length_in_bp>
```

Example for SpCas9 (20 bp spacer, NGG PAM):

```bash
printf 'NNNNNNNNNNNNNNNNNNNNNGG 3\n' > PAMs/20bp-NGG-SpCas9.txt
```

The 20 leading `N`s define the spacer length; `NGG` is the PAM; `3` is the PAM length in bp. For 5' PAM nucleases (e.g. Cas12a), the PAM bases come first and the position index is negative (e.g. `TTTV` + 20 `N`s with index `-4`). See `docs/INPUT_FORMATS.md` for the full PAM-file specification, including the filename convention (`<length>-<motif>-<Cas>.txt`) and the one-motif-per-file rule.

> **Note (bulges + partially-degenerate odd-length PAMs).** A partially-degenerate IUPAC motif of **odd** length (e.g. `WTN`) combined with bulges can trigger a crash in older CRISPRitz engines. CRISPRme prints a non-fatal warning in this case; even-length degenerate motifs such as `TTTV` (Cas12a) are unaffected. The underlying issue is fixed in CRISPRitz ≥ 2.7.1 (bundled with current CRISPRme), so it is a caution rather than a blocker — but if you hit an odd-length degenerate PAM with bulges on an older stack, either add one base to make the motif even-length or reduce bulges to 0.

Once the reference index for a new PAM does not yet exist, the first search with it builds one automatically (see Section 3.5); you can also pre-build it with `crisprme.py build-index-only`.

### 5b. Run the search

No manual reindexing is needed — CRISPRme detects the new PAM at runtime and builds a fresh index automatically. The first search with a new PAM takes as long as a full run from scratch; subsequent searches with the same PAM reuse the index. See Appendix A1 for details on what gets cached and when.

```bash
# EDIT THESE
PAM_FILE="PAMs/20bp-NGG-SpCas9.txt"   # replace with your new PAM file
GUIDE_FILE="my_guide.txt"
OUTPUT_NAME="my_search_1"

crisprme.py complete-search \
  --genome Genomes/hg38 \
  --pam "$PAM_FILE" \
  --guide "$GUIDE_FILE" \
  --vcf vcf.config.txt \
  --samplesID samplesIDs.config.txt \
  --annotation Annotations/dhs+encode+gencode.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 4 \
  --bDNA 1 \
  --bRNA 1 \
  --merge 3 \
  --output "$OUTPUT_NAME" \
  --thread 16
```

---

## Appendix

### A1. What CRISPRme writes and when

CRISPRme writes to four directories during a search. Subfolders are created automatically on first use and reused on subsequent runs — nothing needs to be deleted manually. The tables below show exactly which parameters determine each subfolder name.

**`Genomes/`**

| Subfolder | Written when | Depends on |
|---|---|---|
| `{ref}/` | Section 2 setup only | `--genome` |
| `{ref}+{vcf}/` | First search with this ref + VCF | `--genome`, `--vcf` |
| `{ref}+{vcf}_INDELS/` | First search with this ref + VCF | `--genome`, `--vcf` |

**`genome_library/`**

| Subfolder | Written when | Depends on |
|---|---|---|
| `{pam}_{bulge}_{ref}/` | First search with this PAM + bulge + ref | `--pam`, `--bDNA`/`--bRNA`, `--genome` |
| `{pam}_{bulge}_{ref}+{vcf}/` | First search with this PAM + bulge + ref + VCF | `--pam`, `--bDNA`/`--bRNA`, `--genome`, `--vcf` |
| `{pam}_{bulge}_{ref}+{vcf}_INDELS/` | First search with this PAM + bulge + ref + VCF | `--pam`, `--bDNA`/`--bRNA`, `--genome`, `--vcf` |

Where `{pam}` is the PAM sequence extracted from the PAM file (e.g. `NGG` for SpCas9), and `{bulge}` = `max(--bDNA, --bRNA) + 1` (e.g. `--bDNA 1 --bRNA 1` → `2`). The reference-only index (`{pam}_{bulge}_{ref}/`) is shared across all VCF datasets. `--mm` does not affect any folder name and never triggers reindexing.

**`Dictionaries/`**

| Subfolder | Written when | Depends on |
|---|---|---|
| `dictionaries_{vcf}/` | First search with this VCF | `--vcf` |
| `log_indels_{vcf}/` | First search with this VCF | `--vcf` |

Dictionaries depend only on the VCF dataset name and are reused regardless of PAM, bulge count, or mismatch count.

**`Results/`**

| Subfolder | Written when | Depends on |
|---|---|---|
| `{name}/` | Every search | `--output` (must be a new, empty name) |

Everything else (`--guide`, `--mm`, `--annotation`, `--merge`, `--thread`, `--samplesID`) affects only what is written inside `Results/{name}/` and never triggers new folders in the other three directories.

### A2. What crisprme.py setup does, step by step

Each operation the setup command performs, for users who need to debug, customize, or run pieces manually.

### A2a. Create the working directory structure

```bash
mkdir -p Genomes VCFs Annotations PAMs Dictionaries genome_library Results samplesIDs
```

All subsequent commands must be run from inside this working directory. CRISPRme writes all output relative to it.

### A2b. Download the hg38 reference genome

**Full genome (default):** downloads the complete hg38 assembly from UCSC (~900 MB compressed) and unpacks it as one uncompressed FASTA per chromosome. The archive MD5 is verified before extraction to catch truncated downloads. Existing chromosomes that pass the integrity check are skipped on re-runs.

**Single chromosome (`--chrom chrN`):** downloads only `chrN.fa.gz` from the UCSC per-chromosome endpoint and unpacks it to `Genomes/chrN/chrN.fa`. Use this for testing before committing to the full download.

Expected result: `Genomes/hg38/` (full) or `Genomes/chrN/` (single-chrom) containing one `.fa` file per downloaded chromosome.

### A2c. Download variant VCFs (1000G and HGDP)

Both datasets are always downloaded together. For each chromosome in scope:

- **1000 Genomes Phase 3** — biallelic SNV and INDEL calls on GRCh38 from the EBI FTP (~16 GB total). Per-file MD5 verification; corrupt or incomplete files are deleted and re-downloaded on the next run.
- **HGDP** — Human Genome Diversity Project whole-genome sequencing from the Sanger FTP (~390 GB total). Same MD5 verification.

Expected result: `VCFs/hg38_1000G/` and `VCFs/hg38_HGDP/` each containing 23 `.vcf.gz` files (chr1–22 + chrX) for a full run, or one file each for a single-chromosome run.

### A2d. Download samplesIDs metadata

Downloads pre-built sample metadata files (sample ID → population → superpopulation → sex) for both datasets and stores them with standardised names:

- `samplesIDs/hg38_1000G.samplesID.txt`
- `samplesIDs/hg38_HGDP.samplesID.txt`

### A2e. Download and prepare annotation files

Downloads two annotation archives from the CRISPRme GitHub repository, extracts them, and bgzip-compresses the result:

- `Annotations/gencode.protein_coding.bed.gz` — GENCODE protein-coding gene annotations
- `Annotations/dhs+encode+gencode.hg38.bed.gz` — DHS + ENCODE + GENCODE functional annotation

Each file is checked with MD5 before downloading; if it already exists and passes, the download is skipped.

### A2f. Write config files

Writes two config files in the working directory that point to the downloaded datasets. These are always (re-)written on every setup run:

```
vcf.config.txt
    hg38_1000G
    hg38_HGDP

samplesIDs.config.txt
    hg38_1000G.samplesID.txt
    hg38_HGDP.samplesID.txt
```

### A2g. Write PAM files

Writes PAM definition files for all supported nucleases to `PAMs/`. Each file is skipped if it already exists. The full set written by setup is:

- `20bp-NGG-SpCas9.txt`
- `20bp-NGC-SpCas9.txt`
- `20bp-NGK-SpCas9.txt`
- `20bp-NRG-SpCas9.txt`
- `20bp-NRCH-SpCas9.txt`
- `20bp-NNGT-SpCas9.txt`
- `20bp-NAA-iSpyMacCas9.txt`
- `22bp-NNGRRN-SaCas9.txt`
- `TTCN-20bp-CasX.txt`
- `TTTV-20bp-Cas12a.txt`
- `TTTV-21bp-Cas12a.txt`
- `TTTV-23bp-Cas12a.txt`
- `TTTV-25bp-AsCas12a.txt`
- `18bp-NNN-NO_PAM.txt`
- `20bp-NNN-NO_PAM.txt`
- `21bp-NNN-NO_PAM.txt`
- `23bp-NNN-NO_PAM.txt`

---

For full documentation see the [CRISPRme GitHub repository](https://github.com/pinellolab/CRISPRme).
