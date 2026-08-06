# CRISPRme Docker Quickstart — run the web interface in a few commands

This is the fastest way to get CRISPRme running with its point-and-click **web
interface**, with **no conda, no compiling, and no 410 GB download**. You copy a
few commands, open a browser, fill in a short form, and read the results.

Everything runs inside Docker, so the only thing you install is Docker itself.

---

## 1. Install Docker (one time)

- **Mac / Windows:** install **Docker Desktop** — https://docs.docker.com/get-started/
- **Linux:** install **Docker Engine** — https://docs.docker.com/engine/install/

Then, in Docker Desktop → **Settings → Resources**, give Docker at least **16 GB
of memory** (more is better for whole-genome searches). Check it works:

```bash
docker run --rm hello-world
```

## 2. Make a folder to hold your data and results

Everything CRISPRme downloads and produces will live here (so it is kept between
runs). Pick any folder:

```bash
mkdir -p ~/crisprme && cd ~/crisprme
```

> In every command below, `-v "${PWD}:/DATA"` shares *this* folder with the
> container as `/DATA`. That is how your downloads and results are saved to your
> computer instead of disappearing when the container stops.

## 3. Download the reference data (one time, minutes — not hours)

This pulls the human genome, annotations, PAM files and sample lists from the
CRISPRme HuggingFace mirror (a fast CDN). It replaces the old multi-hour `setup`:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme \
  crisprme.py download --what all --path /DATA
```

## 4. Download a prebuilt index (one time, minutes — skips a long build)

Bulge-enabled searches need a genome **index**. Building it yourself takes ~10
minutes of CPU; instead, download the ready-made SpCas9 (NGG) index:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme \
  crisprme.py download --what index --index-name NGG_2_hg38 --path /DATA
```

The web interface picks this up automatically — a search that uses the NGG PAM
with up to 1 bulge will reuse it instead of rebuilding. (Need a different
nuclease? See **"Installing more indexes"** at the bottom.)

## 5. (Optional) Add 1000 Genomes variants for a variant-aware search

CRISPRme's superpower is finding off-targets created by genetic variants. To
enable that, also download the 1000 Genomes variant set (~16 GB):

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme \
  crisprme.py download --what vcf --dataset 1000G --path /DATA
```

You can skip this for a first run and do a plain reference-genome search.

## 6. Start the web interface

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA -p 8080:8080 -it \
  pinellolab/crisprme crisprme.py web-interface
```

`-p 8080:8080` connects the app inside the container to your browser. Leave this
running, and open:

**http://127.0.0.1:8080**

(If you started it on a remote server, use that server's address instead of
`127.0.0.1`.)

## 7. Run a search in the browser

1. Enter a **guide/spacer** sequence (e.g. `CTAACAGTTGCTTTTATCAC`).
2. Choose the **PAM** (e.g. `20bp-NGG-SpCas9`) and the **genome** (`hg38`).
3. Set **mismatches** and **bulges** (e.g. 4 / 1 / 1). If you downloaded the
   variants in step 5, select the **1000G** dataset to make it variant-aware.
4. Give the job a name and click **Submit**.
5. Watch the live status page; when it finishes, open the **Results** to explore
   the off-targets, scores and plots.

Your results are saved on your computer under `~/crisprme/Results/<job name>/`.

To stop the web server, press **Ctrl+C** in the terminal.

---

## Installing more indexes (as you need them)

An index is specific to a **PAM + bulge count + genome**. Download whichever you
need — for example a Cas12a (TTTV) index — the same way:

```bash
docker run --rm -v "${PWD}:/DATA" -w /DATA pinellolab/crisprme \
  crisprme.py download --what index --index-name TTTV_2_hg38 --path /DATA
```

To see which indexes are published, browse the dataset repository
(`lucapinello/crisprme-data`, folder `indexes/`). If an index you need is not
published yet, the web interface will simply build it on the first search that
uses that PAM (slower once, then reused). You can also pre-build it from the
command line with `crisprme.py build-index-only` — see
[`PRECOMPUTED_INDEXES.md`](PRECOMPUTED_INDEXES.md).

## Troubleshooting

- **The genome/PAM dropdowns are empty** → you skipped step 3. Run
  `crisprme.py download --what all --path /DATA` (inside Docker, as above).
- **Browser can't connect to `127.0.0.1:8080`** → make sure the `docker run`
  in step 6 includes `-p 8080:8080` and is still running.
- **Results disappeared after closing the container** → make sure every command
  includes `-v "${PWD}:/DATA"` so files are written to your folder, not the
  throwaway container.
- **A search ran out of memory** → raise Docker's memory limit (step 1), or
  start with a smaller search (fewer mismatches/bulges, or one chromosome).

For the full web-interface reference (every form field, all the result tabs),
see [`crisprme_web_interface_user_guide.md`](crisprme_web_interface_user_guide.md).
