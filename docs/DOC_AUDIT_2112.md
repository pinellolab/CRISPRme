# CRISPRme 2.1.12 — Documentation Audit (fresh-eyes new-user run)

Findings while following **only the README** to install + run from scratch on a
clean machine. Feeds the doc review (#21) and website update (#22).

## README.md issues

1. **Stale version `v2.1.6`** — "Expected output: `v2.1.6`" appears in §1.1.2 (Step 3),
   §3.1 Quick Test (twice). Should reflect 2.1.12 (or be version-agnostic).
2. **complete-search example uses 6/2/2** — §2.2.1 usage shows `--mm 6 --bDNA 2 --bRNA 2`,
   contradicting the 4/1/1 default adopted in 2.1.11+ (web example, complete-test,
   validate benchmarks all use 4/1/1). Update example to 4/1/1 (or note the default).
3. **Wrong annotation filename** — example: `Annotations/dhs+gencode+encode.hg38.bed`;
   actual downloaded/shipped file is `dhs+encode+gencode.hg38.bed` (encode+gencode order).
   Copy-pasting the example path fails.
4. **`--result_dir` vs `--results_dir`** — §2.2.6 generate-personal-card usage example
   uses `--result_dir`; the Input Arguments list uses `--results_dir`. One is wrong;
   reconcile with the code.
5. **Docker typo `-i i pinellolab/crisprme`** — double `-i` in §2.2.4 targets-integration
   and §2.2.6 generate-personal-card Docker examples.
6. **Bioconda channel config adds `defaults`** — §1.1.1 Step 2 does
   `mamba config --add channels defaults` + `channel_priority strict`. Modern bioconda
   guidance is conda-forge > bioconda and NOT adding `defaults`; the shown order can
   cause solver friction. Recommend the current bioconda 3-channel setup.
7. **Personal-card example dir `Results/sg1617.6.2.2`** — encodes old 6/2/2 params;
   with 4/1/1 the dir name differs. Cosmetic but reinforces #2.

## Scalability / requirements docs

S1. **README "32 GB minimum" is too low for whole-genome variant searches.** A measured
   genome-wide 1000G search (sg1617, 4/1/1) peaked at **~55.6 GB RAM** in the search
   stage. README §0 says "Minimum 32 GB / Recommended 64 GB+"; in practice 32 GB is
   insufficient for whole-genome + population VCFs — 64 GB is the real floor. Reword,
   and note that the heavy memory is the search stage (enrichment/index stay ≤13 GB).
   Also worth stating the *time* cost (enrichment ~12.5 h single-threaded) so users
   aren't surprised — and point them at precomputed indexes once available.

## Robustness bugs (code — candidates for 2.1.13)

R1. **`complete-test` full-genome download has no retry and no resume.** In
   `PostProcess/complete_test.py:download_vcf_data`, each VCF is fetched then md5-checked;
   on any mismatch it does `raise ValueError("Download for X failed")` and aborts the
   whole run. `PostProcess/utils.py:download()`/`http_download()` neither skip existing
   valid files nor retry. Consequence: on a slow/flaky source (EBI 1000G at ~1.4 MB/s,
   trans-Atlantic) a single dropped transfer kills a 3+ hour run, and re-running
   re-downloads everything from chr1. **Hit live: chr15 came down truncated → md5 fail →
   whole genome-wide run aborted after ~2.5 h / 15 of 23 chroms.**
   → Fix: wrap download+md5 in a retry loop (N attempts w/ backoff); skip files that
      already exist with a matching md5 (enables resume). Strongly reinforces the
      HF-mirror / precomputed-index plan (reliability, not just speed).

## Environment / install-experience caveats

8. **`mamba create -n crisprme ...` can land in a shared env dir** — on systems with
   custom `envs_dirs` (e.g. a lab `SHARED_SOFTWARE/envs` ahead of `~/.conda/envs`),
   the README's named-env command resolves the new env into shared space. Consider
   recommending a **prefix env** (`mamba create -p ./crisprme-env ...`) or a note to
   check `conda config --show envs_dirs` first. (Hit on ml007.)

## To verify during the run
- [ ] `mamba create -n crisprme python=3.8 crisprme -y` works as written (named env perms)
- [ ] `crisprme.py complete-test --vcf_dataset 1000G` (full genome) downloads all data
      cleanly from the documented sources and completes at 4/1/1
- [ ] version prints `2.1.12`; crispritz bundled is `2.7.1`
- [ ] validate-test passes across all chromosomes
- [ ] web-interface launches; all tabs + images render
