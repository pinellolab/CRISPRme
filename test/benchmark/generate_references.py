#!/usr/bin/env python3
"""
Turnkey driver: (re)generate the brute-force ground-truth reference(s) for every
benchmark in benchmarks.json, given a per-chromosome variant-enriched FASTA.

This is how new benchmark examples are produced/cached: add an entry to
benchmarks.json, run this driver on the enriched genome, host the resulting TSV,
and paste its printed MD5 back into the registry.

Usage:
    python generate_references.py --enriched-fasta chr22.enriched.fa --chrom chr22 \
        [--only cas12a_hbg] [--out-dir .]

Notes:
- The brute-force generator is a pure-Python exhaustive scanner (~minutes per Mb),
  so a full chromosome is a batch job (run on a server / in CI, not interactively).
- Run once per chromosome and concatenate the per-chromosome TSVs to build a
  genome-wide reference, then record the combined MD5 in benchmarks.json.
- The enriched FASTA MUST be the same genome CRISPRme searches (its add-variants
  output); a differently-enriched FASTA yields spurious hits in variant-dense
  regions.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enriched-fasta", required=True,
                    help="Per-chromosome variant-enriched FASTA (CRISPRme add-variants output)")
    ap.add_argument("--chrom", required=True, help="Chromosome name, e.g. chr22")
    ap.add_argument("--only", default=None, help="Generate only this benchmark 'name'")
    ap.add_argument("--out-dir", default=HERE, help="Directory for output TSVs")
    args = ap.parse_args()

    registry = json.load(open(os.path.join(HERE, "benchmarks.json")))
    th = registry["thresholds"]
    generator = os.path.join(HERE, "generate_brute_force.py")

    for bench in registry["benchmarks"]:
        if args.only and bench["name"] != args.only:
            continue
        out = os.path.join(args.out_dir, f"{bench['name']}.{args.chrom}.tsv")
        cmd = [
            sys.executable, generator,
            "--fasta", args.enriched_fasta,
            "--rna", bench["guide_search"],
            "--max-mismatches", str(th["mm"]),
            "--max-dna-gaps", str(th["bDNA"]),
            "--max-rna-gaps", str(th["bRNA"]),
            "--chrom", args.chrom,
            "--output", out,
        ]
        sys.stderr.write(f"[{bench['name']}] {bench['nuclease']} {bench['guide_search']} "
                         f"mm={th['mm']} bDNA={th['bDNA']} bRNA={th['bRNA']} -> {out}\n")
        rc = subprocess.call(cmd)
        if rc != 0:
            sys.stderr.write(f"[{bench['name']}] generation FAILED (exit {rc})\n")
            sys.exit(rc)
        n = sum(1 for _ in open(out)) - 1
        sys.stderr.write(f"[{bench['name']}] done: {n} rows, md5({os.path.basename(out)})={md5(out)}\n")


if __name__ == "__main__":
    main()
