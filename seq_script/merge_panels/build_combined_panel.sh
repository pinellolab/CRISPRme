#!/usr/bin/env bash
#
# build_combined_panel.sh — end-to-end: merge all chromosomes into one panel,
# assemble its samplesID, then enrich hg38 and build the (pamless) index.
#
# Resumable: per-chromosome merges that already produced <panel>/merged.<chr>.vcf.gz.tbi
# are skipped. Aborts before enrichment if the merge is incomplete (never builds a
# variant-less index by accident).
#
# See merge_vcf_panels.sh for the merge/provenance logic and dependencies (notably
# CRISPRitz PR #36 for enrichment). Configure SOURCES in that script.
#
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

CWD=${1:?usage: build_combined_panel.sh <crisprme_working_dir> [panel_name] [pam] [bdna] [brna] [parallel]}
PANEL=${2:-hg38_1000G_HGDP}
PAM=${3:-20bp-NNN-NO-PAM}          # pamless NNN by default (one index serves all PAMs)
BDNA=${4:-2}; BRNA=${5:-2}
PAR=${6:-4}
CHROMS=${CHROMS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"}
NCHR=$(wc -w <<<"$CHROMS")
VCFDIR="$CWD/VCFs/$PANEL"

echo "===== PHASE A: per-chromosome merge ($PAR-way parallel) [$(date +%T)] ====="
printf '%s\n' $CHROMS | xargs -P "$PAR" -I{} bash "$HERE/merge_vcf_panels.sh" "$VCFDIR" {}
NMERGED=$(ls "$VCFDIR"/merged.*.vcf.gz.tbi 2>/dev/null | wc -l)
echo "merged: $NMERGED / $NCHR"
[ "$NMERGED" -eq "$NCHR" ] || { echo "ABORT: merge incomplete — not enriching/indexing"; exit 1; }

echo "===== PHASE B: panel samplesID (union of per-source samplesIDs) [$(date +%T)] ====="
# concat each source's samplesID (one header). Edit for your source set.
SID="$CWD/samplesIDs/$PANEL.samplesID.txt"
DATA=${CRISPRME_SID_DATA:-/srv/local/crisprme/data/samplesIDs}
{ cat "$DATA/hg38_1000G.samplesID.txt"; tail -n +2 "$DATA/hg38_HGDP.samplesID.txt"; } > "$SID"
echo "samplesID rows: $(($(wc -l < "$SID") - 1))"

echo "===== PHASE C: enrich + build index [$(date +%T)] ====="
cd "$CWD"
crisprme.py build-index-only --genome "$CWD/Genomes/hg38" --pam "$CWD/PAMs/$PAM.txt" \
  --bDNA "$BDNA" --bRNA "$BRNA" --vcf "$VCFDIR" --path "$CWD"
echo "build-index-only exit=$? [$(date +%T)]"
echo "indexes: $(ls "$CWD/genome_library/" | grep -iE "$PANEL" | tr '\n' ' ')"
echo "BUILD_COMBINED_DONE [$(date +%T)]"
