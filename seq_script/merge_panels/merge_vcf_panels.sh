#!/usr/bin/env bash
#
# merge_vcf_panels.sh — merge multiple population-VCF resources into a single
# per-chromosome panel for ONE CRISPRme enrichment/scan.
#
# WHY: CRISPRme enriches a reference genome with variants from a VCF and scans
# once. Running N resources (1000G, HGDP, AllOfUs, TOPMed, gnomAD, ...) as N
# separate datasets means N enrichments + N scans + an after-the-fact merge of
# results. Merging the resources into one panel first lets us enrich once and
# scan once — a large time saving as the number of resources grows.
#
# WHAT THIS PRESERVES (provenance — recoverable per variant):
#   INFO/AF          POOLED allele frequency, recomputed across ALL merged
#                    samples (bcftools +fill-tags: AC/AN over the union panel).
#                    This is what CRISPRme's enricher reports.
#   INFO/AF_<SRC>    each source's ORIGINAL AF, kept verbatim (missing "." where
#                    that source did not contain the variant).
#   INFO/SRC         comma-list of the sources that contain the variant, e.g.
#                    "1000G,HGDP" (shared) or "HGDP" (unique to HGDP). Derived
#                    from which AF_<SRC> fields are present.
#
# TRADE-OFF (documented, intentional): phasing is NOT preserved across sources.
# Cohorts are typically disjoint (no shared samples), so cross-source multi-
# variant haplotypes have no real carrier anyway; single-variant off-targets —
# the overwhelming majority — are fully correct. Phase within a source is lost
# once merged. The win is the single scan.
#
# DEPENDENCIES:
#   bcftools >= 1.18 with the +fill-tags plugin (BCFTOOLS_PLUGINS set if needed)
#   tabix / bgzip (htslib)
#   Downstream enrichment REQUIRES CRISPRitz PR #36 (enricher AF/FILTER
#   robustness): it reads INFO/AF by EXACT key, so it selects the pooled AF and
#   not AF_1000G/AFR_AF that appear earlier in the record; it also guards the
#   multiallelic AF-count that bcftools merge can produce.
#
# USAGE:
#   merge_vcf_panels.sh <out_vcf_dir> <chr>
#   Configure SOURCES below (name | dir | filename glob with {C} | rename_chr 0/1)
#   Produces <out_vcf_dir>/merged.<chr>.vcf.gz (+ .tbi) with AF / AF_<SRC> / SRC.
#
set -euo pipefail

OUTDIR=${1:?usage: merge_vcf_panels.sh <out_vcf_dir> <chr>}
CHR=${2:?usage: merge_vcf_panels.sh <out_vcf_dir> <chr>}
BCFTOOLS=${BCFTOOLS:-bcftools}
DATA=${CRISPRME_VCF_DATA:-/srv/local/crisprme/data/VCFs}

# name | directory | per-chr filename glob ({C} = chr token) | needs chr-prefix rename
SOURCES=(
  "1000G|$DATA/hg38_1000G|ALL.{C}.*.vcf.gz|1"
  "HGDP|$DATA/hg38_HGDP|*.{C}.vcf.gz|0"
  # add more, e.g.:
  # "AllOfUs|$DATA/hg38_AllOfUs|*.{C}.vcf.gz|0"
  # "gnomAD|$DATA/hg38_gnomAD|*.{C}.vcf.gz|0"   # sites-only: contributes positions + AF_gnomAD, no samples
)

mkdir -p "$OUTDIR"
OUT="$OUTDIR/merged.$CHR.vcf.gz"
[ -f "$OUT.tbi" ] && { echo "[skip] $CHR already merged"; exit 0; }

# Per-chromosome intermediates (renamed source copies + merged.noSRC) are multi-GB
# each; a whole-genome, N-way-parallel run writes tens–hundreds of GB of temp. The
# default mktemp location is /tmp (often a small root volume) and overruns it, so
# put WORK on the same (large) volume as the output — honor $TMPDIR if the caller
# points it at a big volume, else fall back to the output directory itself.
WORK=$(mktemp -d -p "${TMPDIR:-$OUTDIR}")
trap 'rm -rf "$WORK"' EXIT
# contig rename map (bare -> chr-prefixed), matches the hg38 reference convention
for c in $(seq 1 22) X Y MT; do echo "$c chr$c"; done > "$WORK/rename_chr.txt"

prepped=()
for entry in "${SOURCES[@]}"; do
  IFS='|' read -r name dir glob needs_rename <<<"$entry"
  vcf=$(ls ${dir}/${glob//\{C\}/$CHR} 2>/dev/null | head -1)
  [ -n "$vcf" ] || { echo "[warn] $name has no VCF for $CHR — skipping this source"; continue; }
  echo "[$name] $(basename "$vcf")"
  printf "INFO/AF AF_%s\n" "$name" > "$WORK/renameaf_$name.txt"
  pre="$WORK/$name.$CHR.vcf.gz"
  # optional contig rename -> rename INFO/AF to AF_<name>
  if [ "$needs_rename" = "1" ]; then
    "$BCFTOOLS" annotate --rename-chrs "$WORK/rename_chr.txt" "$vcf" -Ou \
      | "$BCFTOOLS" annotate --rename-annots "$WORK/renameaf_$name.txt" -Oz -o "$pre"
  else
    "$BCFTOOLS" annotate --rename-annots "$WORK/renameaf_$name.txt" "$vcf" -Oz -o "$pre"
  fi
  "$BCFTOOLS" index -t "$pre"
  prepped+=("$pre")
done
[ "${#prepped[@]}" -ge 1 ] || { echo "[error] no sources for $CHR"; exit 1; }

# merge -> recompute pooled INFO/AF over the union panel
"$BCFTOOLS" merge -m none "${prepped[@]}" -Ou \
  | "$BCFTOOLS" +fill-tags -Oz -o "$WORK/merged.noSRC.vcf.gz" -- -t AF
"$BCFTOOLS" index -t "$WORK/merged.noSRC.vcf.gz"

# derive INFO/SRC from which AF_<name> fields are present (sites-only = fast)
srcnames=(); for entry in "${SOURCES[@]}"; do srcnames+=("$(cut -d'|' -f1 <<<"$entry")"); done
qfmt='%CHROM\t%POS\t%REF\t%ALT'; for n in "${srcnames[@]}"; do qfmt="$qfmt\t%INFO/AF_$n"; done
printf '##INFO=<ID=SRC,Number=.,Type=String,Description="Source dataset(s) containing this variant">\n' > "$WORK/srchdr.txt"
"$BCFTOOLS" view -G "$WORK/merged.noSRC.vcf.gz" \
  | "$BCFTOOLS" query -f "$qfmt\n" \
  | awk -F'\t' -v names="${srcnames[*]}" 'BEGIN{n=split(names,NM," ")}
      { s=""; for(i=1;i<=n;i++){ if($(4+i)!="."){ s=s (s?",":"") NM[i] } } print $1"\t"$2"\t"$3"\t"$4"\t"s }' \
  | bgzip > "$WORK/src.tab.gz"
tabix -s1 -b2 -e2 "$WORK/src.tab.gz"
"$BCFTOOLS" annotate -a "$WORK/src.tab.gz" -c CHROM,POS,REF,ALT,INFO/SRC -h "$WORK/srchdr.txt" \
  "$WORK/merged.noSRC.vcf.gz" -Oz -o "$OUT"
"$BCFTOOLS" index -t "$OUT"
echo "[done] $CHR: $("$BCFTOOLS" index -n "$OUT") records, $("$BCFTOOLS" query -l "$OUT" | wc -l) samples -> $OUT"
