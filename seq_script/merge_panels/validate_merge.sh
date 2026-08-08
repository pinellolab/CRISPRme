#!/usr/bin/env bash
#
# validate_merge.sh — sanity-check a merged panel VCF's provenance + pooled AF.
# Usage: validate_merge.sh <merged.chr.vcf.gz> [source_name ...]
#   e.g. validate_merge.sh merged.chr22.vcf.gz 1000G HGDP
#
set -euo pipefail
BCFTOOLS=${BCFTOOLS:-bcftools}
M=${1:?usage: validate_merge.sh <merged.vcf.gz> [src ...]}
shift; SRCS=("$@"); [ "${#SRCS[@]}" -ge 1 ] || SRCS=(1000G HGDP)

echo "== samples: $("$BCFTOOLS" query -l "$M" | wc -l)   records: $("$BCFTOOLS" index -n "$M")"

echo "== INFO/AF header fields present:"
"$BCFTOOLS" view -h "$M" | grep -E "ID=AF(_|,)" | sed 's/^/   /'

echo "== SRC distribution (first 200k records):"
"$BCFTOOLS" query -f '%INFO/SRC\n' "$M" 2>/dev/null | head -200000 | sort | uniq -c | sort -rn | sed 's/^/   /'

echo "== pooled AF correctness — AC/AN must equal INFO/AF (5 sampled records):"
"$BCFTOOLS" +fill-tags "$M" -Ou -- -t AC,AN,AF 2>/dev/null \
  | "$BCFTOOLS" query -f '%POS AC=%INFO/AC AN=%INFO/AN AF=%INFO/AF\n' 2>/dev/null \
  | awk '{split($2,ac,"=");split($3,an,"=");split($4,af,"=");
          calc=(an[2]>0?ac[2]/an[2]:0);
          ok=(sprintf("%.4f",calc)==sprintf("%.4f",af[2])?"OK":"MISMATCH");
          print "   "$0" -> AC/AN="calc" ["ok"]"}' | head -5

echo "== per-source AF preserved + SRC consistent (5 records):"
qf='%POS SRC=%INFO/SRC AF=%INFO/AF'; for s in "${SRCS[@]}"; do qf="$qf AF_$s=%INFO/AF_$s"; done
"$BCFTOOLS" query -f "$qf\n" "$M" 2>/dev/null | head -5 | sed 's/^/   /'
