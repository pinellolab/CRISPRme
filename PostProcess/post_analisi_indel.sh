#!/bin/bash

output_folder=$1
ref_folder=$2
ref_name=$(basename $2)
vcf_folder=$3
vcf_name=$(basename $3)
guide_file=$4
guide_name=$(basename $4)
mm=$5
bDNA=$6
bRNA=$7
annotation_file=$8
annotation_name=$(basename $8)
pam_file=$9
pam_name=$(basename $9)
# sampleID=${10}
dict_folder=${10}

final_res=${11}
final_res_alt=${12}

key=${13}

# process indel targets on chromosome $key
chrom=$key
fakechrom="fake${chrom}"

# reference and alternative crispritz targets files
targets_tsv_ref="$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_alt="${output_folder}/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"


echo "Processing INDELs results for $key, starting post-analysis"

# extract chrom-specific reference and alternative targets
#
# The indel search runs on the per-chromosome FAKE genome (see
# pool_search_indels.py), so every target row's Chromosome column is named
# "fake${chrom}" (e.g. "fakechr22"), NOT "${chrom}". We therefore MUST match on
# "$fakechrom". Matching on "$chrom" (as a broken refactor did) can never hit
# "fakechr22" with -w because the left word boundary fails ("chr22" is preceded
# by the word character "e"), so the subset comes back empty and every indel
# off-target is silently dropped with a clean exit (GitHub issue #172). -w still
# safely prevents "fakechr2" from matching "fakechr22" because the right
# boundary there is a digit.
targets_tsv_ref_chrom="${targets_tsv_ref}.${chrom}"
targets_tsv_alt_chrom="${targets_tsv_alt}.${chrom}"
LC_ALL=C grep -F -w "$fakechrom" "$targets_tsv_ref" > "$targets_tsv_ref_chrom"
LC_ALL=C grep -F -w "$fakechrom" "$targets_tsv_alt" > "$targets_tsv_alt_chrom"

# drop malformed lines, if any (a well-formed CRISPRitz target row has at least
# 10 tab-separated columns: Bulge_type..Total)
awk -F'\t' 'NF >= 10' "$targets_tsv_ref_chrom" > "${targets_tsv_ref_chrom}.tmp"
mv "${targets_tsv_ref_chrom}.tmp" "$targets_tsv_ref_chrom"
awk -F'\t' 'NF >= 10' "$targets_tsv_alt_chrom" > "${targets_tsv_alt_chrom}.tmp"
mv "${targets_tsv_alt_chrom}.tmp" "$targets_tsv_alt_chrom"

# Visible-but-safe guard: if a per-chrom subset is empty, warn loudly on STDOUT
# (never stderr: this pipeline treats any stderr write as fatal via
# `[ -s $logerror ]`). This is exactly the failure mode of issue #172.
if [ ! -s "$targets_tsv_alt_chrom" ]; then
    echo "WARNING: no indel targets matched '$fakechrom' in $targets_tsv_alt -- the per-chromosome indel subset for $chrom is EMPTY (see GitHub issue #172 if unexpected)."
fi

# adjust targets header
header=$(head -1 "$targets_tsv_alt")
sed -i 1i"$header" "$targets_tsv_alt_chrom"

# perform targets analysis by chromosome (scores, annotation, etc.)
targets_chrom_prefix="${output_folder}/${fakechrom}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}"
./analisi_indels_NNN.sh \
    "$targets_tsv_ref_chrom" \
    "$targets_tsv_alt_chrom" \
    "$targets_chrom_prefix" \
    "$annotation_file" \
    "${dict_folder}/log_indels_${vcf_name}" \
    "${ref_folder}/${chrom}.fa" \
    "$mm" "$bDNA" "$bRNA" \
    "$guide_file" \
    "$pam_file" \
    "$output_folder"

# remove chrom-specific targets tsv files
rm $targets_tsv_ref_chrom $targets_tsv_alt_chrom
