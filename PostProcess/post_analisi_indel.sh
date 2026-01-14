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

echo "Processing INDELs results for $key, starting post-analysis"

# define chromosome keys
chrom=$key
fakechrom="fake${chrom}"

# reference targets
targets_tsv_ref="${output_folder}/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_ref_chrom="${output_folder}/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.${chrom}"
awk -v fakechrom="$fakechrom" '$0 ~ fakechrom {print}' "$targets_tsv_ref" > "$targets_tsv_ref_chrom"  # subset to chrom-specific targets
# remove malformed lines, if any
awk -F'\t' 'NF >= 10' "$targets_tsv_ref_chrom" > "${targets_tsv_ref_chrom}.tmp"
mv "${targets_tsv_ref_chrom}.tmp" "$targets_tsv_ref_chrom"

# alternative targets
targets_tsv_alt="${output_folder}/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_alt_chrom="${output_folder}/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.${chrom}"
awk -v fakechrom="$fakechrom" '$0 ~ fakechrom {print}' "$targets_tsv_alt" > "$targets_tsv_alt_chrom"  # subset to chrom-specific targets
# remove malformed lines, if any
awk -F'\t' 'NF >= 10' "$targets_tsv_alt_chrom" > "${targets_tsv_alt_chrom}.tmp"
mv "${targets_tsv_alt_chrom}.tmp" "$targets_tsv_alt_chrom"

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
