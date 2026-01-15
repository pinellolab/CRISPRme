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
targets_tsv_ref_chrom="${targets_tsv_ref}.${chrom}"
targets_tsv_alt_chrom="${targets_tsv_alt}.${chrom}"
LC_ALL=C grep -F -w $chrom "$targets_tsv_ref" > "$targets_tsv_ref_chrom"
LC_ALL=C grep -F -w $chrom "$targets_tsv_alt" > "$targets_tsv_alt_chrom"

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
