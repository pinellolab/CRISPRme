#!/bin/bash

output_folder=$1
ref_folder=$2
ref_name=$(basename $2)
vcf_name=$3
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

echo "Processing SNPs for $key"

# reference targets
targets_tsv_ref="$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_ref_chrom="${targets_tsv_ref}.${key}"
awk -v key="$key" '$0 ~ key {print}' "$targets_tsv_ref" > "$targets_tsv_ref_chrom"  # subset to chrom-specific targets
# remove malformed lines, if any
awk -F'\t' 'NF >= 10' "$targets_tsv_ref_chrom" > "${targets_tsv_ref_chrom}.tmp"
mv "${targets_tsv_ref_chrom}.tmp" "$targets_tsv_ref_chrom"  

# alternative targets
targets_tsv_alt="$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_alt_chrom="${targets_tsv_alt}.${key}"
if [[ ! -f "$targets_tsv_alt" ]]; then
	cp "$targets_tsv_ref_chrom" "$targets_tsv_alt_chrom"  # to keep pipeline consistency
else
	awk -v key="$key" '$0 ~ key {print}' "$targets_tsv_alt" > "$targets_tsv_alt_chrom"  # subset to chrom-specific targets
	# remove malformed lines, if any
	awk -F'\t' 'NF >= 10' "$targets_tsv_alt_chrom" > "${targets_tsv_alt_chrom}.tmp"
    mv "${targets_tsv_alt_chrom}.tmp" "$targets_tsv_alt_chrom"
fi

# perform targets analysis by chromosome (scores, annotation, etc.)
targets_chrom_prefix="$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_${key}"
./scriptAnalisiNNN_v3.sh \
    "$targets_tsv_ref_chrom" \
    "$targets_tsv_alt_chrom" \
    "$targets_chrom_prefix" \
    "$annotation_file" \
    "${dict_folder}/my_dict_${key}.json" \
    "${ref_folder}/${key}.fa" \
    "$mm" "$bDNA" "$bRNA" \
    "$guide_file" \
    "$pam_file" \
    "$output_folder"

# remove chrom-specific targets tsv files
# rm $targets_tsv_ref_chrom $targets_tsv_alt_chrom  
