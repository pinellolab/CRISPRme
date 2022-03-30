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
true_chr=$key
fake_chr="fake$true_chr"

LC_ALL=C grep -P "fake$key\t" "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
LC_ALL=C grep -P "fake$key\t" "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
header=$(head -1 $output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt)
sed -i 1i"$header" "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"

./analisi_indels_NNN.sh "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr" "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr" "$output_folder/${fake_chr}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}" "$annotation_file" "$dict_folder/log_indels_$vcf_name" "$ref_folder/$true_chr.fa" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$output_folder"
rm "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
rm "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
# tail -n +2 "$output_folder/${fake_chr}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
# tail -n +2 "$output_folder/${fake_chr}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res_alt" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
# rm "$output_folder/${fake_chr}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
# rm "$output_folder/${fake_chr}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
