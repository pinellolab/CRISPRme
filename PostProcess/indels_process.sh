#!/bin/bash

vcf=$1
chr=$2
output_folder=$3
vcf_name=$4
ref_folder=$5
ref_name=$(basename $5)
pam_file=$6
pam_name=$(basename $6)
guide_file=$7
bMax=$8
mm=$9
bDNA=${10}
bRNA=${11}
use_thread=${12}

current_working_directory=${13}

starting_folder=$(pwd)
guide_name=$(basename "$guide_file")

fullseqpam=$(cut -f1 -d' ' "$pam_file")
pos=$(cut -f2 -d' ' "$pam_file")
if [ $pos -gt 0 ]; then
	true_pam=${fullseqpam:${#fullseqpam}-$pos}
else
	true_pam=${fullseqpam:0:-$pos}
fi

./indel_extraction.py "$vcf" "$output_folder/$chr" #it outputs a $chr.vcf.indels_only file 

#gzip $1'/'$vcf
cd "$output_folder"
echo 'Post processing the extracted vcf indels'
cat $chr.vcf.indels_only | sed 's/1|0/1/g' | sed 's/0|1/1/g' | sed 's/1|1/1/g' | sed 's/0|0/0/g' > $chr.vcf.processed
cut -f1,2,4,5 $chr.vcf.processed > $chr.indels_only_rownames # | sed 's/\t/,/g'
cut -f1,2,3,4,5,6,7,8,9 --complement $chr.vcf.processed > $chr.indels_only_matrix
rm $chr.vcf.processed

cd "$starting_folder"
./position_for_indels.py "$output_folder/$chr.indels_only_rownames" "$output_folder/$chr" #it ouputs a $chr.pos_indels

cd "$output_folder"
echo 'Retrieving a single fasta for all indels'
bedtools getfasta -fi "$ref_folder/$chr.fa" -bed $chr.pos_indels -fo $chr.fa.indels

awk ' NR % 2 == 0 { print; } NR % 2 == 1 {print $1"_"++c}' $chr.fa.indels | tr ':' '_' > $chr.fa.indels.ok
rm $chr.fa.indels

if [ ! -d "fake$chr" ]; then
        mkdir "fake$chr"
fi
#cd "$(dirname "$0")/$1"
echo 'Split the common fasta into single ones'
#semi=$1
#echo $semi
cd "fake$chr" 
awk '/^>/ {close(F) ; F = substr($1,2)".fa"} {print >> F}' "../$chr.fa.indels.ok"
rm "../$chr.fa.indels.ok"

#cd "$(dirname "$0")/.."
cd "$starting_folder"
./replace_indels.py "$output_folder/$chr.pos_indels" "$output_folder/$chr.indels_only_matrix" "$output_folder/$chr.indels_only_rownames" "$output_folder/fake$chr" $chr
rm "$output_folder/$chr.pos_indels"
rm "$output_folder/$chr.indels_only_rownames"
rm "$output_folder/$chr.indels_only_matrix"

cd "$output_folder"
if ! [ -d "../log_indels_$vcf_name" ]; then
	mkdir "../log_indels_$vcf_name"
fi
mv "fake$chr/log$chr.txt" "../log_indels_$vcf_name/"

cd "$starting_folder"
./integrate_logs.py "$output_folder/$chr.vcf.indels_only" "$output_folder/../log_indels_$vcf_name/log$chr.txt"
rm "$output_folder/$chr.vcf.indels_only"

cd "$output_folder/../log_indels_$vcf_name"
cat "log$chr.txt" | tr -d '>' > "log$chr.txt.tmp"
mv "log$chr.txt.tmp" "log$chr.txt"

cd "$current_working_directory/"
if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS/${true_pam}_2_fake$chr" ]; then
	echo "Indexing genome fake$chr"
	crispritz.py index-genome "fake${chr}_${ref_name}+${vcf_name}" "$output_folder/fake$chr/" "$pam_file" -bMax 2 -th 1 >/dev/null
	if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS" ]; then
		mkdir "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS"
	fi
	mv "$current_working_directory/genome_library/${true_pam}_2_fake${chr}_${ref_name}+${vcf_name}" "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS/${true_pam}_2_fake$chr"
else
	echo "Index already present for fake$chr"
fi
