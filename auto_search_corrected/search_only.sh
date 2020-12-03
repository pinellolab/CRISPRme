#!/bin/bash

#file for automated search of guide+pam in reference and variant genomes

ref_folder=$(realpath $1)
vcf_folder=$(realpath $2)
guide_file=$(realpath $3)
pam_file=$(realpath $4)
bMax=$5
mm=$6
bDNA=$7
bRNA=$8

output_folder=$(realpath $9)

starting_dir=${10}
ncpus=${11}
echo "CPU used: $ncpus" 
#echo $ref_folder
#echo $vcf_folder
#echo $guide_file
#echo $pam_file
#echo $annotation_file
#echo $output_folder
ref_name=$(basename $1)
#folder_of_folders=$(dirname $1)
vcf_name=$(basename $2)
guide_name=$(basename $3)
pam_name=$(basename $4)

if [ "$vcf_name" != "_" ]; then
	log="log_search_only_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.txt"
else
	log="log_search_only_${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.txt"
fi

touch $output_folder/$log

echo "Search Start:" $(date +%F-%T) >> $output_folder/$log
echo "########################################" >> $output_folder/$log
echo "INPUT:" >> $output_folder/$log
echo "Reference - $ref_name" >> $output_folder/$log
echo "VCFs - $vcf_name" >> $output_folder/$log
echo "Guide - $guide_name" >> $output_folder/$log
echo "PAM - $pam_name" >> $output_folder/$log
echo "########################################" >> $output_folder/$log

declare -a real_chroms
for file_chr in "$ref_folder"/*.fa
do
	file_name=$(basename $file_chr)
	chr=$(echo $file_name | cut -f 1 -d'.')
	echo "$chr"
	real_chroms+=("$chr")	
done

if [ "$vcf_name" != "_" ]; then
	declare -a array_fake_chroms
	for file_chr in "$vcf_folder"/*.vcf.gz
	do
		file_name=$(basename $file_chr)
		chr=$(echo $file_name | cut -f 2 -d'.')
		echo "fake$chr"
		array_fake_chroms+=("fake$chr")	
	done
fi

if ! [ -d "$output_folder" ]; then
	mkdir "$output_folder"
fi
cd "$output_folder/"

fullseqpam=$(cut -f1 -d' ' "$pam_file")
pos=$(cut -f2 -d' ' "$pam_file")
if [ $pos -gt 0 ]; then
	true_pam=${fullseqpam:${#fullseqpam}-$pos}
else
	true_pam=${fullseqpam:0:-$pos}
fi

if [ "$vcf_name" != "_" ]; then
	if ! [ -d "dictionaries_$vcf_name/" ]; then
		echo "Dictionaries Start: "$(date +%F-%T) >> $output_folder/$log
		mkdir "dictionaries_$vcf_name"
		cd "$starting_dir"
		./create_dict.py "$vcf_folder" "$vcf_name" "$output_folder" "$log" &
		pid_dicts=$!
		cd "$output_folder/"
	else
		echo "Dictionaries already present"
	fi
	dict_folder="$output_folder/dictionaries_$vcf_name/"
fi

if [ "$vcf_name" != "_" ]; then
	if ! [ -d "fake_chrom_$vcf_name" ]; then
		echo "INDELs Start: "$(date +%F-%T) >> $output_folder/$log
		echo "Generating fake chromosomes for indels"
		mkdir "fake_chrom_$vcf_name"
		cd "$starting_dir"
		./pool_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder/fake_chrom_$vcf_name" "$log" $ncpus &
		pid_indels=$!
		cd "$output_folder"
	else
		echo "INDELs extraction already present"
		# for key in "${!array_fake_chroms[@]}"
		# do
			# touch finished$key.txt
		# done
	fi
	
	if ! [ -d "${ref_name}+${vcf_name}" ]; then
		echo "Add-variants Start: "$(date +%F-%T) >> $output_folder/$log
		echo "Adding variants"
		crispritz.py add-variants "$vcf_folder/" "$ref_folder/"
		cp -r "variants_genome/SNPs_genome/${ref_name}_enriched/" "./${ref_name}+${vcf_name}"
		rm -r "variants_genome/"
		echo "Add-variants End: "$(date +%F-%T) >> $output_folder/$log
	else
		echo "Variants already added"
	fi
fi

while kill "-0" $pid_indels &>/dev/null; do
	echo "Waiting for indels process"
	sleep 600
done

if ! [ -d "genome_library/"$true_pam"_${bMax}_"$ref_name ]; then
	echo "Index-genome Reference Start: "$(date +%F-%T) >> $output_folder/$log	
	echo "Indexing reference genome"
	crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $(expr $ncpus - 2)
	pid_index_ref=$!
	echo "Index-genome Reference End: "$(date +%F-%T) >> $output_folder/$log	
else
	echo "Reference Index already present"
fi


if [ "$vcf_name" != "_" ]; then
	if ! [ -d "genome_library/"$true_pam"_${bMax}_${ref_name}+${vcf_name}" ]; then
		echo "Index-genome Variant Start: "$(date +%F-%T) >> $output_folder/$log	
		echo "Indexing variant genome"
		crispritz.py index-genome "${ref_name}+${vcf_name}" "${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $(expr $ncpus - 2) #${ref_folder%/}+${vcf_name}/
		pid_index_var=$!
		echo "Index-genome Variant Start: "$(date +%F-%T) >> $output_folder/$log	
	else
		echo "Variant Index already present"
	fi	
fi

if ! [ -d "$output_folder/crispritz_targets" ]; then
	mkdir "$output_folder/crispritz_targets"
fi

if ! [ -f "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	echo "Search Reference Start: "$(date +%F-%T) >> $output_folder/$log	
	crispritz.py search "genome_library/"$true_pam"_${bMax}_$ref_name/" "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th $(expr $ncpus / 4) &
	pid_search_ref=$!
else
	echo "Search for reference already done"
fi

if [ "$vcf_name" != "_" ]; then
	if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		echo "Search Variant Start: "$(date +%F-%T) >> $output_folder/$log	
		crispritz.py search "genome_library/"$true_pam"_${bMax}_${ref_name}+${vcf_name}/" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th $(expr $ncpus / 4) -var
		mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
		echo "Search Variant End: "$(date +%F-%T) >> $output_folder/$log	
	else
		echo "Search for variant already done"
	fi
	
	if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		echo "Search INDELs Start"
		echo "Search INDELs Start: "$(date +%F-%T) >> $output_folder/$log	
		cd $starting_dir
		./pool_search_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder" $true_pam $ncpus
		mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
		echo "Search INDELs End"
		echo "Search INDELs End: "$(date +%F-%T) >> $output_folder/$log
	fi
	
fi

while kill "-0" $pid_search_ref &>/dev/null; do
	echo "Waiting for search genome reference"
	sleep 300
done
if [ -f "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	mv "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
fi
echo "Search Reference End: "$(date +%F-%T) >> $output_folder/$log	

while kill "-0" $pid_dicts &>/dev/null; do
	echo "Waiting for dictionary creation before SNP analysis"
	sleep 300
done


if ! [ -d "$output_folder/crispritz_profiles" ]; then
	mkdir $output_folder/crispritz_profiles
fi
mv $output_folder/*profile* $output_folder/crispritz_profiles/ > /dev/null 2>&1

echo "Search End: "$(date +%F-%T) >> $output_folder/$log
echo "SEARCH END"

: '
if [ "$vcf_name" != "_" ]; then
	for key in "${array_fake_chroms[@]}"
	do
		tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >> "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
		header=$(head -1 "$output_folder/${key}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt")
		rm "$output_folder/${key}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
	done
	sed -i 1i"$header" "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
fi
'

