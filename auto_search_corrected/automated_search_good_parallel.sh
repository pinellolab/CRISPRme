#!/bin/bash

#file for automated search of guide+pam in reference and variant genomes

ref_folder=$(realpath $1)
vcf_folder=$(realpath $2)
guide_file=$(realpath $3)
pam_file=$(realpath $4)
annotation_file=$(realpath $5)
sampleID=$(realpath $6)

bMax=$7
mm=$8
bDNA=$9
bRNA=${10}

output_folder=$(realpath ${11})
starting_dir=${12}
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

declare -A pids_fake_chroms
for file_chr in "$ref_folder"/*.fa
do
	file_name=$(basename $file_chr)
	chr=$(echo $file_name | cut -f 1 -d'.')
	echo "$chr"
	pids_fake_chroms+=(["$chr"]="$!")	
done

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
		mkdir "dictionaries_$vcf_name"
		cd "$starting_dir"
		./create_dict.py "$vcf_folder" "$vcf_name" "$output_folder" &
		pid_dicts=$!
		cd "$output_folder/"
	else
		echo "Dictionaries already present"
	fi
	dict_folder="$output_folder/dictionaries_$vcf_name/"
fi

if ! [ -d "genome_library/"$true_pam"_${bMax}_"$ref_name ]; then
	echo "Indexing reference genome"
	crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th 2 &
	pid_index_ref=$!
else
	echo "Reference Index already present"
fi

if [ "$vcf_name" != "_" ]; then
	if ! [ -d "${ref_name}+${vcf_name}" ]; then
		echo "Adding variants"
		crispritz.py add-variants "$vcf_folder/" "$ref_folder/"
		cp -r "variants_genome/SNPs_genome/${ref_name}_enriched/" "./${ref_name}+${vcf_name}"
		rm -r "variants_genome/"
	else
		echo "Variants already added"
	fi
fi

while kill "-0" $pid_index_ref &>/dev/null; do
	echo "Waiting for index-genome for reference genome"
	sleep 10
done

if [ "$vcf_name" != "_" ]; then
	if ! [ -d "genome_library/"$true_pam"_${bMax}_${ref_name}+${vcf_name}" ]; then
		echo "Indexing variant genome"
		crispritz.py index-genome "${ref_name}+${vcf_name}" "${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th 2 & #${ref_folder%/}+${vcf_name}/
		pid_index_var=$!
	else
		echo "Variant Index already present"
	fi
	
	if ! [ -d "fake_chrom_$vcf_name" ]; then
		echo "Generating fake chromosomes for indels"
		mkdir "fake_chrom_$vcf_name"
		cd "$starting_dir"
		./pool_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder/fake_chrom_$vcf_name"
		cd "$output_folder"
	else
		echo "INDELs extraction already present"
		for key in "${!pids_fake_chroms[@]}"
		do
			touch finished$key.txt
		done
	fi
	
fi

while kill "-0" $pid_index_var &>/dev/null; do
	echo "Waiting for index-genome for variant genome"
	sleep 10
done

while kill "-0" $pid_dicts &>/dev/null; do
	echo "Waiting for dictionary creation before SNP analysis"
	sleep 10
done
	

if ! [ -f "${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	crispritz.py search "genome_library/"$true_pam"_${bMax}_$ref_name/" "$pam_file" "$guide_file" "${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th 2 &
	pid_search_ref=$!
else
	echo "Search for reference already done"
fi

if [ "$vcf_name" != "_" ]; then
	if ! [ -f "${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		crispritz.py search "genome_library/"$true_pam"_${bMax}_${ref_name}+${vcf_name}/" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th 2 -var
	else
		echo "Search for variant already done"
	fi
fi

while kill "-0" $pid_search_ref &>/dev/null; do
	echo "Waiting for search genome reference"
	sleep 10
done

cd "$starting_dir"

echo "Start post-analysis"

if [ "$vcf_name" != "_" ]; then
	final_res="$output_folder/final_results_${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt"
	final_res_alt="$output_folder/final_results_${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	fi

	for key in "${!pids_fake_chroms[@]}"
	do
		echo "Processing SNPs for $key"
		LC_ALL=C fgrep $key "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		LC_ALL=C fgrep $key "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		./scriptAnalisiNNN_v3.sh "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" "$annotation_file" "$dict_folder" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder"
		rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		header=$(head -1 "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt")
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt" >> "$final_res" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt" >> "$final_res_alt" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		#cat "$final_res" "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp" > "$final_res.tmp"
		#cat "$final_res_alt" "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp" > "$final_res_alt.tmp"
		#mv "$final_res.tmp" "$final_res"
		#mv "$final_res_alt.tmp" "$final_res_alt"
		rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt"
		rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt"
		#rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		#rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"

	done
else

	final_res="$output_folder/final_results_${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt"
	final_res_alt="$output_folder/final_results_${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	fi
	for key in "${!pids_fake_chroms[@]}"
	do
		echo "Processing $key"
		LC_ALL=C fgrep $key "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		touch "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		./scriptAnalisiNNN_v3.sh "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" "$annotation_file" "_" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder"
		rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		rm "$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		header=$(head -1 "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt")
		tail -n +2 "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt" >> "$final_res" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt" >> "$final_res_alt" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		#cat "$final_res" "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp" > "$final_res.tmp"
		#cat "$final_res_alt" "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp" > "$final_res_alt.tmp"
		#mv "$final_res.tmp" "$final_res"
		#mv "$final_res_alt.tmp" "$final_res_alt"
		rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt"
		rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt"
		#rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		#rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"

	done

fi

echo "Adding header to files"
sed -i 1i"$header" "$final_res"
sed -i 1i"$header" "$final_res_alt"

if [ "$vcf_name" != "_" ]; then
	echo "SNPs analysis ended. Starting INDELs analysis"
	cd "$starting_dir"
	while [ ${#pids_fake_chroms[@]} -ne 0 ];
	do
		for key in "${!pids_fake_chroms[@]}"
		do
			echo "Checking for INDELs results for $key"
			if ! [ -f "$output_folder/finished$key.txt" ]; then #kill "-0" ${pids_fake_chroms[$key]} &>/dev/null
				continue
			else
				echo "Found results for $key, starting post-analysis"
				true_chr=$key
				fake_chr="fake$true_chr"
				
				LC_ALL=C fgrep $true_chr "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
				./analisi_indels_NNN.sh "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr" "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}" "$annotation_file" "$output_folder/log_indels_$vcf_name" "$ref_folder/$true_chr.fa" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder" 
				rm "$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$true_chr"
				tail -n +2 "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt" >> "$final_res" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
				#cat "$final_res" "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp" > "$final_res.tmp"
				#mv "$final_res.tmp" "$final_res"
				tail -n +2 "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt" >> "$final_res_alt" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
				#cat "$final_res_alt" "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp" > "$final_res_alt.tmp"
				#mv "$final_res_alt.tmp" "$final_res_alt"
				rm "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt"
				rm "$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt"
				unset pids_fake_chroms[$key]
				rm "$output_folder/finished$key.txt"
			fi
		done
		sleep 10
	done 	
fi

cd "$starting_dir"
echo "Adjusting final results - sorting"
header=$(head -1 $final_res)
tail -n +2 "$final_res" | LC_ALL=C sort -k5,5 -k7,7n -o "$final_res.sorted"
sed -i 1i"$header" "$final_res.sorted"
mv "$final_res.sorted" "$final_res"
if [ "$vcf_name" != "_" ]; then
	header=$(head -1 $final_res_alt)
	tail -n +2 "$final_res_alt" | LC_ALL=C sort -k5,5 -k7,7n -o "$final_res_alt.sorted"
	sed -i 1i"$header" "$final_res_alt.sorted"
	mv "$final_res_alt.sorted" "$final_res_alt"
fi

echo "AUTOMATED SEARCH END"



: ' 
	declare -A pids_fake_chroms
	if ! [ -d "fake_chrom_$vcf_name" ]; then
		echo "Generating fake chromosomes for indels"
		mkdir "fake_chrom_$vcf_name"
		cd "$starting_dir"
		for file_chr in "$vcf_folder"/*
		do
			chr=$(echo $file_chr | cut -f 2 -d'.')
			./ indels_process.sh "$file_chr" $chr "$output_folder/fake_chrom_$vcf_name" "$vcf_name" "$ref_folder" "$pam_file" "$guide_file" $bMax $mm $bDNA $bRNA &
			pids_fake_chroms+=(["$chr"]="$!")	
		done
		#python3 indel_process.py "$vcf_folder" $vcf_name "$output_folder/fake_chrom_$vcf_name" "$ref_folder" "$pam_file" "$guide_file" $bMax $mm $bDNA $bRNA &
		cd "$output_folder/"
	fi
else
	declare -A pids_fake_chroms
	for file_chr in "$ref_folder"/*.fa
	do
		file_name=$(basename $file_chr)
		chr=$(echo $file_name | cut -f 1 -d'.')
		echo "$chr"
		pids_fake_chroms+=(["$chr"]="$!")	
	done
fi
'
