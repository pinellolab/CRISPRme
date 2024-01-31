#!/bin/bash

set -e # capture any failure

#file for automated search of guide+pam in reference and variant genomes

ref_folder=$(realpath $1)
vcf_list=$(realpath $2)
# IFS=',' read -ra vcf_list <<< $2
guide_file=$(realpath $3)
pam_file=$(realpath $4)
annotation_file=$(realpath $5)
sampleID=$(realpath $6)

bMax=$7
mm=$8
bDNA=$9
bRNA=${10}

merge_t=${11}

output_folder=$(realpath ${12})

starting_dir=$(realpath ${13})
ncpus=${14}
current_working_directory=$(realpath ${15})

gene_proximity=$(realpath ${16})

email=${17}
echo -e "MAIL: $email"
echo -e "CPU used: $ncpus"

#used to solve base editor check in resultintegration phase
base_check_start=${18}
base_check_end=${19}
base_check_set=${20}

# sorting criteria while merging best targets
sorting_criteria_scoring=${21}
sorting_criteria=${22}

log="$output_folder/log.txt"
touch $log
#echo -e 'Job\tStart\t'$(date) > $log
start_time='Job\tStart\t'$(date)

# output=$output_folder/output.txt
# touch $output
##CREATE DUMMY FILE WITH ONE LINE
echo -e "dummy_file" >"${output_folder}/.dummy.txt"
dummy_file="${output_folder}/.dummy.txt"
##CREATE EMPTY FILE
touch "${output_folder}/.empty.txt"
empty_file="${output_folder}/.empty.txt"
##CREATE EMPTY DIR
mkdir "${output_folder}/.empty"
empty_dir="${output_folder}/.empty"

rm -f $output_folder/queue.txt
#for vcf_f in "${vcf_list[@]}";
if [ $2 == "_" ]; then
	echo -e "_" >>$output_folder/tmp_list_vcf.txt
	vcf_list=$output_folder/tmp_list_vcf.txt
fi
echo >>$vcf_list
if [ $6 != "_" ]; then
	echo >>$6
fi

while read vcf_f; do
	if [ -z "$vcf_f" ]; then
		continue
	fi
	vcf_name+=$vcf_f"+"
	vcf_folder="${current_working_directory}/VCFs/${vcf_f}"
	ref_name=$(basename $1)
	#folder_of_folders=$(dirname $1)
	vcf_name=$(basename $vcf_f)
	echo "STARTING ANALYSIS FOR $vcf_name"
	# echo $vcf_name
	guide_name=$(basename $3)
	pam_name=$(basename $4)
	annotation_name=$(basename $5)

	echo -e $start_time >$log
	# echo -e 'Job\tStart\t'$(date) > $log
	# echo -e 'Job\tStart\t'$(date) >&2

	unset real_chroms
	declare -a real_chroms
	for file_chr in "$ref_folder"/*.fa; do
		file_name=$(basename $file_chr)
		# filename=$(basename -- "$fullfile")
		# extension="${filename##*.}"
		chr="${file_name%.*}"
		# chr=$(echo -e $file_name | cut -f 1 -d'.')
		echo -e "$chr"
		real_chroms+=("$chr")
	done

	if [ "$vcf_name" != "_" ]; then
		unset array_fake_chroms
		declare -a array_fake_chroms
		for file_chr in "$vcf_folder"/*.vcf.gz; do
			file_name=$(basename $file_chr)
			# file_name=$(basename $file_chr)
			IFS='.' read -ra ADDR <<<$file_name
			for i in "${ADDR[@]}"; do
				if [[ $i == *"chr"* ]]; then
					chr=$i
				fi
			done
			# chr=$(echo -e $file_name | cut -f 2 -d'.')
			echo -e "fake$chr"
			array_fake_chroms+=("fake$chr")
		done
	fi

	if ! [ -d "$output_folder" ]; then
		mkdir "$output_folder"
	fi

	fullseqpam=$(cut -f1 -d' ' "$pam_file")
	pos=$(cut -f2 -d' ' "$pam_file")
	if [ $pos -gt 0 ]; then
		true_pam=${fullseqpam:${#fullseqpam}-$pos}
	else
		true_pam=${fullseqpam:0:-$pos}
	fi

	# if ! [ -d "$current_working_directory/Results" ]; then
	# 	mkdir "$current_working_directory/Results"
	# fi

	if ! [ -d "$current_working_directory/Dictionaries" ]; then
		mkdir "$current_working_directory/Dictionaries"
	fi

	if ! [ -d "$current_working_directory/Genomes" ]; then
		mkdir "$current_working_directory/Genomes"
	fi

	if ! [ -d "$current_working_directory/genome_library/" ]; then
		mkdir "$current_working_directory/genome_library"
	fi

	if ! [ -d "$output_folder/crispritz_targets" ]; then
		mkdir "$output_folder/crispritz_targets"
	fi

	if [ "$vcf_name" != "_" ]; then

		cd "$current_working_directory/Genomes"
		if ! [ -d "$current_working_directory/Genomes/${ref_name}+${vcf_name}" ]; then
			echo -e 'Add-variants\tStart\t'$(date) >>$log
			# echo -e 'Add-variants\tStart\t'$(date) >&2
			echo -e "Adding variants"
			crispritz.py add-variants "$vcf_folder/" "$ref_folder/" "true" || {
				echo "CRISPRme ERROR: genome enrichment failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			#if ! [ -d "${ref_name}+${vcf_name}" ]; then
			#	mkdir "${ref_name}+${vcf_name}"
			#fi
			mv "$current_working_directory/Genomes/variants_genome/SNPs_genome/${ref_name}_enriched/" "./${ref_name}+${vcf_name}/"
			if ! [ -d "$current_working_directory/Dictionaries/dictionaries_${vcf_name}/" ]; then
				mkdir "$current_working_directory/Dictionaries/dictionaries_${vcf_name}/"
			fi
			if ! [ -d "$current_working_directory/Dictionaries/log_indels_${vcf_name}/" ]; then
				mkdir "$current_working_directory/Dictionaries/log_indels_${vcf_name}/"
			fi
			mv $current_working_directory/Genomes/variants_genome/SNPs_genome/*.json $current_working_directory/Dictionaries/dictionaries_${vcf_name}/
			mv $current_working_directory/Genomes/variants_genome/SNPs_genome/log*.txt $current_working_directory/Dictionaries/log_indels_${vcf_name}/
			cd "$current_working_directory/"
			if ! [ -d "genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS" ]; then
				mkdir "genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS"
			fi

			echo -e 'Add-variants\tEnd\t'$(date) >>$log
			# echo -e 'Add-variants\tEnd\t'$(date) >&2
			echo -e 'Indexing Indels\tStart\t'$(date) >>$log
			# echo -e 'Indexing Indels\tStart\t'$(date) >&2
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/variants_genome/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus || {
				echo "CRISPRme ERROR: indels indexing failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			echo -e 'Indexing Indels\tEnd\t'$(date) >>$log
			# echo -e 'Indexing Indels\tEnd\t'$(date) >&2
			if ! [ -d $current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS ]; then
				mkdir $current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS
			fi
			mv $current_working_directory/Genomes/variants_genome/fake* $current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS
			rm -r "$current_working_directory/Genomes/variants_genome/"
			dict_folder="$current_working_directory/Dictionaries/dictionaries_$vcf_name/"
		else
			echo -e "Variants already added"
			dict_folder="$current_working_directory/Dictionaries/dictionaries_$vcf_name/"
		fi
	fi

	cd "$current_working_directory/"
	if [ "$vcf_name" != "_" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS" ]; then
			echo -e 'Indexing Indels\tStart\t'$(date) >>$log
			# echo -e 'Indexing Indels\tStart\t'$(date) >&2
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus || {
				echo "CRISPRme ERROR: indels indexing failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			echo -e 'Indexing Indels\tEnd\t'$(date) >>$log
			# echo -e 'Indexing Indels\tEnd\t'$(date) >&2
		fi
	fi

	if [ -d "$current_working_directory/Dictionaries/fake_chrom_$vcf_name" ]; then
		rm -r "$current_working_directory/Dictionaries/fake_chrom_$vcf_name"
	fi

	cd "$current_working_directory/"
	if ! [ -d "$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}" ]; then
			if ! [ $bMax -gt 1 ]; then
				if ! [ -d "$current_working_directory/genome_library/${true_pam}_1_${ref_name}" ]; then
					echo -e 'Index-genome Reference\tStart\t'$(date) >>$log
					# echo -e 'Index-genome Reference\tStart\t'$(date) >&2
					# echo -e 'Indexing_Reference' > $output
					echo -e "Indexing reference genome"
					crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus || {
						echo "CRISPRme ERROR: TST-index construction failed (script: ${0} line $((LINENO - 1)))" >&2
						exit 1
					}
					pid_index_ref=$!
					echo -e 'Index-genome Reference\tEnd\t'$(date) >>$log
					# echo -e 'Index-genome Reference\tEnd\t'$(date) >&2
					idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
				else
					echo -e "Reference Index already present"
					idx_ref="$current_working_directory/genome_library/${true_pam}_1_${ref_name}"
				fi
			else
				echo -e 'Index-genome Reference\tStart\t'$(date) >>$log
				# echo -e 'Index-genome Reference\tStart\t'$(date) >&2
				# echo -e 'Indexing_Reference' > $output
				echo -e "Indexing reference genome"
				crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus || {
					echo "CRISPRme ERROR: TST-index construction failed (script: ${0} line $((LINENO - 1)))" >&2
					exit 1
				}
				pid_index_ref=$!
				echo -e 'Index-genome Reference\tEnd\t'$(date) >>$log
				# echo -e 'Index-genome Reference\tEnd\t'$(date) >&2
				idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
			fi
		else
			echo -e "Reference Index already present"
			echo -e 'Index-genome Reference\tEnd\t'$(date) >>$log
			# echo -e 'Index-genome Reference\tEnd\t'$(date) >&2
			idx_ref="$current_working_directory/genome_library/${true_pam}_2_${ref_name}"
		fi
	else
		echo -e "Reference Index already present"
		echo -e 'Index-genome Reference\tEnd\t'$(date) >>$log
		# echo -e 'Index-genome Reference\tEnd\t'$(date) >&2
		idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
	fi

	if [ "$vcf_name" != "_" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}" ]; then
			if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}" ]; then
				if ! [ $bMax -gt 1 ]; then
					if ! [ -d "$current_working_directory/genome_library/${true_pam}_1_${ref_name}+${vcf_name}" ]; then
						echo -e 'Index-genome Variant\tStart\t'$(date) >>$log
						# echo -e 'Index-genome Variant\tStart\t'$(date) >&2
						# echo -e 'Indexing_Enriched' > $output
						echo -e "Indexing variant genome"
						crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus || {
							echo "CRISPRme ERROR: TST-index construction failed (script: ${0} line $((LINENO - 1)))" >&2
							exit 1
						}
						pid_index_var=$!
						echo -e 'Index-genome Variant\tEnd\t'$(date) >>$log
						# echo -e 'Index-genome Variant\tEnd\t'$(date) >&2
						idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
					else
						echo -e "Variant Index already present"
						idx_var="$current_working_directory/genome_library/${true_pam}_1_${ref_name}+${vcf_name}"
					fi
				else
					echo -e 'Index-genome Variant\tStart\t'$(date) >>$log
					# echo -e 'Index-genome Variant\tStart\t'$(date) >&2
					# echo -e 'Indexing_Enriched' > $output
					echo -e "Indexing variant genome"
					crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus || {
						echo "CRISPRme ERROR: TST-index construction failed (script: ${0} line $((LINENO - 1)))" >&2
						exit 1
					}
					pid_index_ref=$!
					echo -e 'Index-genome Variant\tEnd\t'$(date) >>$log
					# echo -e 'Index-genome Variant\tEnd\t'$(date) >&2
					idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
				fi
			else
				echo -e "Variant Index already present"
				echo -e 'Index-genome Variant\tEnd\t'$(date) >>$log
				# echo -e 'Index-genome Variant\tEnd\t'$(date) >&2
				idx_var="$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}"
			fi
		else
			echo -e "Variant Index already present"
			echo -e 'Index-genome Variant\tEnd\t'$(date) >>$log
			# echo -e 'Index-genome Variant\tEnd\t'$(date) >&2
			idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
		fi
	fi

	#ceil npcus to use 1/2 of cpus per search
	ceiling_result=$((($ncpus) / 2))
	#if ceiling is 0, set ceiling to 1
	if [ $ceiling_result -eq 0 ]; then
		ceiling_result=1
	fi

	#start searches
	cd "$output_folder"
	echo $idx_ref
	#TODO ricerca ref lanciata in parallelo con ricerca alternative
	if ! [ -f "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		echo -e 'Search Reference\tStart\t'$(date) >>$log
		# echo -e 'Search Reference\tStart\t'$(date) >&2
		# echo -e 'Search Reference' >  $output
		if [ "$bDNA" -ne 0 ] || [ "$bRNA" -ne 0 ]; then
			crispritz.py search $idx_ref "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t -th $ceiling_result &
			wait || {
				echo "CRISPRme ERROR: off-targets search failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			} # TODO:
			pid_search_ref=$!
		else
			crispritz.py search "$current_working_directory/Genomes/${ref_name}/" "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
			wait || {
				echo "CRISPRme ERROR: off-targets search failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			pid_search_ref=$!
		fi
	else
		echo -e "Search for reference already done"
	fi

	if [ "$vcf_name" != "_" ]; then
		#TODO RICERCA ALTERNATIVE PARALLELA A REF
		if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e 'Search Variant\tStart\t'$(date) >>$log
			# echo -e 'Search Variant\tStart\t'$(date) >&2
			# echo -e 'Search Variant' >  $output
			if [ "$bDNA" -ne 0 ] || [ "$bRNA" -ne 0 ]; then
				crispritz.py search "$idx_var" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t -th $ceiling_result -var || {
					echo "CRISPRme ERROR: off-targets search failed (script: ${0} line $((LINENO - 1)))" >&2
					exit 1
				}
				# mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
				echo -e 'Search Variant\tEnd\t'$(date) >>$log
			else
				crispritz.py search "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
				wait || {
					echo "CRISPRme ERROR: off-targets search failed (script: ${0} line $((LINENO - 1)))" >&2
					exit 1
				}
				echo -e 'Search Variant\tEnd\t'$(date) >>$log
				# mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			fi
		else
			echo -e "Search for variant already done"
		fi

		if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e "Search INDELs Start"
			echo -e 'Search INDELs\tStart\t'$(date) >>$log
			# echo -e 'Search INDELs\tStart\t'$(date) >&2
			cd $starting_dir
			#commented to avoid indels search
			#TODO REMOVE POOL SCRIPT FROM PROCESSING
			./pool_search_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder" $true_pam "$current_working_directory/" "$ncpus" || {
				echo "CRISPRme ERROR: off-targets search on indels failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			# mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			awk '($3 !~ "n") {print $0}' "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp" || {
				echo "CRISPRme ERROR: off-targets report construction failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp" "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
			echo -e "Search INDELs End"
			echo -e 'Search INDELs\tEnd\t'$(date) >>$log
			# echo -e 'Search INDELs\tEnd\t'$(date) >&2
		else
			echo -e "Search INDELs already done"
		fi
	fi

	while kill "-0" $pid_search_ref &>/dev/null; do
		echo -e "Waiting for search genome reference"
		sleep 100
	done
	echo -e 'Search Reference\tEnd\t'$(date) >>$log
	# echo -e 'Search Reference\tEnd\t'$(date) >&2

	# move all targets into targets directory
	mv $output_folder/*.targets.txt $output_folder/crispritz_targets

	if ! [ -d "$output_folder/crispritz_prof" ]; then
		mkdir $output_folder/crispritz_prof
	fi
	mv $output_folder/*profile* $output_folder/crispritz_prof/ &>/dev/null

	cd "$starting_dir"

	echo -e "Start post-analysis"

	# echo -e 'Post analysis' >  $output
	if [ "$vcf_name" != "_" ]; then
		echo -e 'Post-analysis SNPs\tStart\t'$(date) >>$log
		# echo -e 'Post-analysis SNPs\tStart\t'$(date) >&2
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then
			touch "$final_res_alt"
		fi

		#TODO ANALISI DEGLI SNP IN PARALLELO
		./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $dict_folder $final_res $final_res_alt $ncpus || {
			echo "CRISPRme ERROR: SNP analysis failed (script: ${0} line $((LINENO - 1)))" >&2
			exit 1
		}

		#CONCATENATE REF&VAR RESULTS
		for key in "${real_chroms[@]}"; do
			echo "Concatenating $key"
			#touch file to avoid inconsistencies when files are broken or deleted
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt"
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt"
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt"
			#concatenate all files into respective final best file and remove them after computation
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt" >>"$final_res.bestCFD.txt"
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt" >>"$final_res.bestmmblg.txt"
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt" >>"$final_res.bestCRISTA.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt"
		done

		echo -e 'Post-analysis SNPs\tEnd\t'$(date) >>$log

	else
		echo -e 'Post-analysis\tStart\t'$(date) >>$log
		# echo -e 'Post-analysis\tStart\t'$(date) >&2
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then
			touch "$final_res_alt"
		fi

		./pool_post_analisi_snp.py $output_folder $ref_folder "_" $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "_" $final_res $final_res_alt $ncpus || {
			echo "CRISPRme ERROR: SNP analysis failed (script: ${0} line $((LINENO - 1)))" >&2
			exit 1
		}

		#CONCATENATE REF&VAR RESULTS
		for key in "${real_chroms[@]}"; do
			echo "Concatenating $key"
			#touch file to avoid inconsistencies when files are broken or deleted
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt"
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt"
			touch "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt"
			#concatenate all files into respective final best file and remove them after computation
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt" >>"$final_res.bestCFD.txt"
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt" >>"$final_res.bestmmblg.txt"
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt" >>"$final_res.bestCRISTA.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCFD.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestmmblg.txt"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestCRISTA.txt"
		done
		echo -e 'Post-analysis\tEnd\t'$(date) >>$log

	fi

	if [ "$vcf_name" != "_" ]; then
		echo -e "SNPs analysis ended. Starting INDELs analysis"
		cd "$starting_dir"

		echo -e 'Post-analysis INDELs\tStart\t'$(date) >>$log
		#SKIP INDELS ANALYSIS IF NO RESULTS FOUND
		if [ $(wc -l <"$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt") -gt 1 ]; then

			./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "$current_working_directory/Dictionaries/" $final_res $final_res_alt $ncpus || {
				echo "CRISPRme ERROR: indels analysis failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}

			#CONCATENATE INDELS RESULTS
			for key in "${array_fake_chroms[@]}"; do
				echo "Concatenating $key"
				#MERGE BEST INDEL TARGETS
				#create file if non-existent to avoid errors in tail processing
				touch "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCFD_INDEL.txt"
				touch "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCRISTA_INDEL.txt"
				touch "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestmmblg_INDEL.txt"
				tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCFD_INDEL.txt" >>"$final_res.bestCFD.txt"       #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
				tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCRISTA_INDEL.txt" >>"$final_res.bestCRISTA.txt" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
				tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestmmblg_INDEL.txt" >>"$final_res.bestmmblg.txt"
				#rm BEST INDEL TARGETS
				rm -f "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCFD_INDEL.txt"
				rm -f "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestCRISTA_INDEL.txt"
				rm -f "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestmmblg_INDEL.txt"
			done
		fi
		echo -e 'Post-analysis INDELs\tEnd\t'$(date) >>$log

	fi
done <$vcf_list

echo -e "Adding header to files"

while read samples; do
	if [ -z "$samples" ]; then
		continue
	fi
	# tail -n +2 $samples >> "$output_folder/.sampleID.txt"
	grep -v '#' "${current_working_directory}/samplesIDs/$samples" >>"$output_folder/.sampleID.txt"
done <$sampleID
# if [ "$vcf_name" != "_" ]; then
touch "$output_folder/.sampleID.txt"
sed -i 1i"#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX" "$output_folder/.sampleID.txt" || {
	echo "CRISPRme ERROR: Samples report construction failed (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
# fi

sampleID=$output_folder/.sampleID.txt

# echo -e 'Merging targets' >  $output

#create result file for each scoring method
# echo "header" >$final_res.bestCFD.txt
# echo "header" >$final_res.bestmmblg.txt
# echo "header" >$final_res.bestCRISTA.txt

#header into final_res best
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCFD.txt"
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestmmblg.txt"
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCRISTA.txt"
#header into final_res alt
echo "header" >$final_res_alt.bestCFD.txt
echo "header" >$final_res_alt.bestmmblg.txt
echo "header" >$final_res_alt.bestCRISTA.txt
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestCFD.txt"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestmmblg.txt"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestCRISTA.txt"

echo -e 'Merging Targets\tStart\t'$(date) >>$log
#SORT FILE TO HAVE CHR AND POS IN PROXIMITY TO MERGE THEM
#sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCFD.txt >$final_res.tmp
tail -n +2 $final_res.bestCFD.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCFD.txt
#sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCRISTA.txt >$final_res.tmp
tail -n +2 $final_res.bestCRISTA.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCRISTA.txt
#sort using guide_seq,chr,cluster_pos,total(mm+bul)
head -1 $final_res.bestmmblg.txt >$final_res.tmp
tail -n +2 $final_res.bestmmblg.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestmmblg.txt

# cp $final_res.bestCFD.txt $final_res.sorted.bestCFD.txt
#MERGE BEST FILES TARGETS TO REMOVE CONTIGOUS
#TODO CHECK MERGE
#SCORE CFD
./merge_close_targets_cfd.sh $final_res.bestCFD.txt $final_res.bestCFD.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
#TOTAL (MM+BUL)
./merge_close_targets_cfd.sh $final_res.bestmmblg.txt $final_res.bestmmblg.txt.trimmed $merge_t 'total' $sorting_criteria_scoring $sorting_criteria &
#SCORE CRISTA
./merge_close_targets_cfd.sh $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
wait
#CHANGE NAME TO BEST AND ALT FILES
mv $final_res.bestCFD.txt.trimmed $final_res.bestCFD.txt
mv $final_res.bestCFD.txt.trimmed.discarded_samples $final_res_alt.bestCFD.txt
mv $final_res.bestmmblg.txt.trimmed $final_res.bestmmblg.txt
mv $final_res.bestmmblg.txt.trimmed.discarded_samples $final_res_alt.bestmmblg.txt
mv $final_res.bestCRISTA.txt.trimmed $final_res.bestCRISTA.txt
mv $final_res.bestCRISTA.txt.trimmed.discarded_samples $final_res_alt.bestCRISTA.txt

#sort ALT files to avoid as much as possible duplicates
#REMOVED SINCE CAN RETURN DIFFERENT DIMENSION FILES
# sort -T $output_folder -u $final_res_alt.bestCFD.txt -o $final_res_alt.bestCFD.txt
# sort -T $output_folder -u $final_res_alt.bestmmblg.txt -o $final_res_alt.bestmmblg.txt
# sort -T $output_folder -u $final_res_alt.bestCRISTA.txt -o $final_res_alt.bestCRISTA.txt

echo -e 'Merging Targets\tEnd\t'$(date) >>$log

echo -e 'Annotating results\tStart\t'$(date) >>$log

#ANNOTATE BEST TARGETS
#TODO SISTEMARE ANNOTAZIONE (DIVISIONE INTERVAL TREE / PARALLEL SEARCH)
./annotate_final_results.py $final_res.bestCFD.txt $annotation_file $final_res.bestCFD.txt.annotated &
wait || {
	echo "CRISPRme ERROR: CFD annotation failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./annotate_final_results.py $final_res.bestmmblg.txt $annotation_file $final_res.bestmmblg.txt.annotated &
wait || {
	echo "CRISPRme ERROR: CRISTA annotation failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./annotate_final_results.py $final_res.bestCRISTA.txt $annotation_file $final_res.bestCRISTA.txt.annotated &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges annotation failed - reference(script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait
mv $final_res.bestCFD.txt.annotated $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.annotated $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.annotated $final_res.bestCRISTA.txt
#ANNOTATE ALT TARGETS
./annotate_final_results.py $final_res_alt.bestCFD.txt $annotation_file $final_res_alt.bestCFD.txt.annotated &
wait || {
	echo "CRISPRme ERROR: CFD annotation failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./annotate_final_results.py $final_res_alt.bestmmblg.txt $annotation_file $final_res_alt.bestmmblg.txt.annotated &
wait || {
	echo "CRISPRme ERROR: CRISTA annotation failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./annotate_final_results.py $final_res_alt.bestCRISTA.txt $annotation_file $final_res_alt.bestCRISTA.txt.annotated &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges annotation failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait
mv $final_res_alt.bestCFD.txt.annotated $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.annotated $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.annotated $final_res_alt.bestCRISTA.txt

#SCORING BEST RESULTS
./add_risk_score.py $final_res.bestCFD.txt $final_res.bestCFD.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: CFD risk score analysis failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./add_risk_score.py $final_res.bestmmblg.txt $final_res.bestmmblg.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: CRISTA risk score analysis failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./add_risk_score.py $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges risk score analysis failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait
mv $final_res.bestCFD.txt.risk $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.risk $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.risk $final_res.bestCRISTA.txt
#SCORING ALT RESULTS
./add_risk_score.py $final_res_alt.bestCFD.txt $final_res_alt.bestCFD.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: CFD risk score analysis failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./add_risk_score.py $final_res_alt.bestmmblg.txt $final_res_alt.bestmmblg.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: CRISTA risk score analysis failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
./add_risk_score.py $final_res_alt.bestCRISTA.txt $final_res_alt.bestCRISTA.txt.risk "False" &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges risk score analysis failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait
mv $final_res_alt.bestCFD.txt.risk $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.risk $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.risk $final_res_alt.bestCRISTA.txt

#remove N's and dots from rsID from BEST FILES
python remove_n_and_dots.py $final_res.bestCFD.txt &
wait || {
	echo "CRISPRme ERROR: CFD reports cleaning failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
python remove_n_and_dots.py $final_res.bestmmblg.txt &
wait || {
	echo "CRISPRme ERROR: CRISTA reports cleaning failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
python remove_n_and_dots.py $final_res.bestCRISTA.txt &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges reports cleaning failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait
#remove N's and dots from rsID from ALT FILES
python remove_n_and_dots.py $final_res_alt.bestCFD.txt &
wait || {
	echo "CRISPRme ERROR: CFD reports cleaning failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
python remove_n_and_dots.py $final_res_alt.bestmmblg.txt &
wait || {
	echo "CRISPRme ERROR: CRISTA reports cleaning failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
python remove_n_and_dots.py $final_res_alt.bestCRISTA.txt &
wait || {
	echo "CRISPRme ERROR: mismatch+bulges reports cleaning failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
wait

#join targets by columns for BEST and ALT files
pr -m -t -J $final_res.bestCFD.txt $final_res.bestmmblg.txt $final_res.bestCRISTA.txt >$final_res || {
	echo "CRISPRme ERROR: final report generation failed - reference (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
pr -m -t -J $final_res_alt.bestCFD.txt $final_res_alt.bestmmblg.txt $final_res_alt.bestCRISTA.txt >$final_res_alt || {
	echo "CRISPRme ERROR: final report generation failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}

#MERGE ALTERNATIVE CHR IF SAME SEQUENCE OF ALIGNED CHR
# ./merge_alt_chr.sh $final_res $final_res.chr_merged
# mv $final_res.chr_merged $final_res

#update header for final_res and final_res_alt
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res_alt"

echo -e 'Annotating results\tEnd\t'$(date) >>$log

# echo -e 'Creating images' >  $output
echo -e 'Creating images\tStart\t'$(date) >>$log

cd $output_folder
#FIX FILES NAMES AND REMOVE UNUSED FILES
echo -e "Cleaning directory"
rm -f *.CFDGraph.txt
rm -f indels.CFDGraph.txt
rm -r "crispritz_prof"
# rm -r "crispritz_targets" #remove targets in online version to avoid memory saturation ##CHECK THIS TO REMOVE ORIGINAL TARGETS
#change name to best and alt files
mv $final_res "${output_folder}/$(basename ${output_folder}).bestMerge.txt"
mv $final_res_alt "${output_folder}/$(basename ${output_folder}).altMerge.txt"

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	# ./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}" "var"
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "CFD" || {
		echo "CRISPRme ERROR: CFD report summary failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "fewest" || {
		echo "CRISPRme ERROR: mismatch+bulges report summary failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "CRISTA" || {
		echo "CRISPRme ERROR: CRISTA report summary failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
else
	# ./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}" "ref"
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "CFD" || {
		echo "CRISPRme ERROR: CFD report summary failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "fewest" || {
		echo "CRISPRme ERROR: mismatch+bulges report summary failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "CRISTA" || {
		echo "CRISPRme ERROR: CRISTA report summary failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
fi

if ! [ -d "$output_folder/imgs" ]; then
	mkdir "$output_folder/imgs"
fi

if [ "$vcf_name" != "_" ]; then
	cd "$output_folder/imgs"
	while IFS= read -r line || [ -n "$line" ]; do
		for total in $(seq 0 $(expr $mm + $bMax)); do
			# python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution.txt" $total $line
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CFD.txt" $total $line "CFD" || {
				echo "CRISPRme ERROR: CFD population report failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CRISTA.txt" $total $line "CRISTA" || {
				echo "CRISPRme ERROR: CRISTA population report failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_fewest.txt" $total $line "fewest" || {
				echo "CRISPRme ERROR: mismatch+bulges population report failed (script: ${0} line $((LINENO - 1)))" >&2
				exit 1
			}
		done

	done <$guide_file
fi

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	# ./radar_chart_dict_generator.py $guide_file "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CFD" || {
		echo "CRISPRme ERROR: CFD radar chart report failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA" || {
		echo "CRISPRme ERROR: CRISTA radar chart report failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "fewest" || {
		echo "CRISPRme ERROR: mismatch+bulges radar chart report failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
else
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "CFD" || {
		echo "CRISPRme ERROR: CFD radar chart report failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA" || {
		echo "CRISPRme ERROR: CRISTA radar chart report failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "fewest" || {
		echo "CRISPRme ERROR: mismatch+bulges radar chart report failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
fi
echo -e 'Creating images\tEnd\t'$(date) >>$log

echo $gene_proximity
echo -e 'Integrating results\tStart\t'$(date) >>$log
echo >>$guide_file

if [ $gene_proximity != "_" ]; then
	genome_version=$(echo ${ref_name} | sed 's/_ref//' | sed -e 's/\n//') #${output_folder}/Params.txt | awk '{print $2}' | sed 's/_ref//' | sed -e 's/\n//')
	echo $genome_version
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${gene_proximity}" $empty_file "${guide_file}" $genome_version "${output_folder}" $empty_dir $starting_dir/ $base_check_start $base_check_end $base_check_set || {
		echo "CRISPRme ERROR: postprocessing failed - reference (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).altMerge.txt" "${gene_proximity}" $empty_file "${guide_file}" $genome_version "${output_folder}" $empty_dir $starting_dir/ $base_check_start $base_check_end $base_check_set || {
		echo "CRISPRme ERROR: postprocessing failed - alternative (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	python $starting_dir/CRISPRme_plots.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/imgs/" &>"${output_folder}/warnings.txt" || {
		echo "CRISPRme ERROR: plots generation failed (script: ${0} line $((LINENO - 1)))" >&2
		exit 1
	}
	rm -f "${output_folder}/warnings.txt" #delete warnings file

fi
echo -e 'Integrating results\tEnd\t'$(date) >>$log
truncate -s -1 $guide_file
truncate -s -1 $vcf_list
if [ $6 != "_" ]; then
	truncate -s -1 $6
fi

echo -e 'Building database'
echo -e 'Creating database\tStart\t'$(date) >>$log
# echo -e 'Creating database\tStart\t'$(date) >&2
if [ -f "${output_folder}/$(basename ${output_folder}).db" ]; then
	rm -f "${output_folder}/$(basename ${output_folder}).db"
fi
#python $starting_dir/db_creation.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${output_folder}/$(basename ${output_folder})"
python $starting_dir/db_creation.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/.$(basename ${output_folder})" || {
	echo "CRISPRme ERROR: database creation failed (script: ${0} line $((LINENO - 1)))" >&2
	exit 1
}
echo -e 'Creating database\tEnd\t'$(date) >>$log
# echo -e 'Creating database\tEnd\t'$(date) >&2

# python $starting_dir/change_headers_bestMerge.py "${output_folder}/$(basename ${output_folder}).altMerge.txt" "${output_folder}/$(basename ${output_folder}).altMerge.new.header.txt"
# mv "${output_folder}/$(basename ${output_folder}).altMerge.new.header.txt" "${output_folder}/$(basename ${output_folder}).altMerge.txt"
#hide bestmerge file
# mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${output_folder}/.$(basename ${output_folder}).bestMerge.txt"

# echo -e 'Integrating results\tEnd\t'$(date) >&2
echo -e 'Job\tDone\t'$(date) >>$log
# echo -e 'Job\tDone\t'$(date) >&2
# echo -e 'Job End' >  $output

#change name of empirical not found
# mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.empirical_not_found.tsv" "${output_folder}/.$(basename ${output_folder}).bestMerge.txt.empirical_not_found.tsv"

if [ $(wc -l <"$guide_file") -gt 1 ]; then
	mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	mv "${output_folder}/$(basename ${output_folder}).altMerge.txt.integrated_results.tsv" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	#generate zipped version for file
	zip -j "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.zip" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	zip -j "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
else
	guide_elem=$(head -1 $guide_file)
	mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	mv "${output_folder}/$(basename ${output_folder}).altMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	#generate zipped version for file
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
fi
echo -e "JOB END"

if [ "$email" != "_" ]; then
	python $starting_dir/../pages/send_mail.py $output_folder
fi

#keep log_error but no block visualization
mv $output_folder/log_error.txt $output_folder/log_error_no_check.txt
#removing single best files after use and clean merged file to save space
#keep the two integrated files with all the targets
#save these files to test
rm $final_res.bestCFD.txt
rm $final_res.bestmmblg.txt
rm $final_res.bestCRISTA.txt
rm $final_res_alt.bestCFD.txt
rm $final_res_alt.bestmmblg.txt
rm $final_res_alt.bestCRISTA.txt
#save bestMerge and altMerge
# rm "${output_folder}/$(basename ${output_folder}).bestMerge.txt"
# rm "${output_folder}/$(basename ${output_folder}).altMerge.txt"
