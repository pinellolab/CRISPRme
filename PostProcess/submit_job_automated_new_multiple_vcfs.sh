#!/bin/bash

# This script automates the search for guide RNA and PAM sequences in reference 
# and variant genomes. It processes input files, manages directories, and executes 
# various analyses to identify potential CRISPR targets.
# The script handles both reference and variant genomes, performs indexing, 
# searching, and post-analysis, and generates results in specified output formats.
#
# Args:
#   $1: Reference genome folder path.
#   $2: List of VCF files.
#   $3: Guide RNA file path.
#   $4: PAM file path.
#   $5: Annotation file path.
#   $6: Sample ID file path.
#   $7: Maximum bulge size.
#   $8: Mismatches allowed.
#   $9: Bulge DNA size.
#   ${10}: Bulge RNA size.
#   ${11}: Merge threshold.
#   ${12}: Output folder path.
#   ${13}: Starting directory path.
#   ${14}: Number of CPUs to use.
#   ${15}: Current working directory path.
#   ${16}: Gene proximity file path.
#   ${17}: Email address for notifications.
#
# Returns:
#   Generates various output files including target lists, merged results, and a database.
#   Sends an email notification upon completion if an email address is provided.

ref_folder=$(realpath $1)  # reference genome folder
vcf_list=$(realpath $2)  # vcf folders list
guide_file=$(realpath $3)  # guide 
pam_file=$(realpath $4)  # pam
annotation_file=$(realpath $5)  # annotation bed
sampleID=$(realpath $6)  # sample ids 
bMax=$7  # max number of bulges
bMax_=$((bMax+1))  # superset of current maximum bulge number
mm=$8  # mismatches
bDNA=$9  # dna bulges
bRNA=${10}  # rna bulges 
merge_t=${11}  # targets merge threshold (bp)
output_folder=$(realpath ${12})  # output folder
starting_dir=$(realpath ${13})  # root dir 
ncpus=${14}  # number of threads
current_working_directory=$(realpath ${15})  # current working directory
gene_proximity=$(realpath ${16})  # gene annotation bed
email=${17}  # email address (website only)

echo -e "MAIL: $email"
echo -e "CPU used: $ncpus"

# used to solve base editor check in result-integration phase
base_check_start=${18}
base_check_end=${19}
base_check_set=${20}

# sorting criteria while merging best targets
sorting_criteria_scoring=${21}
sorting_criteria=${22}

# CI/CD test mode
cicd_test=${23}

# log files
log="$output_folder/log.txt"
touch $log
logerror="$output_folder/log_error.txt"
start_time='Job\tStart\t'$(date)

# check if there are input variant datasets
rm -f $output_folder/queue.txt
if [ $2 == "_" ]; then
	echo -e "_" >>$output_folder/tmp_list_vcf.txt
	vcf_list=$output_folder/tmp_list_vcf.txt
fi
echo >>$vcf_list
if [ $6 != "_" ]; then
	echo >>$6
fi

# iterate over each variant dataset to compute the enriched genomes, index each
# enriched genome, search off-targets, reconstruct haplotypes (if input VCFs are
# phased), and score off-targets
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

	# create fake chormosomes fasta for indels
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

	# create output folder in /Results/
	if ! [ -d "$output_folder" ]; then
		mkdir "$output_folder"
	fi

	# retrieve full pam sequence (guide + pam) and the position where the pam 
	# occurs within the full sequence (guide + pam)
	fullseqpam=$(cut -f1 -d' ' "$pam_file")
	pos=$(cut -f2 -d' ' "$pam_file")
	if [ $pos -gt 0 ]; then
		true_pam=${fullseqpam:${#fullseqpam}-$pos}
	else
		true_pam=${fullseqpam:0:-$pos}
	fi

	# ensure that fundamental folders are available 
	if ! [ -d "$current_working_directory/Results" ]; then
		mkdir "$current_working_directory/Results"
	fi
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

	# START STEP 1 - Genome enrichment
	if [ "$vcf_name" != "_" ]; then
		enriched_folder="$current_working_directory/Genomes/${ref_name}+${vcf_name}"
		variants_tmp="$current_working_directory/Genomes/variants_genome"
		dict_folder="$current_working_directory/Dictionaries/dictionaries_${vcf_name}/"
		indel_dict_folder="$current_working_directory/Dictionaries/log_indels_${vcf_name}/"
		indels_out="$current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS"

		if ! [ -d "$enriched_folder" ]; then
			echo -e 'Add-variants\tStart\t'$(date) >>$log
			echo -e "Adding variants"
			
			# enrich genome using crispritz
			cd "$current_working_directory/Genomes"
			crispritz.py add-variants "$vcf_folder/" "$ref_folder/" "true" 
			if [ -s $logerror ]; then
				printf "ERROR: Genome enrichment failed on %s\n" "$vcf_name" >&2
				# since failure happened, force genome enrichment to be repeated
				rm -r "$enriched_folder"* "$variants_tmp"
				exit 1
			fi

			# create indels folder
			mkdir -p $indels_out

			# move enriched snp genome and indels
			mv "$variants_tmp/SNPs_genome/${ref_name}_enriched/" "$enriched_folder/"
			mv $variants_tmp/fake* $indels_out

			# create dictionaries folders if needed
			mkdir -p "$dict_folder" "$indel_dict_folder"

			# move variant annotation files
			mv $variants_tmp/SNPs_genome/*.json "$dict_folder"
			mv $variants_tmp/SNPs_genome/log*.txt "$indel_dict_folder"

			# remove temporary variant genome folder
			rm -r "$variants_tmp"

			echo -e 'Add-variants\tEnd\t'$(date) >>"$log"
		else
			echo -e "Variants already added"
		fi

		cd $current_working_directory
	fi

	# check if anything odd happened during enrichment
	if [ -s $logerror ]; then
		printf "ERROR: Genome enrichment failed on %s\n" "$vcf_name" >&2
		exit 1
	fi
	# END STEP 1 - Genome enrichment

	if [ -d "$current_working_directory/Dictionaries/fake_chrom_$vcf_name" ]; then
		rm -r "$current_working_directory/Dictionaries/fake_chrom_$vcf_name"
	fi

	# START STEP 2 - genome indexing
	cd "$current_working_directory/"

	# START STEP 2.1 Reference genome indexing
	# candidate index folders, ordered by priority
	idx_folder1="${current_working_directory}/genome_library/${true_pam}_${bMax}_${ref_name}"  # index for requested number of bulges
	idx_folder2="${current_working_directory}/genome_library/${true_pam}_${bMax_}_${ref_name}"  # index for requested number of bulges + 1 (superset for required index)
	idx_folder3="${current_working_directory}/genome_library/${true_pam}_1_${ref_name}"  # index for number of bulges = 1 (used also for 0 bulges)

	# try to use an existing index
	if [ -d "$idx_folder1" ]; then
		echo "Reference Index already present"
		echo -e 'Index-genome Reference\tEnd\t'$(date) >>"$log"
		idx_ref="$idx_folder1"
	elif [ -d "$idx_folder2" ]; then
		echo "Reference Index already present"
		echo -e 'Index-genome Reference\tEnd\t'$(date) >>"$log"
		idx_ref="$idx_folder2"
	elif [ $bMax -le 1 ] && [ -d "$idx_folder3" ]; then
		echo "Reference Index already present"
		echo -e 'Index-genome Reference\tEnd\t'$(date) >>"$log"
		idx_ref="$idx_folder3"
	else
		# no valid index found, compute it
		echo -e 'Index-genome Reference\tStart\t'$(date) >>$log
		echo -e "Indexing reference genome"
		# index reference genome using crispritz
		crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus
		if [ -s $logerror ]; then
			printf "ERROR: reference genome indexing failed\n" >&2
			[ -d "$idx_folder1" ] && rm -r "$idx_folder1"
			exit 1
		fi
		pid_index_ref=$!
		echo -e 'Index-genome Reference\tEnd\t'$(date) >>$log
		idx_ref="$idx_folder1"
	fi
	# END STEP 2.1 Reference genome indexing

	# START STEP 2.2 Variant genome indexing
	if [ "$vcf_name" != "_" ]; then  # index alternative genomes (snps only)
		# candidate index folders, ordered by priority
		basedir="${current_working_directory}/genome_library"
		idx_folder1="${basedir}/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
		idx_folder2="${basedir}/${true_pam}_${bMax_}_${ref_name}+${vcf_name}"
		idx_folder3="${basedir}/${true_pam}_1_${ref_name}+${vcf_name}"

		# try existing indexes in priority order
		if [ -d "$idx_folder1" ]; then
			echo "Variant Index already present"
			echo -e 'Index-genome Variant\tEnd\t'$(date) >>"$log"
			idx_var="$idx_folder1"
		elif [ -d "$idx_folder2" ]; then
			echo "Variant Index already present"
			echo -e 'Index-genome Variant\tEnd\t'$(date) >>"$log"
			idx_var="$idx_folder2"
		elif [ "$bMax" -le 1 ] && [ -d "$idx_folder3" ]; then
			echo "Variant Index already present"
			echo -e 'Index-genome Variant\tEnd\t'$(date) >>"$log"
			idx_var="$idx_folder3"
		else
			# no index found, compute it
			echo -e 'Index-genome Variant\tStart\t'$(date) >>$log
			echo -e "Indexing variant genome"
			# index alternative genome using crispritz
			crispritz.py index-genome \
				"${ref_name}+${vcf_name}" \
				"$current_working_directory/Genomes/${ref_name}+${vcf_name}/" \
				"$pam_file" \
				-bMax $bMax -th $ncpus
			if [ -s "$logerror" ]; then
				printf "ERROR: alternative genome indexing failed on %s\n" "$vcf_name" >&2
				[ -d "$idx_folder1" ] && rm -r "$idx_folder1"
				exit 1
			fi
			pid_index_var=$!
			echo -e 'Index-genome Variant\tEnd\t'$(date) >>"$log"
			idx_var="$idx_folder1"
		fi

		# START STEP 2.3 - indels indexing
		indels_index_dir="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}_INDELS"
		indels_out="$current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS"

		if ! [ -d "$indels_index_dir" ]; then
			echo -e 'Indexing Indels\tStart\t'$(date) >>"$log"
			"$starting_dir/pool_index_indels.py" \
				"$indels_out/" \
				"$pam_file" \
				"$true_pam" \
				"$ref_name" \
				"$vcf_name" \
				"$bMax" \
				"$ncpus"

			if [ -s "$logerror" ]; then
				printf "ERROR: indels indexing failed on %s\n" "$vcf_name" >&2
				[ -d "$indels_index_dir" ] && rm -r "$indels_index_dir"
				exit 1
			fi

			echo -e 'Indexing Indels\tEnd\t'$(date) >>"$log"
		else
			echo "Indels Index already present"
		fi
		# END STEP 2.3 - indels indexing
	fi

	if [ -s $logerror ]; then
		printf "ERROR: Genome indexing failed" >&2
		exit 1
	fi
	# END STEP 2 - genome indexing

	#ceil npcus to use 1/2 of cpus per search
	ceiling_result=$((($ncpus) / 2))
	#if ceiling is 0, set ceiling to 1
	if [ $ceiling_result -eq 0 ]; then
		ceiling_result=1
	fi

	# START STEP 3 - off-targets search
	cd "$output_folder"
	pids=()  # reset process ids
	names=()
	# TODO: ricerca ref lanciata in parallelo con ricerca alternative
	echo -e 'Off-targets search\tStart\t'$(date) >>$log
	if ! [ -f "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		echo -e 'Search Reference Start'  # off-targets search on reference genome
		if [ "$bDNA" -ne 0 ] || [ "$bRNA" -ne 0 ]; then  # no bulges 
			crispritz.py search $idx_ref "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t -th $ceiling_result &
			pid_search_ref=$!
			pids+=("$pid_search_ref")  # add reference search pid
			names+=("Reference")  # add pid identifier
		else  # consider dna/rna bulges (not combined)
			crispritz.py search "$current_working_directory/Genomes/${ref_name}/" "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
			pid_search_ref=$!
			pids+=("$pid_search_ref")  # add reference search pid
			names+=("Reference")  # add pid identifier
		fi
		echo -e 'Search Reference completed'
	else
		echo -e "Search for reference already done"
	fi

	if [ "$vcf_name" != "_" ]; then
		# TODO: search in parallel on ref and alt
		if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e 'Search Variant Start'  # search off-targets on alternative genomes (snps only)
			if [ "$bDNA" -ne 0 ] || [ "$bRNA" -ne 0 ]; then  # no bulge
				crispritz.py search "$idx_var" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t -th $ceiling_result -var &
				pid_search_var=$!
				pids+=("$pid_search_var")  # add variants search pid
				names+=("Variant")  # add pid identifier
			else  # consider bulges
				crispritz.py search "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
				pid_search_var=$!
				pids+=("$pid_search_var")  # add variants search pid
				names+=("Variant")  # add pid identifier
			fi
		else
			echo -e "Search for variant already done"
		fi
		echo -e 'Search Variant completed'

		if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e "Search INDELs Start"
			cd $starting_dir
			# TODO: REMOVE POOL SCRIPT FROM PROCESSING
			./pool_search_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder" $true_pam "$current_working_directory/" "$ncpus"
			awk '($3 !~ "n") {print $0}' "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp"
			mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp" "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
		else
			echo -e "Search INDELs already done"
		fi
		echo -e "Search INDELs completed"

	fi
	
	# wait for jobs completion
	for i in "${!pids[@]}"; do
		pid="${pids[$i]}"
		name="${names[$i]}"

		if wait "$pid"; then
			if [ -s $logerror ]; then
				echo "ERROR: off-targets search ${name} failed\n" >&2
				rm -f $output_folder/*.targets.txt $output_folder/*profile*  # delete results folder
				exit 1
			fi
			echo -e "Off-targets search $name\tEnd\t"$(date) >>$log  # off-targets search on reference/variant genome
		else			
			echo "ERROR: Off-targets search $name failed" >&2
			exit 1
		fi
	done
	echo -e 'Off-targets search\tEnd\t'$(date) >>$log
	# move all targets into targets directory
	if [ -d "${output_folder}/crispritz_targets" ]; then
		mv $output_folder/*.targets.txt $output_folder/crispritz_targets &>/dev/null
	fi
	# move profiles into profile folder
	if ! [ -d "$output_folder/crispritz_prof" ]; then
		mkdir $output_folder/crispritz_prof
	fi
	if [ -d "${output_folder}/crispritz_prof" ]; then
		mv $output_folder/*profile* $output_folder/crispritz_prof/ &>/dev/null
	fi
	# END STEP 3 - off-targets search

	if [[ "$cicd_test" == "True" ]]; then  # if CI/CD test stop execution here
		exit 0
	fi

	# START STEP 4 - off-targets post-analysis
	cd "$starting_dir"
	echo -e "Start post-analysis"
	# START STEP 4.1 - off-targets post analysis snps
	if [ "$vcf_name" != "_" ]; then
		echo -e 'Post-analysis SNPs\tStart\t'$(date) >>$log
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then
			touch "$final_res_alt"
		fi
		# TODO: snp analysis in parallel
		./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $dict_folder $final_res $final_res_alt $ncpus
		if [ -s $logerror ]; then
			printf "ERROR: off-targets post-analysis (snps) failed on variants in %s\n" "$vcf_name" >&2
			rm -r $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt  # delete results folder
			exit 1
		fi
		# CONCATENATE REF&VAR RESULTS
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
		if [ -s $logerror ]; then
			printf "ERROR: off-targets post-analysis (snps) file concatenation failed on variants in %s\n" "$vcf_name" >&2
			rm -r $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt  # delete results folder
			exit 1
		fi
		echo -e 'Post-analysis SNPs\tEnd\t'$(date) >>$log
	else  # no variants -> analyze reference off-tagets
		echo -e 'Post-analysis\tStart\t'$(date) >>$log  # 
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then  # mock required to avoid crashes 
			touch "$final_res_alt"
		fi
		./pool_post_analisi_snp.py $output_folder $ref_folder "_" $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "_" $final_res $final_res_alt $ncpus
		if [ -s $logerror ]; then
			printf "ERROR: off-targets post-analysis (reference) failed\n" >&2
			rm -r $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt  # delete results folder
			exit 1
		fi
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
		if [ -s $logerror ]; then
			printf "ERROR: off-targets post-analysis (reference) file concatenation failed\n" >&2
			rm -r $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt  # delete results folder
			exit 1
		fi
		echo -e 'Post-analysis\tEnd\t'$(date) >>$log
	fi
	# END STEP 4.1 - off-targets post analysis snps

	# START STEP 4.2 - off-targets post analysis indels
	if [ "$vcf_name" != "_" ]; then
		echo -e "SNPs analysis ended. Starting INDELs analysis"
		cd "$starting_dir"
		echo -e 'Post-analysis INDELs\tStart\t'$(date) >>$log
		#SKIP INDELS ANALYSIS IF NO RESULTS FOUND
		if [ $(wc -l <"$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt") -gt 1 ]; then
			./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "$current_working_directory/Dictionaries/" $final_res $final_res_alt $ncpus
			if [ -s $logerror ]; then
				printf "ERROR: off-targets post-analysis (indels) failed on variants in %s\n" "$vcf_name" >&2
				rm -r $output_folder/*.bestCFD*.txt $output_folder/*.bestmmblg*.txt $output_folder/*.bestCRISTA*.txt  # delete results folder
				exit 1
			fi
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
			if [ -s $logerror ]; then
				printf "ERROR: off-targets post-analysis (indels) file concatenation failed on variants in %s\n" "$vcf_name" >&2
				rm -r $output_folder/*.bestCFD*.txt $output_folder/*.bestmmblg*.txt $output_folder/*.bestCRISTA*.txt  # delete results folder
				exit 1
			fi
		fi
		echo -e 'Post-analysis INDELs\tEnd\t'$(date) >>$log
	fi
	# END STEP 4.2 - off-targets post analysis indels
	# END STEP 4 - off-targets post-analysis
done <$vcf_list

echo -e "Adding header to files"
while read samples; do
	if [ -z "$samples" ]; then
		continue
	fi
	# copy samples ids file content to temporary sample ids file (skip header)
	grep -v '#' "${current_working_directory}/samplesIDs/$samples" >>"$output_folder/.sampleID.txt"
	if [ -s $logerror ]; then
		printf "ERROR: samples IDs processing failed on dataset %s\n" "$samples" >&2
		rm $output_folder/.sampleID.txt
		exit 1
	fi
done <$sampleID
touch "$output_folder/.sampleID.txt"  # if not created, create it to avoid potentil crashes in only reference searches
sed -i 1i"#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX" "$output_folder/.sampleID.txt"  # add header to sampleID file
sampleID=$output_folder/.sampleID.txt  # point to new sample ID file
if [ -s $logerror ]; then
	printf "ERROR: samples IDs files processing failed\n" >&2
	rm $sampleID
	exit 1
fi

# add header to primary results files 
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCFD.txt"
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestmmblg.txt"
sed -i '1i #Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCRISTA.txt"
if [  -s $logerror ]; then 
	printf "ERROR: failed adding headers to primary results files\n" >&2
	rm $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
# add header to alternative results files
echo "header" >$final_res_alt.bestCFD.txt
echo "header" >$final_res_alt.bestmmblg.txt
echo "header" >$final_res_alt.bestCRISTA.txt
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestCFD.txt"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestmmblg.txt"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref/' "$final_res_alt.bestCRISTA.txt"
if [  -s $logerror ]; then 
	printf "ERROR: failed adding headers to alternative results files\n" >&2
	rm $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi

# START STEP 5 - targets merge
# sort files to have chromosomes and positions in proximity -> simplify merge
echo -e 'Merging Targets\tStart\t'$(date) >>$log
# sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCFD.txt >$final_res.tmp
tail -n +2 $final_res.bestCFD.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCFD.txt
#sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCRISTA.txt >$final_res.tmp
tail -n +2 $final_res.bestCRISTA.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCRISTA.txt
#sort using guide_seq,chr,cluster_pos,total(mm+bul)
head -1 $final_res.bestmmblg.txt >$final_res.tmp
tail -n +2 $final_res.bestmmblg.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestmmblg.txt
if [ -s $logerror ]; then 
	printf "ERROR: results file preprocessing failed for targets merge step\n" >&2
	rm $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
# merge contiguous targets (proximity threshold defined by merge_t)
# CFD score
./merge_close_targets_cfd.sh $final_res.bestCFD.txt $final_res.bestCFD.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
# MM+BUL counts
./merge_close_targets_cfd.sh $final_res.bestmmblg.txt $final_res.bestmmblg.txt.trimmed $merge_t 'total' $sorting_criteria_scoring $sorting_criteria &
# CRISTA score
./merge_close_targets_cfd.sh $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
wait
if [ -s $logerror ]; then
	printf "ERROR: merging contiguous targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt* $output_folder/*.bestmmblg.txt* $output_folder/*.bestCRISTA.txt* $output_folder/*.bestMerge.txt* $output_folder/*.altMerge.txt*
	exit 1
fi
# rename primary and alternative results files
mv $final_res.bestCFD.txt.trimmed $final_res.bestCFD.txt
mv $final_res.bestCFD.txt.trimmed.discarded_samples $final_res_alt.bestCFD.txt
mv $final_res.bestmmblg.txt.trimmed $final_res.bestmmblg.txt
mv $final_res.bestmmblg.txt.trimmed.discarded_samples $final_res_alt.bestmmblg.txt
mv $final_res.bestCRISTA.txt.trimmed $final_res.bestCRISTA.txt
mv $final_res.bestCRISTA.txt.trimmed.discarded_samples $final_res_alt.bestCRISTA.txt
echo -e 'Merging Targets\tEnd\t'$(date) >>$log
# END STEP 5 - targets merge

# START STEP 6 - targets annotation
echo -e 'Annotating results\tStart\t'$(date) >>$log
# annotate primary targets 
python annotation.py $final_res.bestCFD.txt $annotation_file $final_res.bestCFD.txt.annotated &
python annotation.py $final_res.bestmmblg.txt $annotation_file $final_res.bestmmblg.txt.annotated &
python annotation.py $final_res.bestCRISTA.txt $annotation_file $final_res.bestCRISTA.txt.annotated &
wait
mv $final_res.bestCFD.txt.annotated $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.annotated $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.annotated $final_res.bestCRISTA.txt
if [ -s $logerror ]; then
	printf "ERROR: primary targets annotation failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt* $output_folder/*.bestmmblg.txt* $output_folder/*.bestCRISTA.txt* $output_folder/*.bestMerge.txt* $output_folder/*.altMerge.txt*
	exit 1
fi
# annotate alternative targets
python annotation.py $final_res_alt.bestCFD.txt $annotation_file $final_res_alt.bestCFD.txt.annotated &
python annotation.py $final_res_alt.bestmmblg.txt $annotation_file $final_res_alt.bestmmblg.txt.annotated &
python annotation.py $final_res_alt.bestCRISTA.txt $annotation_file $final_res_alt.bestCRISTA.txt.annotated &
wait
mv $final_res_alt.bestCFD.txt.annotated $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.annotated $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.annotated $final_res_alt.bestCRISTA.txt
if [ -s $logerror ]; then
	printf "ERROR: alternative targets annotation failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt* $output_folder/*.bestmmblg.txt* $output_folder/*.bestCRISTA.txt* $output_folder/*.bestMerge.txt* $output_folder/*.altMerge.txt*
	exit 1
fi
# compute risk scores for primary targets
./add_risk_score.py $final_res.bestCFD.txt $final_res.bestCFD.txt.risk "False" &
./add_risk_score.py $final_res.bestmmblg.txt $final_res.bestmmblg.txt.risk "False" &
./add_risk_score.py $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.risk "False" &
wait
mv $final_res.bestCFD.txt.risk $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.risk $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.risk $final_res.bestCRISTA.txt
if [ -s $logerror ]; then
	printf "ERROR: computing risk scores on primary targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt* $output_folder/*.bestmmblg.txt* $output_folder/*.bestCRISTA.txt* $output_folder/*.bestMerge.txt* $output_folder/*.altMerge.txt*
	exit 1
fi
# compute risk scores for alternative targets
./add_risk_score.py $final_res_alt.bestCFD.txt $final_res_alt.bestCFD.txt.risk "False" &
./add_risk_score.py $final_res_alt.bestmmblg.txt $final_res_alt.bestmmblg.txt.risk "False" &
./add_risk_score.py $final_res_alt.bestCRISTA.txt $final_res_alt.bestCRISTA.txt.risk "False" &
wait
mv $final_res_alt.bestCFD.txt.risk $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.risk $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.risk $final_res_alt.bestCRISTA.txt
if [ -s $logerror ]; then
	printf "ERROR: computing risk scores on alternative targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt* $output_folder/*.bestmmblg.txt* $output_folder/*.bestCRISTA.txt* $output_folder/*.bestMerge.txt* $output_folder/*.altMerge.txt*
	exit 1
fi
# remove Ns and dots from rsID from primary targets files
./remove_n_and_dots.py $final_res.bestCFD.txt &
./remove_n_and_dots.py $final_res.bestmmblg.txt &
./remove_n_and_dots.py $final_res.bestCRISTA.txt &
wait
if [ -s $logerror ]; then
	printf "ERROR: rsids NaN values replacement on primary targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
# remove Ns and dots from rsID from alternative targets files
./remove_n_and_dots.py $final_res_alt.bestCFD.txt &
./remove_n_and_dots.py $final_res_alt.bestmmblg.txt &
./remove_n_and_dots.py $final_res_alt.bestCRISTA.txt &
wait
if [ -s $logerror ]; then
	printf "ERROR: rsids NaN values replacement on alternative targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
# join targets by columns for primary and alternative results
pr -m -t -J $final_res.bestCFD.txt $final_res.bestmmblg.txt $final_res.bestCRISTA.txt >$final_res
if [ -s $logerror ]; then
	printf "ERROR: targets join on primary targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
pr -m -t -J $final_res_alt.bestCFD.txt $final_res_alt.bestmmblg.txt $final_res_alt.bestCRISTA.txt >$final_res_alt
if [ -s $logerror ]; then
	printf "ERROR: targets join on alternative targets failed\n" >&2
	rm -f $output_folder/*.bestCFD.txt $output_folder/*.bestmmblg.txt $output_folder/*.bestCRISTA.txt $output_folder/*.bestMerge.txt $output_folder/*.altMerge.txt
	exit 1
fi
# update headers for primary and alternative results 
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res"
if [ -s $logerror ]; then
	printf "ERROR: header update on primary targets failed\n" >&2
	rm $final_res* $final_res_alt*
	exit 1
fi	
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res_alt"
if [ -s $logerror ]; then
	printf "ERROR: header update on alternative targets failed\n" >&2
	rm $final_res* $final_res_alt*
	exit 1
fi	
echo -e 'Annotating results\tEnd\t'$(date) >>$log
# END STEP 6 - targets annotation

# START STEP 7 - graphical reports
echo -e 'Creating images\tStart\t'$(date) >>$log
cd $output_folder
# adjust filenames and remove redundant files
echo -e "Cleaning directory"
rm -f *.CFDGraph.txt
rm -f indels.CFDGraph.txt
rm -r "crispritz_prof"
rm -r "crispritz_targets" # remove targets in online version to avoid memory saturation
# change primary and alt filenames
mv $final_res "${output_folder}/$(basename ${output_folder}).bestMerge.txt"
mv $final_res_alt "${output_folder}/$(basename ${output_folder}).altMerge.txt"
# create result summaries for primary and alternative results
cd $starting_dir
if [ "$vcf_name" != "_" ]; then  # variants available
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "var" "CFD"
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "var" "fewest"
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "var" "CRISTA"
	if [ -s $logerror ]; then
		printf "ERROR: summary processing failed (variants pipeline)\n" >&2
		rm -f $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt
		exit 1
	fi	
else  # only reference search
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "ref" "CFD"
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "ref" "fewest"
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bDNA $bRNA "${output_folder}" "ref" "CRISTA"
	if [ -s $logerror ]; then
		printf  "ERROR: summary processing failed (reference genome pipeline)\n" >&2
		rm -f $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt
		exit 1
	fi
fi
if ! [ -d "$output_folder/imgs" ]; then  # create images folder
	mkdir "$output_folder/imgs"
fi
if [ "$vcf_name" != "_" ]; then  # variants available -> create population distribution plots
	cd "$output_folder/imgs"
	while IFS= read -r line || [ -n "$line" ]; do
		for total in $(seq 0 $(expr $mm + $bMax)); do
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CFD.txt" $total $line "CFD"
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CRISTA.txt" $total $line "CRISTA"
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_fewest.txt" $total $line "fewest"
			if [ -s $logerror ]; then
				printf "ERROR: population distribution plots creation failed for guide %s (mm+bulges: %d)\n" "$line" "$total" >&2
				rm -r "${output_folder}/imgs"
				exit 1
			fi
		done
	done <$guide_file
fi
cd $starting_dir
# generate radar charts
if [ "$vcf_name" != "_" ]; then
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CFD"
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA"
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "fewest"
	if [ -s $logerror ]; then
		printf  "ERROR: summary processing failed (variants pipeline)\n" >&2
		rm -r "${output_folder}/imgs"
		rm -f $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt
		exit 1
	fi
else
	echo -e "dummy_file" >dummy.txt
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt dummy.txt $annotation_file "$output_folder" $ncpus $mm $bMax "CFD"
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt dummy.txt $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA"
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt dummy.txt $annotation_file "$output_folder" $ncpus $mm $bMax "fewest"
	rm dummy.txt
	if [ -s $logerror ]; then
		printf  "ERROR: summary processing failed (reference genome pipeline)\n" >&2
		rm -r "${output_folder}/imgs"
		rm -f $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt
		exit 1
	fi
fi
echo -e 'Creating images\tEnd\t'$(date) >>$log
# END STEP 7 - graphical reports

# START STEP 8 - results integration
echo -e 'Integrating results\tStart\t'$(date) >>$log
echo >>$guide_file
if [ $gene_proximity != "_" ]; then
	touch "${output_folder}/dummy.txt"
	genome_version=$(echo ${ref_name} | sed 's/_ref//' | sed -e 's/\n//') #${output_folder}/Params.txt | awk '{print $2}' | sed 's/_ref//' | sed -e 's/\n//')
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${gene_proximity}" "${output_folder}/dummy.txt" "${guide_file}" $genome_version "${output_folder}" "vuota" $starting_dir/ $base_check_start $base_check_end $base_check_set
	if [ -s $logerror ]; then
		printf  "ERROR: targets integration failed on primary results\n" >&2
		rm $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt $output_folder/*.tsv
		exit 1
	fi
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).altMerge.txt" "${gene_proximity}" "${output_folder}/dummy.txt" "${guide_file}" $genome_version "${output_folder}" "vuota" $starting_dir/ $base_check_start $base_check_end $base_check_set
	if [ -s $logerror ]; then
		printf  "ERROR: targets integration failed on primary results\n" >&2
		rm $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt $output_folder/*.tsv
		exit 1
	fi
	rm "${output_folder}/dummy.txt"
	python $starting_dir/CRISPRme_plots.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/imgs/" &>"${output_folder}/warnings.txt"
	if [ -s $logerror ]; then
		printf  "ERROR: score plots generation failed\n" >&2
		rm -r "${output_folder}/imgs"
		rm $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt $output_folder/*.tsv
		exit 1
	fi
	rm -f "${output_folder}/warnings.txt" # delete warnings file
fi
echo -e 'Integrating results\tEnd\t'$(date) >>$log
truncate -s -1 $guide_file
truncate -s -1 $vcf_list
if [ $6 != "_" ]; then
	truncate -s -1 $6
fi
# END STEP 8 - results integration

# START STEP 9 - database creation
echo -e 'Building database'
echo -e 'Creating database\tStart\t'$(date) >>$log
if [ -f "${output_folder}/$(basename ${output_folder}).db" ]; then
	rm -f "${output_folder}/$(basename ${output_folder}).db"
fi
python $starting_dir/db_creation.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/.$(basename ${output_folder})"
if [ -s $logerror ]; then
	printf "ERROR: database creation failed\n" >&2 
	rm $final_res* $final_res_alt* $output_folder/*.altMerge.txt $output_folder/*.bestMerge.txt $output_folder/*_CFD.txt $output_folder/*_fewest.txt $output_folder/*_CRISTA.txt $output_folder/.*_CFD.txt $output_folder/.*_fewest.txt $output_folder/.*_CRISTA.txt $output_folder/*.tsv
	exit 1
fi
echo -e 'Creating database\tEnd\t'$(date) >>$log
# END STEP 9 - database creation

echo -e 'Job\tDone\t'$(date) >>$log

if [ $(wc -l <"$guide_file") -gt 1 ]; then
	mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	mv "${output_folder}/$(basename ${output_folder}).altMerge.txt.integrated_results.tsv" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	# zip primary and alternative results
	zip -j "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.zip" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	zip -j "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
else
	guide_elem=$(head -1 $guide_file)
	mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	mv "${output_folder}/$(basename ${output_folder}).altMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	# zip primary and alternative results
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
fi
echo -e "JOB END"

if [ "$email" != "_" ]; then
	python $starting_dir/../pages/send_mail.py $output_folder
fi

# keep log_error but no block visualization
mv $output_folder/log_error.txt $output_folder/log_error_no_check.txt
# removing single best files after use and clean merged file to save space
# keep the two targets files
rm $final_res.bestCFD.txt
rm $final_res.bestmmblg.txt
rm $final_res.bestCRISTA.txt
rm $final_res_alt.bestCFD.txt
rm $final_res_alt.bestmmblg.txt
rm $final_res_alt.bestCRISTA.txt
