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

#file for automated search of guide+pam in reference and variant genomes

ref_folder=$(realpath $1)  # reference genome folder
vcf_list=$(realpath $2)  # vcf folders list
guide_file=$(realpath $3)  # guide 
pam_file=$(realpath $4)  # pam
annotation_file=$(realpath $5)  # annotation bed
sampleID=$(realpath $6)  # sample ids 
bMax=$7  # max number of bulges
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

#used to solve base editor check in resultintegration phase
base_check_start=${18}
base_check_end=${19}
base_check_set=${20}

# sorting criteria while merging best targets
sorting_criteria_scoring=${21}
sorting_criteria=${22}

# create log files
log="${output_folder}/log.txt"
touch $log  
logerror="${output_folder}/log_error.txt"  # log error -> trace errors
#echo -e 'Job\tStart\t'$(date) > $log
start_time='Job\tStart\t'$(date)

# output=$output_folder/output.txt
# touch $output
## CREATE DUMMY FILE WITH ONE LINE
echo -e "dummy_file" >"${output_folder}/.dummy.txt"
dummy_file="${output_folder}/.dummy.txt"
## CREATE EMPTY FILE
touch "${output_folder}/.empty.txt"
empty_file="${output_folder}/.empty.txt"
## CREATE EMPTY DIR
mkdir -p "${output_folder}/.empty"
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

# perform target search and processing for each variants dataset
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

	# STEP 1: Enrich genome adding SNPs and INDELs from the input VCF files
	# track haplotypes if input VCF is phased
	if [ "$vcf_name" != "_" ]; then

		cd "$current_working_directory/Genomes"
		if ! [ -d "$current_working_directory/Genomes/${ref_name}+${vcf_name}" ]; then
			echo -e 'Add-variants\tStart\t'$(date) >>$log
			# echo -e 'Add-variants\tStart\t'$(date) >&2
			echo -e "Adding variants"
			crispritz.py add-variants "$vcf_folder/" "$ref_folder/" "true" 
			# check for add-variants failures
			if [ -s $logerror ]; then 
				printf "ERROR: Genome enrichment failed!\n" >&2
				exit 1
			fi
			mv "$current_working_directory/Genomes/variants_genome/SNPs_genome/${ref_name}_enriched/" "./${ref_name}+${vcf_name}/"
			if ! [ -d "$current_working_directory/Dictionaries/dictionaries_${vcf_name}/" ]; then
				mkdir "$current_working_directory/Dictionaries/dictionaries_${vcf_name}/"
			fi
			# check for snp dictionary failures
			if [ -s $logerror ]; then 
				printf "ERROR: SNP dictionary construction failed!\n" >&2
				exit 1
			fi
			if ! [ -d "$current_working_directory/Dictionaries/log_indels_${vcf_name}/" ]; then
				mkdir "$current_working_directory/Dictionaries/log_indels_${vcf_name}/"
			fi
			# check for indel dictionary failures
			if [ -s $logerror ]; then 
				printf "ERROR: Indel dictionary construction failed!\n" >&2
				exit 1
			fi
			mv $current_working_directory/Genomes/variants_genome/SNPs_genome/*.json $current_working_directory/Dictionaries/dictionaries_${vcf_name}/
			mv $current_working_directory/Genomes/variants_genome/SNPs_genome/log*.txt $current_working_directory/Dictionaries/log_indels_${vcf_name}/
			cd "$current_working_directory/"
			if ! [ -d "genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS" ]; then
				mkdir "genome_library/${true_pam}_2_${ref_name}+${vcf_name}_INDELS"
			fi
			# check for genome library failures
			if [ -s $logerror ]; then 
				printf "ERROR: Genome library construction failed!\n" >&2
				exit 1
			fi

			echo -e 'Add-variants\tEnd\t'$(date) >>$log
			# echo -e 'Add-variants\tEnd\t'$(date) >&2

			# STEP 2: indels indexing 
			echo -e 'Indexing Indels\tStart\t'$(date) >>$log
			# echo -e 'Indexing Indels\tStart\t'$(date) >&2
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/variants_genome/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus 
			# check for indels indexing failures
			if [ -s $logerror ]; then 
				printf "ERROR: Indels indexing failed!\n" >&2
				exit 1
			fi
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
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus 
			echo -e 'Indexing Indels\tEnd\t'$(date) >>$log
			# echo -e 'Indexing Indels\tEnd\t'$(date) >&2
			# check for indels indexing failures
			if [ -s $logerror ]; then 
				printf "ERROR: Indels indexing failed!\n" >&2
				exit 1
			fi
		fi
	fi

	if [ -d "$current_working_directory/Dictionaries/fake_chrom_$vcf_name" ]; then
		rm -r "$current_working_directory/Dictionaries/fake_chrom_$vcf_name"
	fi

	# STEP 3: index reference genome
	cd "$current_working_directory/"
	if ! [ -d "$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}" ]; then
			if ! [ $bMax -gt 1 ]; then
				if ! [ -d "$current_working_directory/genome_library/${true_pam}_1_${ref_name}" ]; then
					echo -e 'Index-genome Reference\tStart\t'$(date) >>$log
					# echo -e 'Index-genome Reference\tStart\t'$(date) >&2
					# echo -e 'Indexing_Reference' > $output
					echo -e "Indexing reference genome"
					crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus
					pid_index_ref=$!
					# check for reference genome indexing failures
					if [ -s $logerror ]; then 
						printf "ERROR: Reference genome indexing failed!\n" >&2
						exit 1
					fi
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
				crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus 
				pid_index_ref=$!
				# check for reference genome indexing failures
				if [ -s $logerror ]; then 
					printf "ERROR: Reference genome indexing failed!\n" >&2
					exit 1
				fi
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

	# STEP 4: index variant genome
	if [ "$vcf_name" != "_" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}" ]; then
			if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}" ]; then
				if ! [ $bMax -gt 1 ]; then
					if ! [ -d "$current_working_directory/genome_library/${true_pam}_1_${ref_name}+${vcf_name}" ]; then
						echo -e 'Index-genome Variant\tStart\t'$(date) >>$log
						# echo -e 'Index-genome Variant\tStart\t'$(date) >&2
						# echo -e 'Indexing_Enriched' > $output
						echo -e "Indexing variant genome"
						crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus 
						pid_index_var=$!
						# check for variant genome indexing failures
						if [ -s $logerror ]; then 
							printf "ERROR: Variant genome indexing failed!\n" >&2
							exit 1
						fi
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
					crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus 
					pid_index_ref=$!
					# check for variant genome indexing failures
					if [ -s $logerror ]; then 
						printf "ERROR: Variant genome indexing failed!\n" >&2
						exit 1
					fi
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

	# STEP 5: reference genome search
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
			pid_search_ref=$!
			# check for reference genome search (brute-force) failures
			if [ -s $logerror ]; then 
				printf "ERROR: Reference genome search (brute-force) failed!\n" >&2
				exit 1
			fi
			echo -e 'Search Reference\tEnd\t'$(date) >>$log
		else
			crispritz.py search "$current_working_directory/Genomes/${ref_name}/" "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
			pid_search_ref=$!
			# check for reference genome search (TST) failures
			if [ -s $logerror ]; then 
				printf "ERROR: Reference genome search (TST) failed!\n" >&2
				exit 1
			fi
			echo -e 'Search Reference\tEnd\t'$(date) >>$log
		fi
	else
		echo -e "Search for reference already done"
	fi

	# STEP 6: variant genome search
	if [ "$vcf_name" != "_" ]; then
		if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e 'Search Variant\tStart\t'$(date) >>$log
			# echo -e 'Search Variant\tStart\t'$(date) >&2
			# echo -e 'Search Variant' >  $output
			if [ "$bDNA" -ne 0 ] || [ "$bRNA" -ne 0 ]; then
				crispritz.py search "$idx_var" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t -th $ceiling_result -var 
				# mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
				echo -e 'Search Variant\tEnd\t'$(date) >>$log
				# check for variant genome search (brute-force) failures
				if [ -s $logerror ]; then 
					printf "ERROR: Variant genome search (brute-force) failed!\n" >&2
					exit 1
				fi
			else
				crispritz.py search "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -mm $mm -r -th $ceiling_result &
				# check for variant genome search (TST) failures
				if [ -s $logerror ]; then 
					printf "ERROR: Variant genome search (TST) failed!\n" >&2
					exit 1
				fi	
				echo -e 'Search Variant\tEnd\t'$(date) >>$log
				# mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			fi
		else
			echo -e "Search for variant already done"
		fi

		# STEP 7: search on indels
		if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e "Search INDELs Start"
			echo -e 'Search INDELs\tStart\t'$(date) >>$log
			# echo -e 'Search INDELs\tStart\t'$(date) >&2
			cd $starting_dir
			#commented to avoid indels search
			#TODO REMOVE POOL SCRIPT FROM PROCESSING
			./pool_search_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder" $true_pam "$current_working_directory/" "$ncpus" 
			# check for indels genome search failures
			if [ -s $logerror ]; then 
				printf "ERROR: Indels genome search failed!\n" >&2
				exit 1
			fi
			# mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			awk '($3 !~ "n") {print $0}' "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp"
			mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.tmp" "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
			echo -e "Search INDELs End"
			echo -e 'Search INDELs\tEnd\t'$(date) >>$log
			# echo -e 'Search INDELs\tEnd\t'$(date) >&2
		else
			echo -e "Search INDELs already done"
		fi
	fi

	while kill "-0" $pid_search_ref &>/dev/null; do
		echo -e "Waiting for search genome"
		sleep 100
	done
	echo -e 'Search\tEnd\t'$(date) >>$log
	# echo -e 'Search Reference\tEnd\t'$(date) >&2

	# move all targets into targets directory
	mv $output_folder/*.targets.txt $output_folder/crispritz_targets
	# check for targets folder creation failures
	if [ -s $logerror ]; then 
		printf "ERROR: Targets folder creation failed!\n" >&2
		exit 1
	fi

	if ! [ -d "$output_folder/crispritz_prof" ]; then
		mkdir $output_folder/crispritz_prof
	fi
	mv $output_folder/*profile* $output_folder/crispritz_prof/ &>/dev/null
	# check for profile folder creation failures
	if [ -s $logerror ]; then 
		printf "ERROR: Profile folder creation failed!\n" >&2
		exit 1
	fi

	cd "$starting_dir"

	# STEP 8: snp analysis
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
		./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $dict_folder $final_res $final_res_alt $ncpus 
		# check for snp analysis failures
		if [ -s $logerror ]; then 
			printf "ERROR: SNP analysis failed!\n" >&2
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
		# check for reports creation failures
		if [ -s $logerror ]; then 
			printf "ERROR: Temporary reports (snp) creation failed!\n" >&2
			exit 1
		fi

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

		./pool_post_analisi_snp.py $output_folder $ref_folder "_" $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "_" $final_res $final_res_alt $ncpus
		# check for targets analysis failures
		if [ -s $logerror ]; then 
			printf "ERROR: Targets analysis failed!\n" >&2
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
		# check for reports creation failures
		if [ -s $logerror ]; then 
			printf "ERROR: Temporary reports (reference) creation failed!\n" >&2
			exit 1
		fi
		echo -e 'Post-analysis\tEnd\t'$(date) >>$log
	fi

	# STEP 9: indels analysis
	if [ "$vcf_name" != "_" ]; then
		echo -e "SNPs analysis ended. Starting INDELs analysis"
		cd "$starting_dir"

		echo -e 'Post-analysis INDELs\tStart\t'$(date) >>$log
		#SKIP INDELS ANALYSIS IF NO RESULTS FOUND
		if [ $(wc -l <"$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt") -gt 1 ]; then

			./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "$current_working_directory/Dictionaries/" $final_res $final_res_alt $ncpus
			# check for indels analysis failures
			if [ -s $logerror ]; then 
				printf "ERROR: Indels analysis failed!\n" >&2
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
			# check for reports creation failures
			if [ -s $logerror ]; then 
				printf "ERROR: Temporary reports (indels) creation failed!\n" >&2
				exit 1
			fi
		fi
		echo -e 'Post-analysis INDELs\tEnd\t'$(date) >>$log

	fi
done <$vcf_list

echo -e "Adding header to files"

while read samples; do
	if [ -z "$samples" ]; then
		continue
	fi
	awk '!/^#/ { print }' "${current_working_directory}/samplesIDs/$samples" >>"$output_folder/.sampleID.txt"
done <"$sampleID"
# check for samples ids reading failures
if [ -s $logerror ]; then 
	printf "ERROR: Reading sample IDs failed!\n" >&2
	exit 1
fi
# done <$sampleID
# if [ "$vcf_name" != "_" ]; then
touch "$output_folder/.sampleID.txt"
sed -i 1i"#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX" "$output_folder/.sampleID.txt"
# fi

sampleID=$output_folder/.sampleID.txt

# echo -e 'Merging targets' >  $output

#create result file for each scoring method
header="#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref"
#header into final_res best
sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCFD.txt"
sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestmmblg.txt"
sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref' "$final_res.bestCRISTA.txt"
#header into final_res alt
printf $header >$final_res_alt.bestCFD.txt
printf $header >$final_res_alt.bestmmblg.txt
printf $header >$final_res_alt.bestCRISTA.txt
# check for reports creation failures
if [ -s $logerror ]; then 
	printf "ERROR: Report files creation failed!\n" >&2
	exit 1
fi

# STEP 10: Merging contiguous targets
echo -e 'Merging Targets\tStart\t'$(date) >>$log
#SORT FILE TO HAVE CHR AND POS IN PROXIMITY TO MERGE THEM
#sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCFD.txt >$final_res.tmp
tail -n +2 $final_res.bestCFD.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCFD.txt
# check for CFD report sorting failures
if [ -s $logerror ]; then 
	printf "ERROR: Sorting CFD report failed!\n" >&2
	exit 1
fi
#sort using guide_seq,chr,cluster_pos,score,total(mm+bul)
head -1 $final_res.bestCRISTA.txt >$final_res.tmp
tail -n +2 $final_res.bestCRISTA.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k21,21rg -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestCRISTA.txt
# check for CRISTA report sorting failures
if [ -s $logerror ]; then 
	printf "ERROR: Sorting CRISTA report failed!\n" >&2
	exit 1
fi
#sort using guide_seq,chr,cluster_pos,total(mm+bul)
head -1 $final_res.bestmmblg.txt >$final_res.tmp
tail -n +2 $final_res.bestmmblg.txt | LC_ALL=C sort -k16,16 -k5,5 -k7,7n -k11,11n -T $output_folder >>$final_res.tmp && mv $final_res.tmp $final_res.bestmmblg.txt
# check for mm+bulges report sorting failures
if [ -s $logerror ]; then 
	printf "ERROR: Sorting mm+bulges report failed!\n" >&2
	exit 1
fi

# cp $final_res.bestCFD.txt $final_res.sorted.bestCFD.txt
#MERGE BEST FILES TARGETS TO REMOVE CONTIGOUS
#TODO CHECK MERGE
#SCORE CFD
./merge_close_targets_cfd.sh $final_res.bestCFD.txt $final_res.bestCFD.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
# check for targets merge on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: merging targets in CFD report failed!\n" >&2
	exit 1
fi
#TOTAL (MM+BUL)
./merge_close_targets_cfd.sh $final_res.bestmmblg.txt $final_res.bestmmblg.txt.trimmed $merge_t 'total' $sorting_criteria_scoring $sorting_criteria &
# check for targets merge on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: merging targets in mm+bulges report failed!\n" >&2
	exit 1
fi
#SCORE CRISTA
./merge_close_targets_cfd.sh $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.trimmed $merge_t 'score' $sorting_criteria_scoring $sorting_criteria &
# check for targets merge on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: merging targets in CRISTA report failed!\n" >&2
	exit 1
fi
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

# STEP 11: targets annotation 
#ANNOTATE BEST TARGETS
#TODO SISTEMARE ANNOTAZIONE (DIVISIONE INTERVAL TREE / PARALLEL SEARCH)
./annotate_final_results.py $final_res.bestCFD.txt $annotation_file $final_res.bestCFD.txt.annotated &
# check for annotation on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in CFD report failed!\n" >&2
	exit 1
fi
./annotate_final_results.py $final_res.bestmmblg.txt $annotation_file $final_res.bestmmblg.txt.annotated &
# check for annotation on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in mm+bulges report failed!\n" >&2
	exit 1
fi
./annotate_final_results.py $final_res.bestCRISTA.txt $annotation_file $final_res.bestCRISTA.txt.annotated &
# check for annotation on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in CRISTA report failed!\n" >&2
	exit 1
fi
wait
mv $final_res.bestCFD.txt.annotated $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.annotated $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.annotated $final_res.bestCRISTA.txt
#ANNOTATE ALT TARGETS
./annotate_final_results.py $final_res_alt.bestCFD.txt $annotation_file $final_res_alt.bestCFD.txt.annotated &
# check for annotation on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in CFD alternative report failed!\n" >&2
	exit 1
fi
./annotate_final_results.py $final_res_alt.bestmmblg.txt $annotation_file $final_res_alt.bestmmblg.txt.annotated &
# check for annotation on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in mm+bulges alternative report failed!\n" >&2
	exit 1
fi
./annotate_final_results.py $final_res_alt.bestCRISTA.txt $annotation_file $final_res_alt.bestCRISTA.txt.annotated &
# check for annotation on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: Targets annotation in CRISTA alternative report failed!\n" >&2
	exit 1
fi
wait
mv $final_res_alt.bestCFD.txt.annotated $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.annotated $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.annotated $final_res_alt.bestCRISTA.txt

# STEP 12: compute risk scores 
#SCORING BEST RESULTS
./add_risk_score.py $final_res.bestCFD.txt $final_res.bestCFD.txt.risk "False" &
# check for risk score computing on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in CFD report failed!\n" >&2
	exit 1
fi
./add_risk_score.py $final_res.bestmmblg.txt $final_res.bestmmblg.txt.risk "False" &
# check for risk score computing on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in mm+bulges report failed!\n" >&2
	exit 1
fi
./add_risk_score.py $final_res.bestCRISTA.txt $final_res.bestCRISTA.txt.risk "False" &
# check for risk score computing on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in CRISTA report failed!\n" >&2
	exit 1
fi
wait
mv $final_res.bestCFD.txt.risk $final_res.bestCFD.txt
mv $final_res.bestmmblg.txt.risk $final_res.bestmmblg.txt
mv $final_res.bestCRISTA.txt.risk $final_res.bestCRISTA.txt
#SCORING ALT RESULTS
./add_risk_score.py $final_res_alt.bestCFD.txt $final_res_alt.bestCFD.txt.risk "False" &
# check for risk score computing on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in CFD alternative report failed!\n" >&2
	exit 1
fi
./add_risk_score.py $final_res_alt.bestmmblg.txt $final_res_alt.bestmmblg.txt.risk "False" &
# check for risk score computing on mm_bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in mm+bulges alternative report failed!\n" >&2
	exit 1
fi
./add_risk_score.py $final_res_alt.bestCRISTA.txt $final_res_alt.bestCRISTA.txt.risk "False" &
# check for risk score computing on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: Risk score in CRISTA alternative report failed!\n" >&2
	exit 1
fi
wait
mv $final_res_alt.bestCFD.txt.risk $final_res_alt.bestCFD.txt
mv $final_res_alt.bestmmblg.txt.risk $final_res_alt.bestmmblg.txt
mv $final_res_alt.bestCRISTA.txt.risk $final_res_alt.bestCRISTA.txt

# STEP 13: clean reports from dots and NaN values
#remove N's and dots from rsID from BEST FILES
python remove_n_and_dots.py $final_res.bestCFD.txt &
# check for NaN values cleaning on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in CFD report failed!\n" >&2
	exit 1
fi
python remove_n_and_dots.py $final_res.bestmmblg.txt &
# check for NaN values cleaning on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in mm+bulges report failed!\n" >&2
	exit 1
fi
python remove_n_and_dots.py $final_res.bestCRISTA.txt &
# check for NaN values cleaning on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in CRISTA report failed!\n" >&2
	exit 1
fi
wait
#remove N's and dots from rsID from ALT FILES
python remove_n_and_dots.py $final_res_alt.bestCFD.txt &
# check for NaN values cleaning on CFD failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in CFD alternative report failed!\n" >&2
	exit 1
fi
python remove_n_and_dots.py $final_res_alt.bestmmblg.txt &
# check for NaN values cleaning on mm+bulges failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in mm+bulges alternative report failed!\n" >&2
	exit 1
fi
python remove_n_and_dots.py $final_res_alt.bestCRISTA.txt &
# check for NaN values cleaning on CRISTA failures
if [ -s $logerror ]; then 
	printf "ERROR: NaN values cleaning in CRISTA alternative report failed!\n" >&2
	exit 1
fi
wait

#join targets by columns for BEST and ALT files
pr -m -t -J $final_res.bestCFD.txt $final_res.bestmmblg.txt $final_res.bestCRISTA.txt >$final_res 
pr -m -t -J $final_res_alt.bestCFD.txt $final_res_alt.bestmmblg.txt $final_res_alt.bestCRISTA.txt >$final_res_alt 

#MERGE ALTERNATIVE CHR IF SAME SEQUENCE OF ALIGNED CHR
# ./merge_alt_chr.sh $final_res $final_res.chr_merged
# mv $final_res.chr_merged $final_res

#update header for final_res and final_res_alt
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref\tCRISTA_CFD_Risk_Score\tCRISTA_CFD_Absolute_Risk_Score/' "$final_res_alt"
# check for report headers update failures
if [ -s $logerror ]; then 
	printf "ERROR: Updating report headers failed!\n" >&2
	exit 1
fi

echo -e 'Annotating results\tEnd\t'$(date) >>$log

# STEP 14: figures creation
# echo -e 'Creating images' >  $output
echo -e 'Creating images\tStart\t'$(date) >>$log

cd $output_folder
#FIX FILES NAMES AND REMOVE UNUSED FILES
echo -e "Cleaning directory"
rm -f *.CFDGraph.txt
rm -f indels.CFDGraph.txt
rm -r "crispritz_prof"
rm -r "crispritz_targets" #remove targets in online version to avoid memory saturation
#change name to best and alt files
mv $final_res "${output_folder}/$(basename ${output_folder}).bestMerge.txt"
mv $final_res_alt "${output_folder}/$(basename ${output_folder}).altMerge.txt"

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	# ./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}" "var"
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "CFD" 
	# check for variant summary processing on CFD failures
	if [ -s $logerror ]; then 
		printf "ERROR: Variant summary process on CFD report failed!\n" >&2
		exit 1
	fi
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "fewest" 
	# check for summary processing on mm+bulges failures
	if [ -s $logerror ]; then 
		printf "ERROR: Variant summary process on mm+bulges report failed!\n" >&2
		exit 1
	fi
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bMax "${output_folder}" "var" "CRISTA" 
	# check for summary processing on CRISTA failures
	if [ -s $logerror ]; then 
		printf "ERROR: Variant summary process on CRISTA report failed!\n" >&2
		exit 1
	fi
else
	# ./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}" "ref"
	./process_summaries.py $final_res.bestCFD.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "CFD" 
	# check for reference summary processing on CFD failures
	if [ -s $logerror ]; then 
		printf "ERROR: Reference summary process on CFD report failed!\n" >&2
		exit 1
	fi
	./process_summaries.py $final_res.bestmmblg.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "fewest"
	# check for reference summary processing on mm+bulges failures
	if [ -s $logerror ]; then 
		printf "ERROR: Reference summary process on mm+bulges report failed!\n" >&2
		exit 1
	fi 
	./process_summaries.py $final_res.bestCRISTA.txt $guide_file $sampleID $mm $bMax "${output_folder}" "ref" "CRISTA"
	# check for reference summary processing on CRISTA failures
	if [ -s $logerror ]; then 
		printf "ERROR: Reference summary process on CRISTA report failed!\n" >&2
		exit 1
	fi
fi

if ! [ -d "$output_folder/imgs" ]; then
	mkdir "$output_folder/imgs"
fi

if [ "$vcf_name" != "_" ]; then
	cd "$output_folder/imgs"
	while IFS= read -r line || [ -n "$line" ]; do
		for total in $(seq 0 $(expr $mm + $bMax)); do
			# python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution.txt" $total $line
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CFD.txt" $total $line "CFD"
			# check for population distribution on CFD failures
			if [ -s $logerror ]; then 
				printf "ERROR: Population distribution on CFD report failed!\n" >&2
				exit 1
			fi
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_CRISTA.txt" $total $line "CRISTA"
			# check for population distribution on CRISTA failures
			if [ -s $logerror ]; then 
				printf "ERROR: Population distribution on CRISTA report failed!\n" >&2
				exit 1
			fi
			python $starting_dir/populations_distribution.py "${output_folder}/.$(basename ${output_folder}).PopulationDistribution_fewest.txt" $total $line "fewest"
			# check for population distribution on mm+bulges failures
			if [ -s $logerror ]; then 
				printf "ERROR: Population distribution on mm+bulges report failed!\n" >&2
				exit 1
			fi
		done

	done <$guide_file
fi

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	# ./radar_chart_dict_generator.py $guide_file "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CFD"
	# check for radar chart generation on CFD failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on variant CFD report failed!\n" >&2
		exit 1
	fi
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA"
	# check for radar chart generation on CRISTA failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on variant CRISTA report failed!\n" >&2
		exit 1
	fi
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt $sampleID $annotation_file "$output_folder" $ncpus $mm $bMax "fewest" 
	# check for radar chart generation on mm+bulges failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on variant mm+bulges report failed!\n" >&2
		exit 1
	fi
else
	./radar_chart_dict_generator.py $guide_file $final_res.bestCFD.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "CFD"
	# check for radar chart generation on CFD failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on reference CFD report failed!\n" >&2
		exit 1
	fi
	./radar_chart_dict_generator.py $guide_file $final_res.bestCRISTA.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "CRISTA"
	# check for radar chart generation on CRISTA failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on reference CRISTA report failed!\n" >&2
		exit 1
	fi
	./radar_chart_dict_generator.py $guide_file $final_res.bestmmblg.txt $empty_file $annotation_file "$output_folder" $ncpus $mm $bMax "fewest" 
	# check for radar chart generation on mm+bulges failures
	if [ -s $logerror ]; then 
		printf "ERROR: Radar chart generation on reference mm+bulges report failed!\n" >&2
		exit 1
	fi
fi
echo -e 'Creating images\tEnd\t'$(date) >>$log

# STEP 15: targets gene annotation
echo $gene_proximity
echo -e 'Integrating results\tStart\t'$(date) >>$log
echo >>$guide_file

if [ $gene_proximity != "_" ]; then
	genome_version=$(echo ${ref_name} | sed 's/_ref//' | sed -e 's/\n//') #${output_folder}/Params.txt | awk '{print $2}' | sed 's/_ref//' | sed -e 's/\n//')
	echo $genome_version
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${gene_proximity}" $empty_file "${guide_file}" $genome_version "${output_folder}" $empty_dir $starting_dir/ $base_check_start $base_check_end $base_check_set
	# check for gene annotation of primary targets failures
	if [ -s $logerror ]; then 
		printf "ERROR: Gene annotation on primary targets failed!\n" >&2
		exit 1
	fi
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).altMerge.txt" "${gene_proximity}" $empty_file "${guide_file}" $genome_version "${output_folder}" $empty_dir $starting_dir/ $base_check_start $base_check_end $base_check_set
	# check for gene annotation of alternative targets failures
	if [ -s $logerror ]; then 
		printf "ERROR: Gene annotation on alternative targets failed!\n" >&2
		exit 1
	fi
	python $starting_dir/CRISPRme_plots.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/imgs/" &>"${output_folder}/warnings.txt"
	# check for plot failures
	if [ -s $logerror ]; then 
		printf "ERROR: Plots generation failed!\n" >&2
		exit 1
	fi
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
python $starting_dir/db_creation.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/.$(basename ${output_folder})"
# check for database generation failures
if [ -s $logerror ]; then 
	printf "ERROR: Database generation failed!\n" >&2
	exit 1
fi
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
	# check for compression on multiguide integrated results failures
	if [ -s $logerror ]; then 
		printf "ERROR: File compression for multiguide primary targets report failed!\n" >&2
		exit 1
	fi
	zip -j "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/Multiple_spacers+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	# check for compression on multiguide alternative results failures
	if [ -s $logerror ]; then 
		printf "ERROR: File compression for multiguide alternative targets report failed!\n" >&2
		exit 1
	fi
else
	guide_elem=$(head -1 $guide_file)
	mv "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	mv "${output_folder}/$(basename ${output_folder}).altMerge.txt.integrated_results.tsv" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	#generate zipped version for file
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_integrated_results.tsv"
	# check for compression on single guide integrated results failures
	if [ -s $logerror ]; then 
		printf "ERROR: File compression for single guide primary targets report failed!\n" >&2
		exit 1
	fi
	zip -j "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.zip" "${output_folder}/${guide_elem}+${true_pam}_$(basename ${ref_folder})+${vcf_name}_${mm}+${bMax}_all_results_with_alternative_alignments.tsv"
	# check for compression on single alternative results failures
	if [ -s $logerror ]; then 
		printf "ERROR: File compression for single guide alternative targets report failed!\n" >&2
		exit 1
	fi
fi
echo -e "JOB END"

if [ "$email" != "_" ]; then
	python $starting_dir/../pages/send_mail.py $output_folder
fi

#keep log_error but no block visualization
mv $logerror $output_folder/log_error_no_check.txt
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
