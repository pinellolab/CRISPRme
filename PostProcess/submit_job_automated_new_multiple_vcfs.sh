#!/bin/bash

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
echo -e "CPU used: $ncpus" 


log=$output_folder/log.txt
touch $log

output=$output_folder/output.txt
touch $output

rm $output_folder/queue.txt
#for vcf_f in "${vcf_list[@]}"; 
if [ $2 == "_" ]; then
	echo -e "_" >> $output_folder/tmp_list_vcf.txt
	vcf_list=$output_folder/tmp_list_vcf.txt
fi
echo "" >> $vcf_list
if [ $6 != "_" ]; then
	echo "" >> $6
fi

while read vcf_f;
do
	if [ -z "$vcf_f" ]; then
		continue
	fi
	vcf_folder=$(realpath ${vcf_f})
	ref_name=$(basename $1)
	#folder_of_folders=$(dirname $1)
	vcf_name=$(basename $vcf_f)
	echo "STARTING ANALYSIS FOR $vcf_name"
	# echo $vcf_name
	guide_name=$(basename $3)
	pam_name=$(basename $4)
	annotation_name=$(basename $5)


	echo -e 'Job\tStart\t'$(date) > $log

	unset real_chroms
	declare -a real_chroms
	for file_chr in "$ref_folder"/*.fa
	do
		file_name=$(basename $file_chr)
		chr=$(echo -e $file_name | cut -f 1 -d'.')
		echo -e "$chr"
		real_chroms+=("$chr")	
	done

	if [ "$vcf_name" != "_" ]; then
		unset array_fake_chroms
		declare -a array_fake_chroms
		for file_chr in "$vcf_folder"/*.vcf.gz
		do
			file_name=$(basename $file_chr)
			# file_name=$(basename $file_chr)
			IFS='.' read -ra ADDR <<< $file_name
			for i in "${ADDR[@]}"; 
			do
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
			echo -e 'Add-variants\tStart\t'$(date) >> $log
			echo -e "Adding variants"
			crispritz.py add-variants "$vcf_folder/" "$ref_folder/" "true"
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
			echo -e 'Add-variants\tEnd\t'$(date) >> $log
			echo -e 'Indexing Indels\tStart\t'$(date) >> $log
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/variants_genome/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus
			echo -e 'Indexing Indels\tEnd\t'$(date) >> $log
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
			echo -e 'Indexing Indels\tStart\t'$(date) >> $log
			${starting_dir}/./pool_index_indels.py "$current_working_directory/Genomes/${ref_name}+${vcf_name}_INDELS/" "$pam_file" $true_pam $ref_name $vcf_name $ncpus
			echo -e 'Indexing Indels\tEnd\t'$(date) >> $log
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
					echo -e 'Index-genome Reference\tStart\t'$(date) >> $log	
					echo -e 'Indexing_Reference' > $output
					echo -e "Indexing reference genome"
					crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus
					pid_index_ref=$!
					echo -e 'Index-genome Reference\tEnd\t'$(date) >> $log	
					idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
				else
					echo -e "Reference Index already present"
					idx_ref="$current_working_directory/genome_library/${true_pam}_1_${ref_name}"
				fi
			else
				echo -e 'Index-genome Reference\tStart\t'$(date) >> $log	
				echo -e 'Indexing_Reference' > $output
				echo -e "Indexing reference genome"
				crispritz.py index-genome "$ref_name" "$ref_folder/" "$pam_file" -bMax $bMax -th $ncpus
				pid_index_ref=$!
				echo -e 'Index-genome Reference\tEnd\t'$(date) >> $log	
				idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
			fi
		else
			echo -e "Reference Index already present"
			echo -e 'Index-genome Reference\tEnd\t'$(date) >> $log	
			idx_ref="$current_working_directory/genome_library/${true_pam}_2_${ref_name}"
		fi
	else
		echo -e "Reference Index already present"
		echo -e 'Index-genome Reference\tEnd\t'$(date) >> $log	
		idx_ref="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}"
	fi


	if [ "$vcf_name" != "_" ]; then
		if ! [ -d "$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}" ]; then
			if ! [ -d "$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}" ]; then
				if ! [ $bMax -gt 1 ]; then
					if ! [ -d "$current_working_directory/genome_library/${true_pam}_1_${ref_name}+${vcf_name}" ]; then
						echo -e 'Index-genome Variant\tStart\t'$(date) >> $log	
						echo -e 'Indexing_Enriched' > $output
						echo -e "Indexing variant genome"
						crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus #${ref_folder%/}+${vcf_name}/
						pid_index_var=$!
						echo -e 'Index-genome Variant\tEnd\t'$(date) >> $log	
						idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
					else
						echo -e "Variant Index already present"
						idx_var="$current_working_directory/genome_library/${true_pam}_1_${ref_name}+${vcf_name}"
					fi
				else
					echo -e 'Index-genome Variant\tStart\t'$(date) >> $log	
					echo -e 'Indexing_Enriched' > $output
					echo -e "Indexing variant genome"
					crispritz.py index-genome "${ref_name}+${vcf_name}" "$current_working_directory/Genomes/${ref_name}+${vcf_name}/" "$pam_file" -bMax $bMax -th $ncpus
					pid_index_ref=$!
					echo -e 'Index-genome Variant\tEnd\t'$(date) >> $log	
					idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
				fi
			else
				echo -e "Variant Index already present"
				echo -e 'Index-genome Variant\tEnd\t'$(date) >> $log
				idx_var="$current_working_directory/genome_library/${true_pam}_2_${ref_name}+${vcf_name}"
			fi
		else
			echo -e "Variant Index already present"
			echo -e 'Index-genome Variant\tEnd\t'$(date) >> $log
			idx_var="$current_working_directory/genome_library/${true_pam}_${bMax}_${ref_name}+${vcf_name}"
		fi	
	fi

	cd "$output_folder"
	echo $idx_ref
	if ! [ -f "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		echo -e 'Search Reference\tStart\t'$(date) >> $log	
		echo -e 'Search Reference' >  $output
		crispritz.py search $idx_ref "$pam_file" "$guide_file" "${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th $(expr $ncpus / 4) &
		pid_search_ref=$!
	else
		echo -e "Search for reference already done"
	fi

	if [ "$vcf_name" != "_" ]; then
		if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e 'Search Variant\tStart\t'$(date) >> $log	
			echo -e 'Search Variant' >  $output
			crispritz.py search "$idx_var" "$pam_file" "$guide_file" "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}" -index -mm $mm -bDNA $bDNA -bRNA $bRNA -t  -th $(expr $ncpus / 4) -var
			mv "${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			echo -e 'Search Variant\tEnd\t'$(date) >> $log	
		else
			echo -e "Search for variant already done"
		fi
		
		if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
			echo -e "Search INDELs Start"
			echo -e 'Search INDELs\tStart\t'$(date) >> $log	
			cd $starting_dir
			./pool_search_indels.py "$ref_folder" "$vcf_folder" "$vcf_name" "$guide_file" "$pam_file" $bMax $mm $bDNA $bRNA "$output_folder" $true_pam "$current_working_directory/"
			mv "$output_folder/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
			echo -e "Search INDELs End"
			echo -e 'Search INDELs\tEnd\t'$(date) >> $log
		else
			echo -e "Search INDELs already done"
		fi
	fi

	while kill "-0" $pid_search_ref &>/dev/null; do
		echo -e "Waiting for search genome reference"
		sleep 100
	done
	echo -e 'Search Reference\tEnd\t'$(date) >> $log	

	if [ -f "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
		mv "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
	fi

	if ! [ -d "$output_folder/crispritz_prof" ]; then
		mkdir $output_folder/crispritz_prof
	fi
	mv $output_folder/*profile* $output_folder/crispritz_prof/ &>/dev/null

	cd "$starting_dir"

	echo -e "Start post-analysis"

	echo -e 'Post analysis' >  $output
	if [ "$vcf_name" != "_" ]; then
		echo -e 'Post-analysis SNPs\tStart\t'$(date) >> $log	
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then
			touch "$final_res_alt"
		fi
		
		./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $dict_folder $final_res $final_res_alt $ncpus
		
		echo -e 'Post-analysis SNPs\tEnd\t'$(date) >> $log	
		for key in "${real_chroms[@]}"
		do
			echo "Concatenating $key"
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
			# tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res_alt" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
			# rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
		done
	else
		echo -e 'Post-analysis\tStart\t'$(date) >> $log	
		final_res="$output_folder/final_results_$(basename ${output_folder}).bestMerge.txt"
		final_res_alt="$output_folder/final_results_$(basename ${output_folder}).altMerge.txt"
		if ! [ -f "$final_res" ]; then
			touch "$final_res"
		fi
		if ! [ -f "$final_res_alt" ]; then
			touch "$final_res_alt"
		fi

		./pool_post_analisi_snp.py $output_folder $ref_folder "_" $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "_" $final_res $final_res_alt $ncpus
		
		for key in "${real_chroms[@]}"
		do
			tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
			# tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res_alt" #"$output_folder/${ref_name}+${vcf_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
			rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
			# rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
		done
		echo -e 'Post-analysis\tEnd\t'$(date) >> $log	

	fi

	if [ "$vcf_name" != "_" ]; then
		echo -e "SNPs analysis ended. Starting INDELs analysis"
		cd "$starting_dir"
		
		echo -e 'Post-analysis INDELs\tStart\t'$(date) >> $log	
		./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file "$current_working_directory/Dictionaries/" $final_res $final_res_alt $ncpus
		echo -e 'Post-analysis INDELs\tEnd\t'$(date) >> $log	
		for key in "${array_fake_chroms[@]}"
		do	
			echo "Concatenating $key"
			tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt" >> "$final_res" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
			tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt" >> "$final_res_alt" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
			rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
			rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
		done

	fi
done < $vcf_list
echo -e "Adding header to files"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res"

while read samples;
do
	if [ -z "$samples" ]; then
		continue
	fi
	tail -n +2 $samples >> "$output_folder/sampleID.txt"
done < $sampleID
sed -i 1i"#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tGENDER" "$output_folder/sampleID.txt"

sampleID=$output_folder/sampleID.txt


echo -e 'Merging targets' >  $output
echo -e 'Merging Close Targets\tStart\t'$(date) >> $log	
./merge_close_targets_cfd.sh $final_res $final_res.trimmed $merge_t
mv $final_res.trimmed $final_res
mv $final_res.trimmed.discarded_samples $final_res_alt

echo -e 'Merging Close Targets\tEnd\t'$(date) >> $log	

echo -e 'Merging Alternative Chromosomes\tStart\t'$(date) >> $log	
./merge_alt_chr.sh $final_res $final_res.chr_merged

echo -e 'Merging Alternative Chromosomes\tEnd\t'$(date) >> $log	

mv $final_res.chr_merged $final_res

sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res_alt"

./annotate_final_results.py $final_res $annotation_file $final_res.annotated
./annotate_final_results.py $final_res_alt $annotation_file $final_res_alt.annotated

mv $final_res.annotated $final_res
mv $final_res_alt.annotated $final_res_alt

echo -e "Cleaning directory"

if ! [ -d "$output_folder/cfd_graphs" ]; then
	mkdir $output_folder/cfd_graphs
fi
./assemble_cfd_graphs.py $output_folder
mv $output_folder/snps.CFDGraph.txt $output_folder/cfd_graphs
mv $output_folder/indels.CFDGraph.txt $output_folder/cfd_graphs

echo -e 'Creating images' >  $output
echo -e 'Creating images\tStart\t'$(date) >> $log	
echo -e "Adding risk score"
./add_risk_score.py $final_res $final_res.risk
mv "$final_res.risk" "${output_folder}/$(basename ${output_folder}).bestMerge.txt"
./add_risk_score.py $final_res_alt $final_res_alt.risk
mv "$final_res_alt.risk" "${output_folder}/$(basename ${output_folder}).altMerge.txt"
echo -e "Risk score added"

cd $output_folder
rm -r "cfd_graphs"
rm -r "crispritz_prof"
# rm -r "crispritz_targets"
rm $final_res
rm $final_res_alt

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}/$(basename ${output_folder})" "var"
else
	./process_summaries.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $guide_file $sampleID $mm $bMax "${output_folder}/$(basename ${output_folder})" "ref"
fi	


if ! [ -d "$output_folder/imgs" ]; then
	mkdir "$output_folder/imgs"
fi

if [ "$vcf_name" != "_" ]; then
	cd "$output_folder/imgs"
	while IFS= read -r line || [ -n "$line" ]; do
		for total in $(seq 0 $(expr $mm + $bMax))
		do
			python $starting_dir/populations_distribution.py "${output_folder}/$(basename ${output_folder}).PopulationDistribution.txt" $total $line
		done
		
	done < $guide_file
fi

cd $starting_dir
if [ "$vcf_name" != "_" ]; then
	./radar_chart.py $guide_file "${output_folder}/$(basename ${output_folder}).bestMerge.txt" $sampleID $annotation_file "$output_folder/imgs" $ncpus
else
	echo -e "dummy_file" > dummy.txt
	./radar_chart.py $guide_file "${output_folder}/$(basename ${output_folder}).bestMerge.txt" dummy.txt $annotation_file "$output_folder/imgs" $ncpus
	rm dummy.txt
fi
echo -e 'Creating images\tEnd\t'$(date) >> $log	

# if [ "$vcf_name" != "_" ]; then
# 	cp $sampleID $output_folder/sampleID.txt
# fi

echo -e 'Building database'
echo -e 'Creating database\tStart\t'$(date) >> $log
if [ -f "${output_folder}/$(basename ${output_folder}).db" ]; then
	rm "${output_folder}/$(basename ${output_folder}).db"
fi
python $starting_dir/db_creation.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${output_folder}/$(basename ${output_folder})"
echo -e 'Creating database\tEnd\t'$(date) >> $log

echo $gene_proximity
echo -e 'Integrating results\tStart\t'$(date) >> $log
if [ $gene_proximity != "_" ]; then
	touch "${output_folder}/dummy.txt"
	genome_version=$(echo ${ref_name} | sed 's/_ref//' | sed -e 's/\n//') #${output_folder}/Params.txt | awk '{print $2}' | sed 's/_ref//' | sed -e 's/\n//')
	echo $genome_version
	bash $starting_dir/post_process.sh "${output_folder}/$(basename ${output_folder}).bestMerge.txt" "${gene_proximity}" "${output_folder}/dummy.txt" "${guide_file}" $genome_version "${output_folder}" "vuota"
	rm "${output_folder}/dummy.txt"
	python $starting_dir/CRISPRme_plots.py "${output_folder}/$(basename ${output_folder}).bestMerge.txt.integrated_results.tsv" "${output_folder}/imgs/"
fi
echo -e 'Integrating results\tEnd\t'$(date) >> $log
echo -e 'Job\tDone\t'$(date) >> $log
echo -e 'Job End' >  $output
echo -e "JOB END"

if [ "$email" != "" ]; then
	python $starting_dir/../pages/send_mail.py $output_folder
fi
