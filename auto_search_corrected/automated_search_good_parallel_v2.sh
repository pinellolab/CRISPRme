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

merge_t=${11}

output_folder=$(realpath ${12})

starting_dir=${13}
ncpus=${14}
echo "CPU used: $ncpus" 

ref_name=$(basename $1)
#folder_of_folders=$(dirname $1)
vcf_name=$(basename $2)
guide_name=$(basename $3)
pam_name=$(basename $4)
annotation_name=$(basename $5)

if [ "$vcf_name" != "_" ]; then
	log="log_complete_search_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.txt"
else
	log="log_complete_search_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.txt"
fi

touch $output_folder/$log

echo "Complete search Start:" $(date +%F-%T) >> $output_folder/$log
echo "########################################" >> $output_folder/$log
echo "INPUT:" >> $output_folder/$log
echo "Reference - $ref_name" >> $output_folder/$log
echo "VCFs - $vcf_name" >> $output_folder/$log
echo "Guide - $guide_name" >> $output_folder/$log
echo "PAM - $pam_name" >> $output_folder/$log
echo "Annotation - $annotation_name" >> $output_folder/$log
echo "########################################" >> $output_folder/$log
#echo "Output - $output_folder" >> $output_folder/$log

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
		./create_dict.py "$vcf_folder" "$vcf_name" "$output_folder" "$log"&
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
	else
		echo "Search INDELs already done"
	fi
fi

while kill "-0" $pid_search_ref &>/dev/null; do
	echo "Waiting for search genome reference"
	sleep 300
done
if [ -f "$output_folder/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	mv "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" "$output_folder/crispritz_targets"
fi
if ! [ -d "$output_folder/crispritz_profiles" ]; then
	mkdir $output_folder/crispritz_profiles
fi
mv $output_folder/*profile* $output_folder/crispritz_profiles/ > /dev/null 2>&1

echo "Search Reference End: "$(date +%F-%T) >> $output_folder/$log	

while kill "-0" $pid_dicts &>/dev/null; do
	echo "Waiting for dictionary creation before SNP analysis"
	sleep 300
done


cd "$starting_dir"

echo "Start post-analysis"

if [ "$vcf_name" != "_" ]; then
	echo "Post-analysis SNPs Start: "$(date +%F-%T) >> $output_folder/$log	
	final_res="$output_folder/final_results_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
	final_res_alt="$output_folder/final_results_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	else
		echo "This result already exists. $final_res"
		exit
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	else
		echo "This result already exists. $final_res_alt"
		exit
	fi
	
	./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID $dict_folder $final_res $final_res_alt $ncpus
	
	echo "Post-analysis SNPs End: "$(date +%F-%T) >> $output_folder/$log	
	for key in "${real_chroms[@]}"
	do
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" 
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res"
		rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
		rm "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
	done
else

	final_res="$output_folder/final_results_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
	final_res_alt="$output_folder/final_results_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	if ! [ -f "$final_res" ]; then
		touch "$final_res"
	else
		echo "This result already exists. $final_res"
		exit
	fi
	if ! [ -f "$final_res_alt" ]; then
		touch "$final_res_alt"
	else
		echo "This result already exists. $final_res_alt"
		exit
	fi
	for key in "${real_chroms[@]}"
	do
		echo "Processing $key"
		LC_ALL=C grep -P "$key\t" "$output_folder/crispritz_targets/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" > "$output_folder/crispritz_targets/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		touch "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		./scriptAnalisiNNN_v3.sh "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key" "$annotation_file" "_" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder"
		rm "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		rm "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		header=$(head -1 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt")
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >> "$final_res" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >> "$final_res" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"
		
	done

fi

echo "Adding header to files"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res"
#sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res_alt"

if [ "$vcf_name" != "_" ]; then
	echo "SNPs analysis ended. Starting INDELs analysis"
	cd "$starting_dir"
	
	echo "Post-analysis INDELs Start: "$(date +%F-%T) >> $output_folder/$log	
	./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID "$output_folder/log_indels_$vcf_name" $final_res $final_res_alt $ncpus
	echo "Post-analysis INDELs End: "$(date +%F-%T) >> $output_folder/$log	
	for key in "${array_fake_chroms[@]}"
	do
		tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt" >> "$final_res" 
		tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt" >> "$final_res" 
		rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
		rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	done

fi

echo "Merging Close Targets Start: "$(date +%F-%T) >> $output_folder/$log	
./merge_close_targets_cfd.sh $final_res $final_res.trimmed $merge_t
mv $final_res.trimmed $final_res
mv $final_res.trimmed.discarded_samples $final_res_alt

#./merge_close_targets_cfd.sh $final_res_alt $final_res_alt.trimmed $merge_t
#rm $final_res_alt
echo "Merging Close Targets End: "$(date +%F-%T) >> $output_folder/$log	

echo "Merging Alternative Chromosomes Start: "$(date +%F-%T) >> $output_folder/$log	
./merge_alt_chr.sh $final_res $final_res.chr_merged
#rm $final_res.trimmed

#./merge_alt_chr.sh $final_res_alt.trimmed $final_res_alt.trimmed.chr_merged
#rm $final_res_alt.trimmed
echo "Merging Alternative Chromosomes End: "$(date +%F-%T) >> $output_folder/$log	

mv $final_res.chr_merged $final_res
#mv $final_res_alt.trimmed.chr_merged $final_res_alt

sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res_alt"

echo "Cleaning directory"

if ! [ -d "$output_folder/cfd_graphs" ]; then
	mkdir $output_folder/cfd_graphs
fi
./assemble_cfd_graphs.py $output_folder
mv $output_folder/snps.CFDGraph.txt $output_folder/cfd_graphs
mv $output_folder/indels.CFDGraph.txt $output_folder/cfd_graphs

echo "Automated Search End: "$(date +%F-%T) >> $output_folder/$log
echo "AUTOMATED SEARCH END"

: '

echo "Adjusting Results Start: "$(date +%F-%T) >> $output_folder/$log	
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
echo "Adjusting Results End: "$(date +%F-%T) >> $output_folder/$log	

'