#!/bin/bash

set -e  # trace all failures

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
annotation_name=$(basename $5)

if [ "$vcf_name" != "_" ]; then
	log="log_post_analysis_only_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.txt"
else
	log="log_post_analysis_only_${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.txt"
fi

touch $output_folder/$log

echo "Post-analysis Start:" $(date +%F-%T) >>$output_folder/$log
echo "########################################" >>$output_folder/$log
echo "INPUT:" >>$output_folder/$log
echo "Reference - $ref_name" >>$output_folder/$log
echo "VCFs - $vcf_name" >>$output_folder/$log
echo "Guide - $guide_name" >>$output_folder/$log
echo "PAM - $pam_name" >>$output_folder/$log
echo "Annotation - $annotation_name" >>$output_folder/$log
echo "########################################" >>$output_folder/$log

declare -a real_chroms
for file_chr in "$ref_folder"/*.fa; do
	file_name=$(basename $file_chr)
	chr=$(echo $file_name | cut -f 1 -d'.')
	echo "$chr"
	real_chroms+=("$chr")
done

if [ "$vcf_name" != "_" ]; then
	declare -a array_fake_chroms
	for file_chr in "$vcf_folder"/*.vcf.gz; do
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

cd "$starting_dir"

echo "Start post-analysis"

if ! [ -f "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	echo "Ref targets not found, please launch search-only"
	exit
fi
if ! [ -f "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	echo "Variant targets not found, please launch search-only"
	exit
fi
if ! [ -f "$output_folder/crispritz_targets/indels_${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" ]; then
	echo "INDELs targets not found, please launch search-only"
	exit
fi

if [ "$vcf_name" != "_" ]; then
	dict_folder="$output_folder/dictionaries_$vcf_name/"
	echo "Post-analysis SNPs Start: "$(date +%F-%T) >>$output_folder/$log
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

	./pool_post_analisi_snp.py $output_folder $ref_folder $vcf_name $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID $dict_folder $final_res $final_res_alt $ncpus || {
		echo "CRISPRme ERROR: indels postprocessing failed (script: ${0} line $((LINENO-1)))" >&2
		exit 1
	}

	echo "Post-analysis SNPs End: "$(date +%F-%T) >>$output_folder/$log

	for key in "${real_chroms[@]}"; do
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >>"$final_res"
		tail -n +2 "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >>"$final_res"
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
	for key in "${real_chroms[@]}"; do
		echo "Processing $key"
		LC_ALL=C grep -F -w $key "$output_folder/crispritz_targets/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt" >"$output_folder/crispritz_targets/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		touch "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		./scriptAnalisiNNN_v3.sh "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key" "${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key" "$annotation_file" "_" "$ref_folder" $mm $bDNA $bRNA "$guide_file" "$pam_file" "$sampleID" "$output_folder" || {
			echo "CRISPRme ERROR: analysis failed (script: ${0} line $((LINENO-1)))" >&2
			exit 1
		}
		rm "$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		rm "$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt.$key"
		header=$(head -1 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt")
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt" >>"$final_res" #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt" >>"$final_res"  #"$output_folder/${ref_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.bestMerge.txt"
		rm "$output_folder/${ref_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_$key.altMerge.txt"

	done

fi

echo "Adding header to files"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res"
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref" "$final_res_alt"

if [ "$vcf_name" != "_" ]; then
	echo "SNPs analysis ended. Starting INDELs analysis"
	cd "$starting_dir"

	echo "Post-analysis INDELs Start: "$(date +%F-%T) >>$output_folder/$log
	./pool_post_analisi_indel.py $output_folder $ref_folder $vcf_folder $guide_file $mm $bDNA $bRNA $annotation_file $pam_file $sampleID "$output_folder/log_indels_$vcf_name" $final_res $final_res_alt $ncpus || {
		echo "CRISPRme ERROR:indels analysis failed (script: ${0} line $((LINENO-1)))" >&2
		exit 1
	}
	echo "Post-analysis INDELs End: "$(date +%F-%T) >>$output_folder/$log
	for key in "${array_fake_chroms[@]}"; do
		tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt" >>"$final_res" #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.bestCFD.txt.tmp"
		tail -n +2 "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt" >>"$final_res"  #"$output_folder/${fake_chr}_${guide_name}_${mm}_${bDNA}_${bRNA}.altCFD.txt.tmp"
		rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.bestMerge.txt"
		rm "$output_folder/${key}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}.altMerge.txt"
	done

fi

cd "$starting_dir"

echo "Merging Close Targets Start: "$(date +%F-%T) >>$output_folder/$log
./merge_close_targets_cfd.sh $final_res $final_res.trimmed $merge_t || {
	echo "CRISPRme ERROR: CFD targets merge failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
mv $final_res.trimmed $final_res
mv $final_res.trimmed.discarded_samples $final_res_alt

#./merge_close_targets_cfd.sh $final_res_alt $final_res_alt.trimmed $merge_t
#rm $final_res_alt
echo "Merging Close Targets End: "$(date +%F-%T) >>$output_folder/$log

echo "Merging Alternative Chromosomes Start: "$(date +%F-%T) >>$output_folder/$log
./merge_alt_chr.sh $final_res $final_res.chr_merged  || {
	echo "CRISPRme ERROR: alternative targets merge failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
#rm $final_res.trimmed

#./merge_alt_chr.sh $final_res_alt.trimmed $final_res_alt.trimmed.chr_merged
#rm $final_res_alt.trimmed
echo "Merging Alternative Chromosomes End: "$(date +%F-%T) >>$output_folder/$log

mv $final_res.chr_merged $final_res
#mv $final_res_alt.trimmed.chr_merged $final_res_alt

sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res"
sed -i '1 s/^.*$/#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' "$final_res_alt"

echo "Cleaning directory"

if ! [ -d "$output_folder/cfd_graphs" ]; then
	mkdir $output_folder/cfd_graphs
fi
./assemble_cfd_graphs.py $output_folder || {
	echo "CRISPRme ERROR: CFD graph creation failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
mv $output_folder/snps.CFDGraph.txt $output_folder/cfd_graphs
mv $output_folder/indels.CFDGraph.txt $output_folder/cfd_graphs

echo "Post-analysis End: "$(date +%F-%T) >>$output_folder/$log
echo "POST-ANALYSIS END"
