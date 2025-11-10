#!/bin/bash

# Script for analyzing REF and ENR search targets with PAM NNN
# The file containing the reference genome search targets is named $REFtargets          -> INPUT $1
# The file containing the enriched genome search targets is named $ENRtargets           -> INPUT $2
# The output prefix for intermediate files is specified by $jobid                       -> INPUT $3
# Ensure that both pam.txt and guides.txt are present in the same directory as this script.
# The .bed annotation file is provided via the variable $annotationfile                 -> INPUT $4
# The dictionary folder, containing .json files, is specified by $dictionaries          -> INPUT $5
# The $mismatch variable defines the number of mismatches used for the search (default: 6),
# while $bulgesDNA and $bulgesRNA define the number of allowed DNA and RNA bulges (default: 2).
# The reference genome folder is indicated by $referencegenome                          -> INPUT $6
# Make sure that the file samplesID.txt (containing the sample IDs) is located in the same directory as this script.
#
# EXAMPLE USAGE:
# ./scriptAnalisiNNN.sh reference.targets.txt enriched.targets.txt output_name hg38_annotations.bed dictionary_hg38/ hg38_ref/
#
# NOTE: If a file with the suffix .total.txt already exists, steps (1) and (2) can be skipped.

# read input arguments
REFtargets=$1
ENRtargets=$2
jobid=$3
annotations=$4
dictionaries=$5
referencegenome=$6
mismatch=$7
bulgesDNA=$8
bulgesRNA=$9
guide_file=${10}
pam_file=${11}
output_folder=${12}

touch $REFtargets.corrected  # initialize reference targets out fname

# remove duplicate targets, extract unique targets and create final summary file
# OUTPUT    $jobid.common_targets.txt -> NOT USED
#           $jobid.semi_common_targets.txt
#           $jobid.unique_targets.txt
./extraction.sh "$REFtargets.corrected" "$ENRtargets" "$jobid" 
rm "$jobid.common_targets.txt"  # remove unnecessary files
rm "$REFtargets.corrected"  # remove unnecessary files

# add PAM creation column to headers
awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' "$jobid.semi_common_targets.txt" >"$jobid.semi_common_targets.minmaxdisr.txt"
rm "$jobid.semi_common_targets.txt"

awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' "$jobid.unique_targets.txt" >"$jobid.unique_targets.pamcreation.txt" # Add pam creation, variant unique, real guide column
rm "$jobid.unique_targets.txt"

cat "$jobid.unique_targets.pamcreation.txt" "$jobid.semi_common_targets.minmaxdisr.txt" >"$jobid.total.txt"
rm "$jobid.unique_targets.pamcreation.txt"
rm "$jobid.semi_common_targets.minmaxdisr.txt"

# targets clustering
./cluster.dict.py "$jobid.total.txt" 'no' 'True' 'True' "$guide_file" 'total' 'orderChr'  # OUTPUT     $jobid.total.cluster.txt

# adjust header (final version)
sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide\tID\tAF\tSNP/' "$jobid.total.cluster.txt"
rm "$jobid.total.txt"

# extract samples
# For this step, instead of using $jobid.total.cluster.txt and $jobid as a single input,
# the script can be executed multiple times using separate files, each containing
# targets from one or more chromosomes extracted from $jobid.total.cluster.txt.
#
# Each of these files MUST:
#   - Include the same header as in $jobid.total.cluster.txt
#   - Contain all targets for a given chromosome (clusters cannot be split across files)
#
# NOTE: When extracting chromosomes, remember to also include alternate contigs such as
#       chr5_KI270791v1_alt from $jobid.total.cluster.txt.
#
# The .CFDGraph.txt files can be summed to obtain the final CFD graph output.
# -> Currently, this step is NOT required for this analysis, but the files SHOULD BE KEPT for consistency.
./analisi_indels_NNN.py "$annotations" "$jobid.total.cluster.txt" "$jobid" "$dictionaries" "$pam_file" "$mismatch" "$referencegenome" "$guide_file" $bulgesDNA $bulgesRNA 
rm "$jobid.total.cluster.txt"  # remove unnecessary files

# sort and adjust processed targets 
echo 'Sorting and adjusting results'  # write to log_verbose
./adjust_cols.py "$jobid.bestCFD_INDEL.txt" 
./adjust_cols.py "$jobid.bestCRISTA_INDEL.txt" 
./adjust_cols.py "$jobid.bestmmblg_INDEL.txt" 

# polish targets found on indels 
./remove_bad_indel_targets.py "$jobid.bestCFD_INDEL.txt" 
./remove_bad_indel_targets.py "$jobid.bestCRISTA_INDEL.txt" 
./remove_bad_indel_targets.py "$jobid.bestmmblg_INDEL.txt" 
