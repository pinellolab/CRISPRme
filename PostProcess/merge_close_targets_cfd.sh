#!/bin/bash

fileIn=$1
fileOut=$2
thresh=$3 #threshold to use in order to merge near targets
sort_criteria=$4  # pivot column to use while sorting targets before merge
sorting_criteria_scoring=$5  # other sorting criteria (score has highest priority)
sorting_criteria=$6   # other sorting criteria

##########ADJUST THESE PARAMETERS BASED ON INPUT FILE##########
#columns start from 1
chrom=5    #column for chromosome
position=7 #column for cluster_position
total=11
true_guide=16
snp_info=19
cfd=21

echo "Sorting file"
# header=$(head -1 $fileIn)
STARTTIME=$(date +%s)
# tail -n +2 $fileIn | LC_ALL=C sort -T$(dirname $(realpath $fileIn)) -k$true_guide,$true_guide -k$chrom,$chrom -k$position,${position}n -o $fileIn.sorted.tmp
ENDTIME=$(date +%s)
echo "Sorting done in $(($ENDTIME - $STARTTIME)) seconds"
# echo -e $header | cat - $fileIn.sorted.tmp > $fileIn.sorted
# rm $fileIn.sorted.tmp
echo "Merging contiguous targets"
if [[ "${sort_pivot}" == "score" ]]; then
    criteria=$sorting_criteria_scoring
else
    criteria=$sorting_criteria
fi

# python remove_contiguous_samples_cfd.py $fileIn.sorted $fileOut $thresh $chrom $position $total $true_guide $snp_info $cfd
python remove_contiguous_samples.py $fileIn $fileOut $thresh $chrom $position $total $true_guide $snp_info $cfd $sort_criteria $criteria
# rm $fileIn.sorted
