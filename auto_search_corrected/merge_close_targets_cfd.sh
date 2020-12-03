#!/bin/bash


fileIn=$1
fileOut=$2
thresh=$3 #threshold to use in order to merge near targets

##########ADJUST THESE PARAMETERS BASED ON INPUT FILE##########
#columns start from 1
chrom=5 #column for chromosome
position=6 #column for position
total=11
snp_info=19
cfd=21

echo "Sorting file"
header=$(head -1 $fileIn)
STARTTIME=$(date +%s)
tail -n +2 $fileIn | LC_ALL=C sort -k$chrom,$chrom -k$position,${position}n -o $fileIn.sorted.tmp
ENDTIME=$(date +%s)
echo "Sorting done in $(($ENDTIME - $STARTTIME)) seconds"
echo -e $header | cat - $fileIn.sorted.tmp > $fileIn.sorted
rm $fileIn.sorted.tmp
echo "Merging targets"
python3 remove_contiguous_samples_cfd.py $fileIn.sorted $fileOut $thresh $chrom $position $total $snp_info $cfd
rm $fileIn.sorted
