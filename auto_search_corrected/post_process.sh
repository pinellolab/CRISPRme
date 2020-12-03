#!/bin/bash

#$1 is CFD targets file
#$2 is genecode annotation
#$3 is empirical data for the guide
#$4 is the guide (e.g. CTAACAGTTGCTTTTATCACNNN)
#$5 is output name for integrated data

#EXAMPLE CALL bash post_process.sh sg1617.best.only_indels.txt gencode.protein_coding.bed sg1617.empiricalresults.tsv guide.txt hg38 sg1617.test

starting_dir=$7
echo 'Preparing files for post processing'
sed -i '/#/d' $1
LC_ALL=C sort -T . -k21,21rg $1 -o $1
awk '{print $5"\t"$7"\t"$7+length($3)"\t"$21}' $1 > $1.bed
sort-bed $1.bed > $1.bed.sort
mv $1.bed.sort $1.bed
sort-bed $2 > $2.bed.sort
mv $2.bed.sort $2
echo 'Finding genecode annotation in range with targets'
closest-features --closest --delim "\t" --dist $1.bed $2 > $1.found.bed
echo 'Sorting final annotation results to correspond with original results file'
LC_ALL=C sort -T . -k4,4rg $1.found.bed -o $1.found.bed
echo 'Starting integration with empirical data (this may take a while)'
"$starting_dir"./resultIntegrator.py $1 $3 $1.found.bed $4 $6/ true $5
echo 'Removing unnecessary files'
rm $1.bed $1.found.bed
