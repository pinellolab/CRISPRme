#!/bin/bash

dir=$(dirname $1)
fileIn=$1
fileOut=$2

STARTTIME=$(date +%s)

chroms=($(cut -f 5 $fileIn | tail -n +2 | LC_ALL=C grep -v "_" | LC_ALL=C sort -T "$dir" | uniq))

head -1 $fileIn >$fileOut

for chrom in ${chroms[@]}; do

    echo $chrom
    #awk -v "key=$chrom" '$5 == key {print($0)}' $fileIn > $fileIn.$chrom.ref
    grep -F -w "$chrom" $fileIn >$fileIn.$chrom.ref
    cut -f 3 $fileIn.$chrom.ref | LC_ALL=C sort -T "$dir" | uniq >$fileIn.$chrom.ref.targets
    LC_ALL=C grep -F "${chrom}_" $fileIn >$fileIn.$chrom.alt # | awk '$5 ~ "_" {print($0)}'
    LC_ALL=C grep -v -F -f $fileIn.$chrom.ref.targets $fileIn.$chrom.alt >$fileIn.$chrom.merged
    rm $fileIn.$chrom.ref.targets
    rm $fileIn.$chrom.alt
    cat $fileIn.$chrom.ref $fileIn.$chrom.merged >>$fileOut
    rm $fileIn.$chrom.merged
    rm $fileIn.$chrom.ref

done

ENDTIME=$(date +%s)
echo "Merging done in $(($ENDTIME - $STARTTIME)) seconds"
