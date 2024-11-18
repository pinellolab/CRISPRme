#!/bin/bash

dir=$(dirname $1)
fileIn=$1
fileOut=$2

STARTTIME=$(date +%s)

chroms=($(awk 'NR>1 && !/_/ {print $5}' $fileIn | LC_ALL=C sort -T "$dir" | uniq))

head -1 $fileIn >$fileOut

for chrom in ${chroms[@]}; do

    echo $chrom
    # awk "/${chrom}\t/" test.targets.txt >$fileIn.$chrom.ref
    grep -F -w "$chrom" $fileIn >$fileIn.$chrom.ref
    cut -f 3 $fileIn.$chrom.ref | LC_ALL=C sort -T "$dir" | uniq >$fileIn.$chrom.ref.targets
    awk -v chrom="$chrom" '$0 ~ chrom"_" {print($0)}' $fileIn >$fileIn.$chrom.alt
    awk 'NR==FNR{a[$0];next} !($0 in a)' $fileIn.$chrom.ref.targets $fileIn.$chrom.alt >$fileIn.$chrom.merged
    rm $fileIn.$chrom.ref.targets
    rm $fileIn.$chrom.alt
    cat $fileIn.$chrom.ref $fileIn.$chrom.merged >>$fileOut
    rm $fileIn.$chrom.merged
    rm $fileIn.$chrom.ref

done

ENDTIME=$(date +%s)
echo "Merging done in $(($ENDTIME - $STARTTIME)) seconds"
