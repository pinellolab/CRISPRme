#!/bin/bash

#$1 is CFD targets file
#$2 is genecode annotation
#$3 is empirical data for the guide
#$4 is the guide file
#$5 is genome version (hg38)
#$6 is output directory for the data
#$7 is vcf directory to multi-variant haplotype frequence count

#EXTRACT DIRECTORY FROM REF FILE
dir=$(dirname $1)

#EXAMPLE CALL bash post_process.sh sg1617.best.only_indels.txt gencode.protein_coding.bed sg1617.empiricalresults.tsv guide.txt hg38 sg1617.test

starting_dir=$8
echo 'Preparing files for post processing'
sed -i '/#/d' $1
LC_ALL=C sort -T $dir -k21,21rg $1 -o $1
awk '{print $5"\t"$7"\t"$7+length($3)"\t"$21}' $1 > $1.bed
sort-bed $1.bed > $1.bed.sort
mv $1.bed.sort $1.bed
sort-bed $2 > $2.bed.sort
mv $2.bed.sort $2
echo 'Finding genecode annotation in range with targets'
closest-features --closest --delim "\t" --dist $1.bed $2 > $1.found.bed
echo 'Sorting final annotation results to correspond with original results file'
LC_ALL=C sort -T $dir -k4,4rg $1.found.bed -o $1.found.bed
echo 'Starting integration with empirical data (this may take a while)'
"$starting_dir"./resultIntegrator.py $1 $3 $1.found.bed $4 $6/ true $5 $7
echo 'Removing unnecessary files'
rm $1.bed $1.found.bed $1.redirectFile.out $1.temp.bed
sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\tHighest_CFD_Risk_Score\tHighest_CFD_Absolute_Risk_Score\tMMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref\tMMBLG_CFD_Risk_Score\tMMBLG_CFD_Absolute_Risk_Score" "$1"

