#!/bin/bash

# Script per l'analisi dei targets della ricerca REF e ENR con PAM NNN
# Il file dei targets della ricerca sul genoma reference si chiama $REFtargets  -> INPUT $1
# Il file dei targets della ricerca sul genoma enriched si chiama $ENRtargets   -> INPUT $2
# Output name dei file di intermezzo si chiama $jobid                           -> INPUT $3
# Assicurarsi che il file della pam (pam.txt) e delle guide (guides.txt) siano
# presenti all'interno della cartella contentente lo script
# Il filde .bed delle annotazioni, variabile $annotationfile                    -> INPUT $4
# Cartella dei dizionari, che contiene i file .json, variabile $dictionaries    -> INPUT $5
# $mismatch indica il numero di mismatches usati per la ricerca, default 6,
# $bulgesDNA e $bulgesRNA indicano il num di bulges usati nella ricerca, default 2
# cartella del genoma reference, indicata con $referencegenome                  -> INPUT $6
# Assicurarsi che il file con gli ID dei sample (samplesID.txt) sia nella cartella
# dello script

#ESEMPIO CHIAMATA
# ./scriptAnalisiNNN.sh reference.targets.txt enriched.targets.txt nome_output hg38_annotations.bed dictionary_hg38/ hg38_ref/

#NOTA se è già presente un file con .total.txt si possono commentare gli step 1) e 2)

REFtargets=$1
ENRtargets=$2
jobid=$3
annotationfile=$4
dictionaries=$5
referencegenome=$6

mismatch=$7
bulgesDNA=$8
bulgesRNA=$9

guide_file=${10}
pam_file=${11}
# sampleID=${12}

output_folder=${12}

echo $jobid
# 1) Rimozione duplicati, estrazione semicommon e unique e creazione file total
#echo 'Creazione file .total.txt'
./extraction.sh $REFtargets $ENRtargets $jobid # OUTPUT    $jobid.common_targets.txt -> Non usato
#           $jobid.semi_common_targets.txt
#           $jobid.unique_targets.txt
rm $jobid.common_targets.txt

# 2) Creazione colonne PAM creation etc
awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' $jobid.semi_common_targets.txt >$jobid.semi_common_targets.minmaxdisr.txt
rm $jobid.semi_common_targets.txt

awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' $jobid.unique_targets.txt >$jobid.unique_targets.pamcreation.txt #Add pam creation, variant unique, real guide column
rm $jobid.unique_targets.txt

cat $jobid.unique_targets.pamcreation.txt $jobid.semi_common_targets.minmaxdisr.txt >$jobid.total.txt
rm $jobid.unique_targets.pamcreation.txt
rm $jobid.semi_common_targets.minmaxdisr.txt

#echo 'Creazione cluster del file .total.txt'
# 3) Clustering
./cluster.dict.py $jobid.total.txt 'no' 'True' 'True' "$guide_file" 'total' 'orderChr' # OUTPUT     $jobid.total.cluster.txt
rm $jobid.total.txt

#sed -i ':a;N;$!ba;s/\n/\tn\tn\tn\n/g' $jobid.total.cluster.txt
#sed -i '$s/$/\tn\tn\tn/g' $jobid.total.cluster.txt

#sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide\trsID\tAF\tSNP position/' $jobid.total.cluster.txt

# 4) Estrazione del samples

# Per questa parte, al posto di $jobid.total.cluster.txt e $jobid, si può chiamare lo script più volte mettendo i file divisi per uno o più
# cromosomi, estratti partendo dal file $jobid.total.cluster.txt
# Ognuno di questi file DEVE avere l'header preso dal file $jobid.total.cluster.txt
# Ognuno di questi file DEVE contenere tutti i target di un cromosoma (non posso avere cluster spezzati in due file diversi)
# NOTA greppare da $jobid.total.cluster.txt anche i cromosomi del tipo chr5_KI270791v1_alt

# ESEMPIO so ho 3 file del tipo chr1-5.txt, chr 6-15.txt e chr 16-22XYM.txt
# farò queste chiamate
# python AnnotatorAllTargets.py $annotationfile chr1-5.txt chr1-5 $dictionaries pam.txt $mismatch $referencegenome guides.txt $bulgesDNA $bulgesRNA samplesID.txt
# python AnnotatorAllTargets.py $annotationfile chr6-15.txt chr6-15 $dictionaries pam.txt $mismatch $referencegenome guides.txt $bulgesDNA $bulgesRNA samplesID.txt
# python AnnotatorAllTargets.py $annotationfile chr16-22XYM.txt chr16-22XYM $dictionaries pam.txt $mismatch $referencegenome guides.txt $bulgesDNA $bulgesRNA samplesID.txt
# OUTPUT    chr1-5.bestCFD.txt      chr6-15.bestCFD.txt     chr16-22XYM.bestCFD.txt
#           chr1-5.CFDGraph.txt     chr6-15.CFDGraph.txt    chr16-22XYM.CFDGraph.txt
# I file .bestCFD.txt si possono unire (attenzione che avrò 3 header all'interno del file) per ottenere il file completo
# I file .CFDGraph.txt vanno sommati tra loro per ottenere il file finale -> AL MOMENTO NON NECESSARIO PER QUESTA ANALISI (TENERE COMUNQUE I FILE)

#echo 'Estrazione sample dal file .total.cluster.txt'

# ./simpleAnalysis_v3.py "$annotationfile" "$jobid.total.cluster.txt" "$jobid" "$dictionaries" "$pam_file" $mismatch "$referencegenome" "$guide_file" $bulgesDNA $bulgesRNA
./new_simple_analysis.py "$referencegenome" "$dictionaries" "$jobid.total.cluster.txt" "${pam_file}" "$jobid" "$mismatch"
# cp $jobid.bestCFD.txt $jobid.bestCFD.txt.check_analysis
# cp $jobid.bestmmblg.txt $jobid.bestmmblg.txt.check_analysis
# cp $jobid.bestCRISTA.txt $jobid.bestCRISTA.txt.check_analysis

# OUTPUT    $jobid.bestCFD.txt
#           $jobid.CFDGraph.txt     (per fare l'area graph dei CFD REF vs ENR)
# NOTA AnnotatorAllTargets.py salva su disco SOLO il target con CFD più alto nel cluster e tra le scomposizioni esistenti
# Quindi i files jobid.samples.annotation (contentente le scomposizioni del TOP1 esistenti) e jobid.cluster.tmp.txt (contenente il miglior TOP1
# scomposto e i target dei cluster con quell'aplotipo) NON sono creati
rm "$jobid.total.cluster.txt"

echo 'Sorting and adjusting results'
#copy header in tmp file
# head -1 $jobid.bestCFD.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,score,total(mm+bul)
# tail -n +2 $jobid.bestCFD.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k21,21rg -k10,10n -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestCFD.txt
# # cp $jobid.bestCFD.txt $jobid.bestCFD.txt.after_sort

# #copy header in tmp file
# head -1 $jobid.bestmmblg.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,total(mm+bul)
# tail -n +2 $jobid.bestmmblg.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k10,10n -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestmmblg.txt
# # cp $jobid.bestmmblg.txt $jobid.bestmmblg.txt.after_sort

# #copy header in tmp file
# head -1 $jobid.bestCRISTA.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,score,total(mm+bul)
# tail -n +2 $jobid.bestCRISTA.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k21,21rg -k10,10n -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestCRISTA.txt
# cp $jobid.bestCRISTA.txt $jobid.bestCRISTA.txt.after_sort

#adjustin columns to have the correct order and remove uncessary ones
./adjust_cols.py $jobid.bestCFD.txt
./adjust_cols.py $jobid.bestmmblg.txt
./adjust_cols.py $jobid.bestCRISTA.txt
# ./adjust_cols.py $jobid.altmmblg.txt

# sed -i 1i"#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\tDirection\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref\t_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref\tCRISTA_#Bulge_type\tCRISTA_crRNA\tCRISTA_DNA\tCRISTA_Reference\tCRISTA_Chromosome\tCRISTA_Position\tCRISTA_Cluster_Position\tCRISTA_Direction\tCRISTA_Mismatches\tCRISTA_Bulge_Size\tCRISTA_Total\tCRISTA_PAM_gen\tCRISTA_Var_uniq\tCRISTA_Samples\tCRISTA_Annotation_Type\tCRISTA_Real_Guide\tCRISTA_rsID\tCRISTA_AF\tCRISTA_SNP\tCRISTA_#Seq_in_cluster\tCRISTA_CFD\tCRISTA_CFD_ref" "$final_res"
# sed -i 1i"_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref/" $jobid.bestCFD.txt
# sed -i 1i"_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref/" $jobid.bestmmblg.txt
# sed -i 1i"_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref/" $jobid.bestCRISTA.txt

# sed -i '1s/.*/_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref/' $jobid.bestmmblg.txt
# sed -i '1s/.*/_#Bulge_type\t_crRNA\t_DNA\t_Reference\t_Chromosome\t_Position\t_Cluster_Position\t_Direction\t_Mismatches\t_Bulge_Size\t_Total\t_PAM_gen\t_Var_uniq\t_Samples\t_Annotation_Type\t_Real_Guide\t_rsID\t_AF\t_SNP\t_#Seq_in_cluster\t_CFD\t_CFD_ref/' $jobid.altmmblg.txt

# pr -m -t -J $jobid.bestCFD.txt $jobid.bestmmblg.txt >$jobid.bestMerge.txt
# pr -m -t -J $jobid.altCFD.txt $jobid.altmmblg.txt > $jobid.altMerge.txt

# rm $jobid.bestCFD.txt
# rm $jobid.altCFD.txt
# rm $jobid.bestmmblg.txt
# rm $jobid.altmmblg.txt

# mv $jobid.bestMerge.txt $output_folder
# mv $jobid.altMerge.txt $output_folder
# mv $jobid.CFDGraph.txt $output_folder

# echo -e 'Merging Targets\tStart\t'$(date) >>$log
#echo -e 'Merging Close Targets\tStart\t'$(date) >> $log

#merge targets in same chr when they are at distance 3 from each other (inclusive)
# ./merge_close_targets_cfd.sh $jobid.bestCFD.txt $jobid.bestCFD.txt.trimmed 3 'score'
# # cp $jobid.bestCFD.txt.trimmed $jobid.bestCFD.txt.check_merge
# # cp $jobid.bestCFD.txt.trimmed.discarded_samples $jobid.bestCFD.txt.trimmed.discarded_samples.check_merge
# mv $jobid.bestCFD.txt.trimmed $jobid.bestCFD.txt
# mv $jobid.bestCFD.txt.trimmed.discarded_samples $jobid.bestCFD.txt.alt

# #merge targets in same chr when they are at distance 3 from each other (inclusive)
# ./merge_close_targets_cfd.sh $jobid.bestmmblg.txt $jobid.bestmmblg.txt.trimmed 3 'total'
# # cp $jobid.bestmmblg.txt.trimmed $jobid.bestmmblg.txt.check_merge
# # cp $jobid.bestmmblg.txt.trimmed.discarded_samples $jobid.bestmmblg.txt.trimmed.discarded_samples.check_merge
# mv $jobid.bestmmblg.txt.trimmed $jobid.bestmmblg.txt
# mv $jobid.bestmmblg.txt.trimmed.discarded_samples $jobid.bestmmblg.txt.alt

# #merge targets in same chr when they are at distance 3 from each other (inclusive)
# ./merge_close_targets_cfd.sh $jobid.bestCRISTA.txt $jobid.bestCRISTA.txt.trimmed 3 'score'
# # cp $jobid.bestCRISTA.txt.trimmed $jobid.bestCRISTA.txt.check_merge
# # cp $jobid.bestCRISTA.txt.trimmed.discarded_samples $jobid.bestCRISTA.txt.trimmed.discarded_samples.check_merge
# mv $jobid.bestCRISTA.txt.trimmed $jobid.bestCRISTA.txt
# mv $jobid.bestCRISTA.txt.trimmed.discarded_samples $jobid.bestCRISTA.txt.alt

# echo -e 'Annotating results\tStart\t'$(date) >>$log
#annotate bestCFD
# ./annotate_final_results.py $jobid.bestCFD.txt $annotationfile $jobid.bestCFD.txt.annotated
# ./annotate_final_results.py $jobid.bestCFD.txt.alt $annotationfile $jobid.bestCFD.txt.alt.annotated
# mv $jobid.bestCFD.txt.annotated $jobid.bestCFD.txt
# mv $jobid.bestCFD.txt.alt.annotated $jobid.bestCFD.txt.alt
# #annotate bestmmblg
# ./annotate_final_results.py $jobid.bestmmblg.txt $annotationfile $jobid.bestmmblg.txt.annotated
# ./annotate_final_results.py $jobid.bestmmblg.txt.alt $annotationfile $jobid.bestmmblg.txt.alt.annotated
# mv $jobid.bestmmblg.txt.annotated $jobid.bestmmblg.txt
# mv $jobid.bestmmblg.txt.alt.annotated $jobid.bestmmblg.txt.alt
# #annotate bestCRISTA
# ./annotate_final_results.py $jobid.bestCRISTA.txt $annotationfile $jobid.bestCRISTA.txt.annotated
# ./annotate_final_results.py $jobid.bestCRISTA.txt.alt $annotationfile $jobid.bestCRISTA.txt.alt.annotated
# mv $jobid.bestCRISTA.txt.annotated $jobid.bestCRISTA.txt
# mv $jobid.bestCRISTA.txt.alt.annotated $jobid.bestCRISTA.txt.alt
# #correct files names
# #bestCFD

# #bestmmblg

# #bestCRISTA

# # echo -e 'Annotating results\tEnd\t'$(date) >>$log

# # echo -e "Adding risk score"
# #scoring bestCFD
# ./add_risk_score.py $jobid.bestCFD.txt $jobid.bestCFD.txt.risk "False"
# ./add_risk_score.py $jobid.bestCFD.txt.alt $jobid.bestCFD.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestCFD.txt.risk $jobid.bestCFD.txt
# mv $jobid.bestCFD.txt.alt.risk $jobid.bestCFD.txt.alt
# #scoring bestmmblg
# ./add_risk_score.py $jobid.bestmmblg.txt $jobid.bestmmblg.txt.risk "False"
# ./add_risk_score.py $jobid.bestmmblg.txt.alt $jobid.bestmmblg.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestmmblg.txt.risk $jobid.bestmmblg.txt
# mv $jobid.bestmmblg.txt.alt.risk $jobid.bestmmblg.txt.alt
# #scoring bestCRISTA
# ./add_risk_score.py $jobid.bestCRISTA.txt $jobid.bestCRISTA.txt.risk "False"
# ./add_risk_score.py $jobid.bestCRISTA.txt.alt $jobid.bestCRISTA.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestCRISTA.txt.risk $jobid.bestCRISTA.txt
# mv $jobid.bestCRISTA.txt.alt.risk $jobid.bestCRISTA.txt.alt
# # echo -e "Risk score added"

# #remove N's and dots from rsID from best
# ./remove_n_and_dots.py $jobid.bestCFD.txt
# ./remove_n_and_dots.py $jobid.bestmmblg.txt
# ./remove_n_and_dots.py $jobid.bestCRISTA.txt
# #remove N's and dots from rsID from alt
# ./remove_n_and_dots.py $jobid.bestCFD.txt.alt
# ./remove_n_and_dots.py $jobid.bestmmblg.txt.alt
# ./remove_n_and_dots.py $jobid.bestCRISTA.txt.alt
