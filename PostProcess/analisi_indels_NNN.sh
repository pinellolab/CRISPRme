#!/bin/bash

set -e  # trace all errors 

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

#awk '{print $0"\tn"}' "$REFtargets" > "$REFtargets.corrected"
#echo 'Adjusting indel positions'
#./correct_positions_targets.py "$ENRtargets"
#./bed_for_ref_seq.py "$ENRtargets.corrected" "$jobid"
#bedtools getfasta -fi "$referencegenome" -bed "$jobid.bed_for_ref" -fo "$jobid.ref_seq.fa"
touch $REFtargets.corrected

# 1) Rimozione duplicati, estrazione semicommon e unique e creazione file total
#echo 'Creazione file .total.txt'
./extraction.sh "$REFtargets.corrected" "$ENRtargets" "$jobid" || {
	echo "CRISPRme ERROR: indels analysis failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
} # OUTPUT    $jobid.common_targets.txt -> Non usato
#           $jobid.semi_common_targets.txt
#           $jobid.unique_targets.txt

rm "$jobid.common_targets.txt"
rm "$REFtargets.corrected"
#rm "$ENRtargets.corrected"

# 2) Creazione colonne PAM creation etc
awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' "$jobid.semi_common_targets.txt" >"$jobid.semi_common_targets.minmaxdisr.txt"
rm "$jobid.semi_common_targets.txt"

awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' "$jobid.unique_targets.txt" >"$jobid.unique_targets.pamcreation.txt" #Add pam creation, variant unique, real guide column
rm "$jobid.unique_targets.txt"

cat "$jobid.unique_targets.pamcreation.txt" "$jobid.semi_common_targets.minmaxdisr.txt" >"$jobid.total.txt"
rm "$jobid.unique_targets.pamcreation.txt"
rm "$jobid.semi_common_targets.minmaxdisr.txt"

#awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide"\tn\tn\tn"}' $ENRtargets.corrected > $jobid.total.txt

#echo 'Creazione cluster del file .total.txt'
# 3) Clustering
./cluster.dict.py "$jobid.total.txt" 'no' 'True' 'True' "$guide_file" 'total' 'orderChr' || {
	echo "CRISPRme ERROR: indels clustering failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
} # OUTPUT     $jobid.total.cluster.txt

#sed -i ':a;N;$!ba;s/\n/\tn\tn\tn\n/g' $jobid.total.cluster.txt
#sed -i '$s/$/\tn\tn\tn/g' $jobid.total.cluster.txt

sed -i '1s/.*/#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tChromosome_fake\tPAM_gen\tVar_uniq\tSamples\tAnnotation Type\tReal Guide\tID\tAF\tSNP/' "$jobid.total.cluster.txt"
rm "$jobid.total.txt"

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
# OUTPUT    chr1-5.bestCFD_INDEL.txt      chr6-15.bestCFD_INDEL.txt     chr16-22XYM.bestCFD_INDEL.txt
#           chr1-5.CFDGraph.txt     chr6-15.CFDGraph.txt    chr16-22XYM.CFDGraph.txt
# I file .bestCFD_INDEL.txt si possono unire (attenzione che avrò 3 header all'interno del file) per ottenere il file completo
# I file .CFDGraph.txt vanno sommati tra loro per ottenere il file finale -> AL MOMENTO NON NECESSARIO PER QUESTA ANALISI (TENERE COMUNQUE I FILE)

#echo 'Estrazione sample dal file .total.cluster.txt'

./analisi_indels_NNN.py "$annotationfile" "$jobid.total.cluster.txt" "$jobid" "$dictionaries" "$pam_file" "$mismatch" "$referencegenome" "$guide_file" $bulgesDNA $bulgesRNA || {
	echo "CRISPRme ERROR: indels analysis failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
# OUTPUT    $jobid.bestCFD_INDEL.txt
#           $jobid.CFDGraph.txt     (per fare l'area graph dei CFD REF vs ENR)
# NOTA AnnotatorAllTargets.py salva su disco SOLO il target con CFD più alto nel cluster e tra le scomposizioni esistenti
# Quindi i files jobid.samples.annotation (contentente le scomposizioni del TOP1 esistenti) e jobid.cluster.tmp.txt (contenente il miglior TOP1
# scomposto e i target dei cluster con quell'aplotipo) NON sono creati
rm "$jobid.total.cluster.txt"
#rm "$jobid.bed_for_ref"
#rm "$jobid.ref_seq.fa"

echo 'Sorting and adjusting results'
#copy header in tmp file
# head -1 $jobid.bestCFD_INDEL.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,score
# tail -n +2 $jobid.bestCFD_INDEL.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k21,21rg -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestCFD_INDEL.txt

# #copy header in tmp file
# head -1 $jobid.bestmmblg_INDEL.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,total(mm+bul)
# tail -n +2 $jobid.bestmmblg_INDEL.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k10,10n -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestmmblg_INDEL.txt

# #copy header in tmp file
# head -1 $jobid.bestCRISTA_INDEL.txt >$jobid.tmp
# #tail file w/o header and sort for realguide,chr,cluster_pos,score
# tail -n +2 $jobid.bestCRISTA_INDEL.txt | LC_ALL=C sort -k15,15 -k4,4 -k6,6n -k21,21rg -T ./ >>$jobid.tmp && mv $jobid.tmp $jobid.bestCRISTA_INDEL.txt

./adjust_cols.py "$jobid.bestCFD_INDEL.txt" || {
	echo "CRISPRme ERROR: CFD indels report failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
./adjust_cols.py "$jobid.bestCRISTA_INDEL.txt" || {
	echo "CRISPRme ERROR: CRISTA indels report failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
./adjust_cols.py "$jobid.bestmmblg_INDEL.txt" || {
	echo "CRISPRme ERROR: mismatch+bulges indels report failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}

# sed -i '1s/.*/MMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' $jobid.bestmmblg_INDEL.txt
# sed -i '1s/.*/MMBLG_#Bulge_type\tMMBLG_crRNA\tMMBLG_DNA\tMMBLG_Reference\tMMBLG_Chromosome\tMMBLG_Position\tMMBLG_Cluster_Position\tMMBLG_Direction\tMMBLG_Mismatches\tMMBLG_Bulge_Size\tMMBLG_Total\tMMBLG_PAM_gen\tMMBLG_Var_uniq\tMMBLG_Samples\tMMBLG_Annotation_Type\tMMBLG_Real_Guide\tMMBLG_rsID\tMMBLG_AF\tMMBLG_SNP\tMMBLG_#Seq_in_cluster\tMMBLG_CFD\tMMBLG_CFD_ref/' $jobid.altmmblg.txt

# pr -m -t -J $jobid.bestCFD_INDEL.txt $jobid.bestmmblg_INDEL.txt >$jobid.bestMerge.txt
# pr -m -t -J $jobid.altCFD.txt $jobid.altmmblg.txt >$jobid.altMerge.txt

./remove_bad_indel_targets.py "$jobid.bestCFD_INDEL.txt" || {
	echo "CRISPRme ERROR: CFD indels report cleaning failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
./remove_bad_indel_targets.py "$jobid.bestCRISTA_INDEL.txt" || {
	echo "CRISPRme ERROR: CRISTA indels report cleaning failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}
./remove_bad_indel_targets.py "$jobid.bestmmblg_INDEL.txt" || {
	echo "CRISPRme ERROR: mismatch+bulges indels report cleaning failed (script: ${0} line $((LINENO-1)))" >&2
	exit 1
}

#merge targets in same chr when they are at distance 3 from each other (inclusive) preserving the highest scoring one
# ./merge_close_targets_cfd.sh $jobid.bestCFD_INDEL.txt $jobid.bestCFD_INDEL.txt.trimmed 3 'score'
# mv $jobid.bestCFD_INDEL.txt.trimmed $jobid.bestCFD_INDEL.txt
# mv $jobid.bestCFD_INDEL.txt.trimmed.discarded_samples $jobid.bestCFD_INDEL.txt.alt

# #merge targets in same chr when they are at distance 3 from each other (inclusive) preserving the lowest total(mm+bul) one
# ./merge_close_targets_cfd.sh $jobid.bestmmblg_INDEL.txt $jobid.bestmmblg_INDEL.txt.trimmed 3 'total'
# mv $jobid.bestmmblg_INDEL.txt.trimmed $jobid.bestmmblg_INDEL.txt
# mv $jobid.bestmmblg_INDEL.txt.trimmed.discarded_samples $jobid.bestmmblg_INDEL.txt.alt

# #merge targets in same chr when they are at distance 3 from each other (inclusive) preserving the highest scoring one
# ./merge_close_targets_cfd.sh $jobid.bestCRISTA_INDEL.txt $jobid.bestCRISTA_INDEL.txt.trimmed 3 'score'
# mv $jobid.bestCRISTA_INDEL.txt.trimmed $jobid.bestCRISTA_INDEL.txt
# mv $jobid.bestCRISTA_INDEL.txt.trimmed.discarded_samples $jobid.bestCRISTA_INDEL.txt.alt

# echo -e 'Annotating results\tStart\t'$(date) >>$log
#annotate bestCFD_INDEL
# ./annotate_final_results.py $jobid.bestCFD_INDEL.txt $annotationfile $jobid.bestCFD_INDEL.txt.annotated
# ./annotate_final_results.py $jobid.bestCFD_INDEL.txt.alt $annotationfile $jobid.bestCFD_INDEL.txt.alt.annotated
# mv $jobid.bestCFD_INDEL.txt.annotated $jobid.bestCFD_INDEL.txt
# mv $jobid.bestCFD_INDEL.txt.alt.annotated $jobid.bestCFD_INDEL.txt.alt
# #annotate bestmmblg_INDEL
# ./annotate_final_results.py $jobid.bestmmblg_INDEL.txt $annotationfile $jobid.bestmmblg_INDEL.txt.annotated
# ./annotate_final_results.py $jobid.bestmmblg_INDEL.txt.alt $annotationfile $jobid.bestmmblg_INDEL.txt.alt.annotated
# mv $jobid.bestmmblg_INDEL.txt.annotated $jobid.bestmmblg_INDEL.txt
# mv $jobid.bestmmblg_INDEL.txt.alt.annotated $jobid.bestmmblg_INDEL.txt.alt
# #annotate bestCRISTA_INDEL
# ./annotate_final_results.py $jobid.bestCRISTA_INDEL.txt $annotationfile $jobid.bestCRISTA_INDEL.txt.annotated
# ./annotate_final_results.py $jobid.bestCRISTA_INDEL.txt.alt $annotationfile $jobid.bestCRISTA_INDEL.txt.alt.annotated
# mv $jobid.bestCRISTA_INDEL.txt.annotated $jobid.bestCRISTA_INDEL.txt
# mv $jobid.bestCRISTA_INDEL.txt.alt.annotated $jobid.bestCRISTA_INDEL.txt.alt
# #correct files names
#bestCFD_INDEL

#bestmmblg_INDEL

#bestCRISTA_INDEL

# echo -e 'Annotating results\tEnd\t'$(date) >>$log

# echo -e "Adding risk score"
#scoring bestCFD_INDEL
# ./add_risk_score.py $jobid.bestCFD_INDEL.txt $jobid.bestCFD_INDEL.txt.risk "False"
# ./add_risk_score.py $jobid.bestCFD_INDEL.txt.alt $jobid.bestCFD_INDEL.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestCFD_INDEL.txt.risk $jobid.bestCFD_INDEL.txt
# mv $jobid.bestCFD_INDEL.txt.alt.risk $jobid.bestCFD_INDEL.txt.alt
# #scoring bestmmblg_INDEL
# ./add_risk_score.py $jobid.bestmmblg_INDEL.txt $jobid.bestmmblg_INDEL.txt.risk "False"
# ./add_risk_score.py $jobid.bestmmblg_INDEL.txt.alt $jobid.bestmmblg_INDEL.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestmmblg_INDEL.txt.risk $jobid.bestmmblg_INDEL.txt
# mv $jobid.bestmmblg_INDEL.txt.alt.risk $jobid.bestmmblg_INDEL.txt.alt
# #scoring bestCRISTA_INDEL
# ./add_risk_score.py $jobid.bestCRISTA_INDEL.txt $jobid.bestCRISTA_INDEL.txt.risk "False"
# ./add_risk_score.py $jobid.bestCRISTA_INDEL.txt.alt $jobid.bestCRISTA_INDEL.txt.alt.risk "False" #"True" change to True if ID_CLUSTER is inserted during merge_phase
# mv $jobid.bestCRISTA_INDEL.txt.risk $jobid.bestCRISTA_INDEL.txt
# mv $jobid.bestCRISTA_INDEL.txt.alt.risk $jobid.bestCRISTA_INDEL.txt.alt
# # echo -e "Risk score added"

# #remove N's and dots from rsID from best
# ./remove_n_and_dots.py $jobid.bestCFD_INDEL.txt
# ./remove_n_and_dots.py $jobid.bestmmblg_INDEL.txt
# ./remove_n_and_dots.py $jobid.bestCRISTA_INDEL.txt
# #remove N's and dots from rsID from alt
# ./remove_n_and_dots.py $jobid.bestCFD_INDEL.txt.alt
# ./remove_n_and_dots.py $jobid.bestmmblg_INDEL.txt.alt
# ./remove_n_and_dots.py $jobid.bestCRISTA_INDEL.txt.alt
