#!/bin/sh
#$1 è cwd

genome_general_dir=$1'/Genomes/'

#Creo file per indicare che è in corso un'azione

touch $genome_general_dir'aggiunta_nuovo_genoma.txt'
out_file=$genome_general_dir'aggiunta_nuovo_genoma.txt'

#Comandi per eseguire le azioni
#Esempio
echo '1 Copy_Genome' > $out_file
sleep 5

echo '2 Indexing' > $out_file
sleep 10

echo '3 Enrichment' > $out_file
sleep 10

#....

echo '4 Done' > $out_file

