echo "unzip gencode+encode annotations"
#unzip annotations
cd Annotations/
tar -xzf gencode_encode.hg38.tar.gz &> error.txt
rm -f gencode_encode.hg38.tar.gz error.txt
cd ../

#unzip gencode
#echo "unzip gencode protein-coding proximity file"
#cd Gencode/
#tar -xzf gencode.protein_coding.tar.gz &> error.txt
#rm -f gencode.protein_coding.tar.gz error.txt
#cd ../

echo "start testing"
crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 20
