echo "download the complete crisprme package and then unzip it"
gdown "download link"
mkdir crisprme_test
tar -xvf complete_test_package.tar.gz --directory crisprme_test/
cd crisprme_test/

echo "generate input test data"
echo "CTAACAGTTGCTTTTATCACNNN" > sg1617.txt
echo "hg38_1000G/" > list_vcf.txt
echo "hg38_1000G.samplesID.txt" > list_samplesID.txt

echo "to start a test you can copy the following line into your terminal"
echo "crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NNN-SpCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 8"
echo "if you want to use the web-interface, just use this command"
echo "crisprme.py web-interface"
