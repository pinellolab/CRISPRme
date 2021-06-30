#!/bin/bash

#md5 value for both pacakages
md5_complete_package=b9111f5e1fc9de77ae6f098e12a6ab4e
md5_VCF_package=7feb40213b329c506d6a5e16b7dc10cf

echo "download the complete crisprme package and then unzip it"
wget -c -q -O complete_test_package_no_VCFs.tar.gz https://www.dropbox.com/s/pzfeb1k9v9ekyhr/complete_test_package_no_VCFs.tar.gz?dl=1
wget -c -q -O VCFs.tar.gz https://www.dropbox.com/s/v1fxyhopjek6ib1/VCFs.tar.gz?dl=1
# axel -a -q -n 1 https://www.dropbox.com/s/pzfeb1k9v9ekyhr/complete_test_package_no_VCFs.tar.gz?dl=1
# axel -a -q -n 1 https://www.dropbox.com/s/v1fxyhopjek6ib1/VCFs.tar.gz?dl=1

#check value for complete package
test_md5=$(md5sum complete_test_package_no_VCFs.tar.gz | awk '{print $1}')
if [ "$md5_complete_package" != "$test_md5" ]; then
    echo "error during file download, please rerun this test"
    exit 1
fi

#check value for VCFs package
test_md5=$(md5sum VCFs.tar.gz | awk '{print $1}')
if [ "$md5_VCF_package" != "$test_md5" ]; then
    echo "error during file download, please rerun this test"
    exit 1
fi

echo "unziping all packages"
mkdir crisprme_complete_test
tar -xzf complete_test_package_no_VCFs.tar.gz --directory crisprme_complete_test/
tar -xzf VCFs.tar.gz --directory crisprme_complete_test/
cd crisprme_complete_test/

echo "generate input test data"
echo "CTAACAGTTGCTTTTATCACNNN" >sg1617.txt
echo "hg38_1000G/" >list_vcf.txt
echo "hg38_1000G.samplesID.txt" >list_samplesID.txt
echo "crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NNN-SpCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 8" >crisprme_complete_auto_test.sh

echo "to test the software, launch crisprme_complete_auto_test.sh"
