#crisprme download and test
echo "starting download and unzip of data"

echo "unzip gencode+encode annotations"
#unzip annotations
cd Annotations/
tar -xzf encode+gencode.hg38.tar.gz
rm encode+gencode.hg38.tar.gz
tar -xzf gencode.protein_coding.tar.gz
rm gencode.protein_coding.tar.gz
cd ../

#unzip gencode
# echo "unzip gencode protein-coding proximity file"
# cd Gencode/
# tar -xzf gencode.protein_coding.tar.gz
# rm gencode.protein_coding.tar.gz
# cd ../

echo "start download VCF data and genome (this may take a long time due to connection speed)"
#download VCFs data
cd VCFs/
#download 1000G
cd hg38_1000G/
echo "download 1000G VCFs"
for i in {1..22}; do
   wget -c -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr$i.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
done
wget -c -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
cd ../../

#download HGDP
#uncomment these lines if you want to download also HGDP VCFs
# cd ../hg38_HGDP
# echo "download HGDP VCFs"
# for i in {1..22}
# do
#    wget -c -q ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr$i.vcf.gz
# done
# wget -c -q ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrX.vcf.gz
# wget -c -q ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrY.vcf.gz
# cd ../../

#download hg38
cd Genomes/
echo "download hg38"
wget -c -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xzf hg38.chromFa.tar.gz
mv chroms hg38
cd ../

echo "start testing"
docker run -v ${PWD}:/DATA -w /DATA -i scancellieri/crisprme crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 4
