# ------------------------------------------------------------------------------
# TEST CRISPRme conda package

echo "Download and extract fundamental data..."

# download hg38 genome assembly
GENOMEDIR="Genomes"
HG38="hg38"
echo "Downloading hg38 genome assembly..."
mkdir -p $GENOMEDIR  # create Genomes folder
cd $GENOMEDIR
# download chromosomes FASTA files
original_md5sum="a5aa5da14ccf3d259c4308f7b2c18cb0"  # see https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt
while true; do  # retry download if caught timeout
   wget -T 15 -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz && break
done
chromsfasta="hg38.chromFa.tar.gz"
local_md5sum="$(md5sum $chromsfasta | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz"
   exit 1
fi
echo "Extracting ${chromsfasta}..."
tar -xzf $chromsfasta
mv chroms $HG38  # create hg38 directory
cd ..

# download 1000G VCFs
echo "Downloading 1000 Genomes VCFs..."
VCF="VCFs"
mkdir -p $VCF  # create VCF folder
cd $VCF
VCF1000G="hg38_1000G"
mkdir -p $VCF1000G  # create 1000 genomes folder
cd $VCF1000G
# download 1000 genomes VCF files
for i in $(seq 1 22; echo "X");
do 
   original_md5sum="$(curl -sL ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | md5sum | cut -d ' ' -f 1)"  # compute original md5sum
   while true; do  # retry download if caught timeout
      wget -T 15 -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz && break
   done
   local_md5sum="$(md5sum ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | cut -d ' ' -f 1)"
   if [ "$original_md5sum" != "$local_md5sum" ]; then  # check download consistency
      echo "ERROR: unexpected failure while downloading ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
      exit 1
   fi
done
cd ../..

# download annotation data
ANNOTATIONDIR="Annotations"
mkdir -p $ANNOTATIONDIR  # create annotation folder
cd $ANNOTATIONDIR
echo "Downloading ENCODE+GENCODE annotation data..."
original_md5sum="$(curl -sL https://www.dropbox.com/s/1n2f0qxdba7u3gb/encode%2Bgencode.hg38.bed.zip?dl=0 | md5sum | cut -d ' ' -f 1)"
encodegencode="encode+gencode.hg38.bed.zip"
while true; do  # retry download if caught timeout
   wget -T 15 -c -O $encodegencode https://www.dropbox.com/s/1n2f0qxdba7u3gb/encode%2Bgencode.hg38.bed.zip?dl=1 && break
done
local_md5sum="$(md5sum $encodegencode | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading ${encodegencode}"
   exit 1
fi
echo "Extracting ${encodegencode}..."
unzip $encodegencode
echo "Downloading GENCODE encoding sequences..."
original_md5sum="$(curl -sL https://www.dropbox.com/s/isqpkg113cr1xea/gencode.protein_coding.bed.zip?dl=0 | md5sum | cut -d ' ' -f 1)"
gencode="gencode.protein_coding.bed.zip"
while true; do  # retry download if caught timeout
   wget -T 15 -c -O $gencode https://www.dropbox.com/s/isqpkg113cr1xea/gencode.protein_coding.bed.zip?dl=1 && break
done
local_md5sum="$(md5sum $gencode | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading ${gencode}"
   exit 1
fi
echo "Extracting ${gencode}..."
unzip $gencode
cd ..

echo "Start CRISPRme test..."
crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 4
