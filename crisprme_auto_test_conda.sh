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
   sleep 2  # allow disk sync for md5sum check
   local_md5sum="$(md5sum ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | cut -d ' ' -f 1)"
   if [ "$original_md5sum" != "$local_md5sum" ]; then  # check download consistency
      echo "ERROR: unexpected failure while downloading ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
      exit 1
   fi
done
cd ../..

# initialize VCF config file
VCFCONFIG="vcf_config.1000G.txt"
printf "${VCF1000G}\n" > $VCFCONFIG

# download 1000G samplesIDs  
SAMPLESIDS="samplesIDs"
mkdir -p $SAMPLESIDS  # create sample ids dir
cd $SAMPLESIDS
# download 1000G samples IDs
echo "Downloading samples ids for 1000G dataset"
SAMPLES1000G="hg38_1000G.samplesID.txt"
original_md5sum="$(curl -sL https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/gnomad-4.1-converter/download_data/${SAMPLES1000G} | md5sum | cut -d ' ' -f 1)"  # compute original md5sum
while true; do # retry download if caught timeout
   wget -T 15 -c https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/gnomad-4.1-converter/download_data/${SAMPLES1000G} && break
done
sleep 2  # allow disk sync for md5sum check
local_md5sum="$(md5sum $SAMPLES1000G | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading ${SAMPLES1000G}"
   exit 1
fi
cd ..

# initialize samples config file
SAMPLESCONFIG="samplesIDs.1000G.txt"
printf "${SAMPLES1000G}\n" > $SAMPLESCONFIG

# download annotation data
ANNOTATIONDIR="Annotations"
mkdir -p $ANNOTATIONDIR  # create annotation folder
cd $ANNOTATIONDIR
echo "Downloading ENCODE+GENCODE annotation data..."
original_md5sum="$(curl -sL https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data/dhs+encode+gencode.hg38.bed.tar.gz | md5sum | cut -d ' ' -f 1)"
DHSENCODE="dhs+encode+gencode.hg38.bed.zip"
while true; do  # retry download if caught timeout
   wget -T 15 -c -O $DHSENCODE https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data/dhs+encode+gencode.hg38.bed.tar.gz && break
done
sleep 2  # allow disk sync for md5sum check
local_md5sum="$(md5sum $DHSENCODE | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading ${DHSENCODE}"
   exit 1
fi
echo "Extracting ${DHSENCODE}..."
tar -xvf $DHSENCODE
DHSENCODE="dhs+encode+gencode.hg38.bed"

echo "Downloading GENCODE encoding sequences..."
original_md5sum="$(curl -sL https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data/gencode.protein_coding.bed.tar.gz | md5sum | cut -d ' ' -f 1)"
GENCODE="gencode.protein_coding.bed.zip"
while true; do  # retry download if caught timeout
   wget -T 15 -c -O $GENCODE https://raw.githubusercontent.com/pinellolab/CRISPRme/gnomad-4.1-converter/download_data/gencode.protein_coding.bed.tar.gz && break
done
sleep 2  # allow disk sync for md5sum check
local_md5sum="$(md5sum $GENCODE | cut -d ' ' -f 1)"
if [ "$original_md5sum" != "$local_md5sum" ]; then
   echo "ERROR: unexpected failure while downloading ${GENCODE}"
   exit 1
fi
echo "Extracting ${GENCODE}..."
tar -xvf $GENCODE
GENCODE="gencode.protein_coding.bed"
cd ..

# create Dictionaries folder
mkdir -p "Dictionaries"

# create sg1617 guide file
GUIDEFILE="sg1617.txt"
printf "CTAACAGTTGCTTTTATCACNNN\n" > $GUIDEFILE

# create NGG PAM file
PAM="PAMs"
mkdir -p $PAM
cd $PAM
NGGPAM="20bp-NGG-spCas9.txt"
printf "NNNNNNNNNNNNNNNNNNNNNGG 3\n" > $NGGPAM
cd ..

echo "Start CRISPRme test..."
crisprme.py complete-search --genome Genomes/$HG38 --vcf $VCFCONFIG --guide $GUIDEFILE --pam PAMs/$NGGPAM --annotation Annotations/$DHSENCODE --samplesID $SAMPLESCONFIG --gene_annotation Annotations/$GENCODE --sorting-criteria-scoring mm+bulges --sorting-criteria mm,bulges --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 4
