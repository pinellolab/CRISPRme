#!/bin/bash

set -e  # Exit on error

# This script is designed to test whether CRISPRme correctly retrieves all
# variant-aware (SNP-based) off-target sites for a given guide RNA.
#
# Specifically, the benchmark focuses on guide sg1617 with PAM NGG on the
# human reference genome hg38, enriched with genetic variation from the
# 1000 Genomes Project. The goal is to verify that CRISPRme captures all
# potential off-targets introduced by population-level single-nucleotide
# variants, including those not present in the reference genome alone.
#
# As a ground truth, we use off-target sites obtained via a brute-force
# pairwise sequence alignment approach, which exhaustively compares the
# guide sequence against all variant-enriched genomic sequences. These
# alignments are precomputed to avoid prohibitive runtime costs.
# Implementation details of the brute-force approach can be found at:
# https://github.com/benjaminvyshedskiy/Dynamic_checker
#
# To ensure scalability and practical runtimes, the analysis is restricted
# to a maximum of 4 mismatches and up to 1 DNA and RNA bulges. These
# parameters reflect commonly adopted thresholds in CRISPR off-target
# analyses while still capturing a broad spectrum of plausible off-target
# events.
#
# This benchmark is intended for validation and reproducibility purposes
# and is not meant to serve as a general-purpose off-target discovery
# pipeline.


# URLs and servers
HG38URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38"
VCF1000GSERVER="ftp.1000genomes.ebi.ac.uk"
VCF1000GURL="/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.{}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
TESTDATAURL="https://raw.githubusercontent.com/pinellolab/CRISPRme/refs/heads/main/test/data"

# Directory structure
GENOMES_DIR="Genomes"
VCFS_DIR="VCFs"
ANNOTATION_DIR="Annotations"
PAMS_DIR="PAMs"
SAMPLESIDS_DIR="samplesIDs"

# MD5 checksums
declare -A MD5GENOME=(
    ["hg38.chromFa.tar.gz"]="a5aa5da14ccf3d259c4308f7b2c18cb0"
)

declare -A MD51000G=(
    ["ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="77f154e53c2b7c36b04d03bab3af8b74"
    ["ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="f9d29c4935e591b2b269eed7cd7e35d8"
    ["ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="6e59d00235de71562b4199e09b7e5934"
    ["ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="70a2c1ede97eceb7baeea06c8e46cf3c"
    ["ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="74d5486c0fd29b0e6add24d3740fc3b4"
    ["ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="8c5d83c1a9253058120368af39baf0c8"
    ["ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="dfaa282712fc1292146173dd2ffeb1d9"
    ["ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="ddf7b370fcee63462037c237f12b4444"
    ["ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="5ade69521dc50d88ad7c91bf4ec6fcd8"
    ["ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="1c409a674426eda2fd29b49078137c5d"
    ["ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="65339bffc61bc97f2130832fe9f84d7c"
    ["ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="9a1bda389121140d30c768ef6a1b1370"
    ["ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="47b0463541be137a8bbfe40f6aade864"
    ["ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="241aedf0792c45d5345d421105c782af"
    ["ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="b48e7c64e35b727d34786faa76467f94"
    ["ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="1ce7d66799cab6718852d78dd2aab765"
    ["ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="ecc22783fd1ee7a1c66b053491873192"
    ["ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="fdf3e460e91cd955a9e8cebf01b5d815"
    ["ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="a2f17e4ec552fc07cbd05c1eac0cf7ec"
    ["ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="155c3b440d7990630132e4756f7fcc85"
    ["ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="52882490028507e5d4e606b0905072b1"
    ["ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="57a1722e6ed7d9df08cb3c0e42b62d53"
    ["ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"]="e6a3d41811faee60de177061edcd6fe6"
)

declare -A MD5ANNOTATION=(
    ["dhs+encode+gencode.hg38.bed.tar.gz"]="d3325e347c731b7c24c579a91b447b1b"
    ["gencode.protein_coding.bed.tar.gz"]="c6747bf2610ff144daafc8b02cef251d"
)

declare -A MD5SAMPLES=(
    ["samplesIDs.1000G.txt"]="720af666c9a938de74a2808033aa4509"
    ["samplesIDs.HGDP.txt"]="f92e14e5317221486f20597560ca3a31"
)

# Valid chromosomes
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Function to compute MD5
compute_md5() {
    local file="$1"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        md5 -q "$file"
    else
        md5sum "$file" | cut -d' ' -f1
    fi
}

# Function to verify MD5
verify_md5() {
    local file="$1"
    local expected_md5="$2"
    local actual_md5=$(compute_md5 "$file")
    
    if [ "$actual_md5" != "$expected_md5" ]; then
        echo "ERROR: MD5 mismatch for $(basename $file)" >&2
        echo "Expected: $expected_md5" >&2
        echo "Got: $actual_md5" >&2
        return 1
    fi
    return 0
}

# Function to ensure directory exists
ensure_directory() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
    echo "$dir"
}

# Main workflow
run_crisprme_test() {    
    # Ensure required directories exist
    ensure_directory "$GENOMES_DIR"
    ensure_directory "$VCFS_DIR"
    ensure_directory "$ANNOTATION_DIR"
    ensure_directory "$PAMS_DIR"
    ensure_directory "$SAMPLESIDS_DIR"
    
    # Download genome data
    echo "Downloading reference genome data"
    local archive="${GENOMES_DIR}/hg38.chromFa.tar.gz"
    wget -O "$archive" "$HG38URL/bigZips/hg38.chromFa.tar.gz"
    verify_md5 "$archive" "${MD5GENOME[hg38.chromFa.tar.gz]}"
    tar -xzf "$archive" -C "$GENOMES_DIR"
    mv "${GENOMES_DIR}/chroms" "${GENOMES_DIR}/hg38"
    
    # Download VCF data    
    echo "Downloading 1000 Genomes Project VCF data"    
    local vcf_dataset_dir=$(ensure_directory "$dest/hg38_1000G")
    local server="$VCF1000GSERVER"
    local url_template="$VCF1000GURL"
    local chroms=("${CHROMS[@]}")
    for c in "${chroms[@]}"; do
        local url=$(echo "$url_template" | sed "s/{}/$c/")
        local filename=$(basename "$url")
        local output="$vcf_dataset_dir/$filename"
        wget -O "$output" "ftp://$server$url"
        verify_md5 "$output" "${MD51000G[$filename]}"
    done
    
    # Write VCF config    
    echo "Creating VCF config file for 1000G dataset"
    local vcf_config="vcf.config.test.txt"
    echo "hg38_1000G" >> "$vcf_config"
    vcf=$(echo "$vcf_config")
    
    # Download samples IDs
    echo "Downloading sample ids for 1000G dataset"
    local samplesids_dir=$(ensure_directory "$SAMPLESIDS_DIR")
    local filename="samplesIDs.1000G.txt"
    local output="$samplesids_dir/$filename"
    wget -O "$output" "$TESTDATAURL/samplesIDs/$filename"
    verify_md5 "$output" "${MD5SAMPLES[$filename]}"
    
    # Write samples IDs config
    echo "Creating samples config file for 1000G dataset"
    local samples_config="samplesIDs.config.test.txt"    
    echo "samplesIDs.1000G.txt" >> "$samples_config"
    samplesids=$(echo "$samples_config")
    
    # Download annotation data
    echo "Downloading ENCODE and GENCODE annotation data"
    local annotation_dir=$(ensure_directory "$ANNOTATION_DIR")
    # Download GENCODE
    local gencode_tar="$annotation_dir/gencode.protein_coding.bed.tar.gz"
    wget -O "$gencode_tar" "$TESTDATAURL/Annotations/gencode.protein_coding.bed.tar.gz"
    verify_md5 "$gencode_tar" "${MD5ANNOTATION[gencode.protein_coding.bed.tar.gz]}"
    tar -xzf "$gencode_tar" -C "$annotation_dir"
    bgzip -f "$annotation_dir/gencode.protein_coding.bed"
    gencode=$(echo "$annotation_dir/gencode.protein_coding.bed.gz")
    
    # Download ENCODE
    local encode_tar="$annotation_dir/dhs+encode+gencode.hg38.bed.tar.gz"
    wget -O "$encode_tar" "$TESTDATAURL/Annotations/dhs+encode+gencode.hg38.bed.tar.gz"
    verify_md5 "$encode_tar" "${MD5ANNOTATION[dhs+encode+gencode.hg38.bed.tar.gz]}"
    tar -xzf "$encode_tar" -C "$annotation_dir"
    bgzip -f "$annotation_dir/dhs+encode+gencode.hg38.bed"
    encode=$(echo "$annotation_dir/dhs+encode+gencode.hg38.bed.gz")
    
    # Write PAM file
    echo "Creating PAM file"
    local pams_dir=$(ensure_directory "$PAMS_DIR")
    local pamfile="$pams_dir/20bp-NGG-SpCas9.txt"    
    echo "NNNNNNNNNNNNNNNNNNNNNGG 3" > "$pamfile"
    pam=$(echo "$pamfile")
    
    # Write guide file
    echo "Creating guide file"
    local guidefile="sg1617_test_guide.txt"
    echo "CTAACAGTTGCTTTTATCACNNN" > "$guidefile"
    guide=$(echo "$guidefile")
        
    # Run CRISPRme
    crisprme.py complete-search \
        --genome "$GENOMES_DIR/hg38" \
        --thread 8 \
        --mm 4 \
        --bDNA 1 \
        --bRNA 1 \
        --merge 3 \
        --pam "$pam" \
        --guide "$guide" \
        --vcf "$vcf" \
        --samplesID "$samplesids" \
        --annotation "$encode" \
        --gene_annotation "$gencode" \
        --output crisprme-ci-cd-test \
        --debug \
        --ci-cd-test
}

# Run the test
run_crisprme_test
