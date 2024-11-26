"""Define static variables and utilities functions used throughout CRISPRme.
"""

import sys
import os


def download_genome(chr: str, directory: str) -> None:
    if not isinstance(chr, str):
        raise TypeError(f"Expected {str.__name__}, got {type(chr).__name__}")
    if not isinstance(directory, str):
        raise TypeError(f"Expected {str.__name__}, got {type(directory).__name__}")
    if chr == "all":
        os.system(
            "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz -O hg38.chromFa.tar.gz"
        )
        os.system("mv hg38.chromFa.tar.gz Genomes/")
        os.system("tar -xzf Genomes/hg38.chromFa.tar.gz")
        os.system(f"mv Genomes/chroms mv Genomes/{directory}")
    else:
        os.system(
            f"wget https://hgdownload2.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chr}.fa.gz -O {chr}.fa.gz"
        )
        os.system(f"mv {chr}.fa.gz Genomes/")
        os.system(f"gunzip Genomes/{chr}.fa.gz")
        os.makedirs(f"Genomes/{directory}", exist_ok=True)
        os.system(f"mv Genomes/{chr}.fa Genomes/{directory}/{chr}.fa")


def download_vcf(chr: str, origin: str) -> None:
    if not isinstance(chr, str):
        raise TypeError(f"Expected {str.__name__}, got {type(chr).__name__}")
    if not isinstance(origin, str):
        raise TypeError(f"Expected {str.__name__}, got {type(origin).__name__}")
    if chr == "all" and origin == "1000G":
        os.makedirs(f"VCFs/hg38_1000G", exist_ok=True)
        for i in range(1, 23):
            os.system(
                f"wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
            )
            os.system(
                f"mv ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz VCFs/hg38_1000G/"
            )
    elif chr != "all" and origin == "1000G":
        os.makedirs(f"VCFs/hg38_1000G", exist_ok=True)
        os.system(
            f"wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        )
        os.system(
            f"mv ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz VCFs/hg38_1000G/"
        )
    elif chr == "all" and origin == "HGDP":
        os.makedirs(f"VCFs/hg38_HGDP", exist_ok=True)
        for i in range(1, 23):
            os.system(
                f"wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr{i}.vcf.gz -O hgdp_wgs.20190516.full.chr{i}.vcf.gz"
            )
            os.system(f"mv hgdp_wgs.20190516.full.chr{i}.vcf.gz VCFs/hg38_HGDP/")
    elif chr != "all" and origin == "HGDP":
        os.makedirs(f"VCFs/hg38_HGDP", exist_ok=True)
        os.system(
            f"wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{chr}.vcf.gz -O hgdp_wgs.20190516.full.{chr}.vcf.gz"
        )
        os.system(f"mv hgdp_wgs.20190516.full.{chr}.vcf.gz VCFs/hg38_HGDP/")


def download_samplesID() -> None:
    os.system(
        "wget https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_1000G.samplesID.txt -O hg38_1000G.samplesID.txt"
    )
    os.system(
        "wget https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_HGDP.samplesID.txt -O hg38_HGDP.samplesID.txt"
    )
    os.system(
        "wget https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_gnomAD.samplesID.txt -O hg38_gnomAD.samplesID.txt"
    )
    os.system("mv hg38_1000G.samplesID.txt samplesIDs/")
    os.system("mv hg38_HGDP.samplesID.txt samplesIDs/")
    os.system("mv hg38_gnomAD.samplesID.txt samplesIDs/")


def download_annotation() -> None:
    os.system(
        "wget https://github.com/pinellolab/CRISPRme/raw/test-function/download_data/gencode.protein_coding.tar.gz -O gencode.protein_coding.tar.gz"
    )
    os.system("tar -xzf gencode.protein_coding.tar.gz")
    os.system("rm gencode.protein_coding.tar.gz")
    os.system("mv gencode.protein_coding.bed Annotations/")
    os.system(
        "wget https://github.com/pinellolab/CRISPRme/raw/test-function/download_data/encode+gencode.hg38.tar.gz -O encode+gencode.hg38.tar.gz"
    )
    os.system("tar -xzf encode+gencode.hg38.tar.gz")
    os.system("rm encode+gencode.hg38.tar.gz")
    os.system("mv encode+gencode.hg38.bed Annotations/")
