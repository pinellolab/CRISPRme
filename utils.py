"""Define static variables and utilities functions used throughout CRISPRme.
"""

import sys
import os

CRISPRME_DIRS = [
    "Genomes",
    "Results",
    "Dictionaries",
    "VCFs",
    "Annotations",
    "PAMs",
    "samplesIDs",
]


def check_directories(basedir: str) -> None:
    """The function checks the consistency of CRISPRme's directory tree.
    If a directory is not found in the tree, it will be created.

    ...

    Parameters
    ----------
    basedir : str
        Base directory

    Returns
    -------
    None
    """

    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


def download_genome(chr: str, directory: str) -> None:
    if not isinstance(chr, str):
        raise TypeError(f"Expected {str.__name__}, got {type(chr).__name__}")
    if not isinstance(directory, str):
        raise TypeError(f"Expected {str.__name__}, got {type(directory).__name__}")
    if chr == "all":
        os.system(
            "curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz -o hg38.chromFa.tar.gz"
        )
        os.system("mv hg38.chromFa.tar.gz Genomes/")
        os.system("tar -xzf Genomes/hg38.chromFa.tar.gz")
        os.system(f"mv Genomes/chroms mv Genomes/{directory}")
    else:
        os.system(
            f"curl https://hgdownload2.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chr}.fa.gz -o {chr}.fa.gz"
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
                f"curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -o ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
            )
            os.system(
                f"mv ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz VCFs/hg38_1000G/"
            )
    elif chr != "all" and origin == "1000G":
        os.makedirs(f"VCFs/hg38_1000G", exist_ok=True)
        os.system(
            f"curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -o ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        )
        os.system(
            f"mv ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz VCFs/hg38_1000G/"
        )
    elif chr == "all" and origin == "HGDP":
        os.makedirs(f"VCFs/hg38_HGDP", exist_ok=True)
        for i in range(1, 23):
            os.system(
                f"curl ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr{i}.vcf.gz -o hgdp_wgs.20190516.full.chr{i}.vcf.gz"
            )
            os.system(f"mv hgdp_wgs.20190516.full.chr{i}.vcf.gz VCFs/hg38_HGDP/")
    elif chr != "all" and origin == "HGDP":
        os.makedirs(f"VCFs/hg38_HGDP", exist_ok=True)
        os.system(
            f"curl ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.{chr}.vcf.gz -o hgdp_wgs.20190516.full.{chr}.vcf.gz"
        )
        os.system(f"mv hgdp_wgs.20190516.full.{chr}.vcf.gz VCFs/hg38_HGDP/")


def download_samplesID() -> None:
    os.system(
        "curl https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_1000G.samplesID.txt -o hg38_1000G.samplesID.txt"
    )
    os.system(
        "curl https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_HGDP.samplesID.txt -o hg38_HGDP.samplesID.txt"
    )
    os.system(
        "curl https://raw.githubusercontent.com/pinellolab/CRISPRme/test-function/download_data/hg38_gnomAD.samplesID.txt -o hg38_gnomAD.samplesID.txt"
    )
    os.system("mv hg38_1000G.samplesID.txt samplesIDs/")
    os.system("mv hg38_HGDP.samplesID.txt samplesIDs/")
    os.system("mv hg38_gnomAD.samplesID.txt samplesIDs/")


def download_annotation() -> None:
    os.system(
        "curl https://github.com/pinellolab/CRISPRme/raw/test-function/download_data/gencode.protein_coding.tar.gz -o gencode.protein_coding.tar.gz"
    )
    os.system("tar -xzf gencode.protein_coding.tar.gz")
    os.system("mv gencode.protein_coding.bed Annotations/")
    os.system(
        "curl https://github.com/pinellolab/CRISPRme/raw/test-function/download_data/encode+gencode.hg38.tar.gz -o encode+gencode.hg38.tar.gz"
    )
    os.system("tar -xzf encode+gencode.hg38.tar.gz")
    os.system("mv encode+gencode.hg38.bed Annotations/")
