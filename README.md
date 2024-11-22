# CRISPRme

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crisprme/README.html)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/pinellolab/crisprme)
![Conda](https://img.shields.io/conda/dn/bioconda/crisprme)
![license](https://img.shields.io/badge/license-AGPL--3.0-lightgrey)

CRISPRme is a comprehensive tool designed for thorough off-target assessment in 
CRISPR-Cas systems. Available as a web application 
([http://crisprme.di.univr.it/](http://crisprme.di.univr.it/)), offline tool, and 
command-line interface, it integrates human genetic variant datasets with orthogonal 
genomic annotations to predict and prioritize potential off-target sites at scale. 
CRISPRme accounts for single-nucleotide variants (SNVs) and indels, considers 
*bona fide* haplotypes, and allows for spacer:protospacer mismatches and bulges, 
making it well-suited for both population-wide and personal genome analyses. CRISPRme 
automates the entire workflow, from data download to executing the search, and delivers 
detailed reports complete with tables and figures through an interactive web-based 
interface.

## Table Of Contents

1 [Installation](#1-installation)
<br>&nbsp;&nbsp;1.1 [Install CRISPRme via Conda/Mamba](#11-install-crisprme-via-condamamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.1 [Installing Conda or Mamba](#111-installing-conda-or-mamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.2 [Installing CRISPRme](#112-installing-crisprme)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.3 [Updating CRISPRme](#113-updating-crisprme)
<br>&nbsp;&nbsp;1.2 [Install CRISPRme via Docker](#12-install-crisprme-via-docker)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.2.1 [Installing Docker](#121-installing-docker)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.2.2 [Building and Pulling CRISPRme Docker Image](#122-building-and-pulling-crisprme-docker-image)

## 1 Installation

This section outlines the steps to install CRISPRme, tailored to suit different 
operating systems. Select the method that best matches your setup:

- [Install CRISPRme via Conda/Mamba (for Linux users)](#11-install-crisprme-via-condamamba)
  
- [Install CRISPRme via Docker (compatible with all operating systems)]()

Each method ensures a streamlined and efficient installation, enabling you to use 
CRISPRme with minimal effort. Follow the detailed instructions provided in the 
respective sections below.

### 1.1 Install CRISPRme via Conda/Mamba
---

This section is organized into three subsections to guide you through the installation 
and maintenance of CRISPRme:

- [Installing Conda or Mamba](#111-installing-conda-or-mamba):
  <br>This subsection provides step-by-step instructions to install either 
  Conda or Mamba. Begin here if you do not have these package managers installed 
  on your machine.

- [Installing CRISPRme](#112-installing-crisprme):
  <br>Once you have Conda or Mamba set up, proceed to this subsection for detailed 
  instructions on creating the CRISPRme environment and installing the necessary 
  dependencies.

- [Updating CRISPRme](#113-updating-crisprme):
  <br>Learn how to update an existing CRISPRme installation to the latest version, 
  ensuring access to new features and bug fixes.

#### 1.1.1 Installing Conda or Mamba
---

Before installing CRISPRme, ensure you have either Conda or Mamba installed on 
your machine. Based on recommendations from the Bioconda community, we highly 
recommend using Mamba over Conda. Mamba is a faster, more efficient drop-in replacement 
for Conda, leveraging a high-performance dependency solver and components 
optimized in C++.

**Step1: Install `Conda` or `Mamba`**

- To install `Conda`, refer to the official installation guide:
  <br>[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- To install `Mamba`, refer to the official installation guide: 
  <br>[Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

**Step 2: Configure Bioconda Channels**

Once `Mamba` is installed, configure it to use Bioconda and related channels by 
running the following one-time setup commands:
```bash
mamba config --add channels bioconda
mamba config --add channels defaults
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```
> **Note:** If you prefer to use `Conda`, replace `mamba` with `conda` in the 
commands above

By completing these steps, your system will be fully prepared for installing CRISPRme.

#### 1.1.2 Installing CRISPRme
---

We strongly recommend using **`Mamba`** to create CRISPRme's `conda` environment 
due to its superior speed and reliability in dependency management. However, if 
you prefer `Conda`, you can replace `mamba` with `conda` in all the commands below.

**Step 1: Create CRISPRme's Environment**  

Open a terminal and execute the following command:  

```bash
mamba create -n crisprme python=3.9 crisprme -y  # Install CRISPRme and its dependencies
```  

This command sets up a dedicated `conda` environment named `crisprme`, installing 
CRISPRme along with all required dependencies.

**Step 2: Activate the Environment**  

To activate the newly created CRISPRme environment, type:  

```bash
mamba activate crisprme  # Enable the CRISPRme environment
```  

**Step 3: Test the Installation**  

To verify that CRISPRme is correctly installed, run the following commands in your terminal:  

```bash
crisprme.py --version  # Display the installed CRISPRme version
crisprme.py            # List CRISPRme functionalities
```  

- The first command will output the version of CRISPRme (e.g., `2.1.6`).  
- The second command should display CRISPRme's functionalities.  

If both commands execute successfully, your installation is complete, and 
CRISPRme is ready to use.

#### 1.1.3 Updating CRISPRme
---

To update an existing CRISPRme installation using `Mamba` or `Conda`, follow the 
steps below:

**Step 1: Check the Latest Version**

Visit the CRISPRme README to identify the latest version of the tool.

**Step 2: Update CRISPRme**

Run the following command in your terminal, replacing <latest_version> with the 
desired version number:
```bash
mamba install crisprme=<latest_version>  # Update CRISPRme to the specified version
```

For example, to update CRISPRme to version `2.1.6`, execute:
```bash
mamba install crisprme=2.1.6
```
If you're using `Conda`, replace `mamba` with `conda` in the commands above.

**Step 3: Verify the Update**

After the update completes, ensure the installation was successful by checking the 
version:
```bash
crisprme.py --version  # Confirm the installed version
```
If the displayed version matches the one you installed, the update was successful.

### 1.2 Install CRISPRme via Docker
---

This section is organized into two subsections to guide you through the setup of 
**CRISPRme** using Docker:  

- [Installing Docker](#121-installing-docker): 
<br>Provides step-by-step instructions for installing Docker on your system, 
ensuring compatibility with all operating systems, including Linux, macOS, and Windows.  

- [Building and Pulling CRISPRme Docker Image](#122-building-and-pulling-crisprme-docker-image): 
  <br>Explains how to create or download the CRISPRme Docker image to set up a 
  containerized environment for seamless execution.  

Follow the subsections in order if Docker is not yet installed on your machine. 
If Docker is already installed, skip to the second subsection.

#### 1.2.1 Installing Docker 
---

MacOS and Windows users are encouraged to install [Docker](https://www.docker.com/get-started) 
to use CRISPRme. Linux users may also choose Docker for convenience and compatibility.

Docker provides tailored distributions for different operating systems. Follow the 
official Docker installation guide specific to your OS:

- [MacOS Installation Guide](https://docs.docker.com/docker-for-mac/install/)
- [Windows Installation Guide](https://docs.docker.com/docker-for-windows/install/)
- [Linux Installation Guide](https://docs.docker.com/engine/install/ubuntu/)

**Linux-Specific Post-Installation Steps**

If you're using Linux, additional configuration steps are required:
1. Create the Docker Group:
  ```bash
  sudo groupadd docker
  ```

2. Add Your User to the Docker Group:
  ```bash
  sudo usermod -aG docker $USER
  ```
  Repeat this command for any additional users you want to include in the Docker 
  Group.

3. Restart Your Machine
  <br>Log out and log back in, or restart your machine to apply the changes.

**Testing Docker Installation**

Once Docker is installed, verify the setup by opening a terminal window and typing:
```bash
docker run hello-world
```

If Docker is installed correctly, you should see output like this:
```
Hello from Docker!
This message shows that your installation appears to be working correctly.

To generate this message, Docker took the following steps:
 1. The Docker client contacted the Docker daemon.
 2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
 3. The Docker daemon created a new container from that image, which runs the executable that produces this output.
 4. The Docker daemon streamed this output to the Docker client, which displayed it on your terminal.

For more examples and ideas, visit:
 https://docs.docker.com/get-started/
```

#### 1.2.2 Building and Pulling CRISPRme Docker Image
---

After installing Docker, you can download and build the CRISPRme Docker image by 
running the following command in a terminal:
```bash
docker pull pinellolab/crisprme
```

This command retrieves the latest pre-built CRISPRme image from Docker Hub and sets 
it up on your system, ensuring all required dependencies and configurations are 
included.

Once the download is complete, the CRISPRme Docker image will be ready for use. 
To confirm the image is successfully installed, you can list all available Docker 
images by typing:
```bash
docker images
```

Look for an entry similar to the following:
```
REPOSITORY          TAG       IMAGE ID       CREATED        SIZE
pinellolab/crisprme latest    <image_id>     <timestamp>    <size>
```

You are now ready to run CRISPRme using Docker.


## Usage

CRISPRme offers a range of functionalities for efficient variant- and haplotype-aware CRISPR off-target searches, along with tools for further analysis of the results. Additionally, it provides a graphical interface that can be deployed locally to simplify usage.

### System requirements

CRISPRme typically requires at least 32 GB of available memory (RAM) for most cases. However, for intensive searches involving complete genomes and large variant datasets, more memory may be necessary. In such cases, we recommend ensuring that the system running CRISPRme has **at least 64 GB of memory available**.

### CRISPRme directory structure

CRISPRme is designed to work within a specific directory tree structure in the current working directory. The tool requires five main directories to be properly organized for it to run correctly:

- **Genomes**: This directory contains subdirectories, each holding the reference genome to be used. Reference genomes must be split by chromosome, with each chromosome sequence in a separate FASTA file.

- **VCFs**: This directory contains subdirectories, each holding the VCF files, which should be compressed using `bgzip` (with a `.gz` extension). Similar to the Genomes directory, VCF files must be split by chromosome, with one VCF file per chromosome.

- **sampleIDs**: This directory contains tab-separated files listing the samples for each VCF dataset.

- **Annotation**: This directory contains annotation files in BED format.

- **PAMs**: This directory contains text files specifying the PAM sequences.

The image below provides a visual representation of the directory structure required by CRISPRme:
![fig1](./docs/readme/crisprme_dirtree.png)

### CRISPRme functionalities

This section provides an overview of CRISPRme's core functionalities, including descriptions of each feature, the required input data and formats, and the resulting output data. Below is a list of CRISPRme's key functionalities:

- [**complete-search**](#complete-search): Conducts a comprehensive search across the entire genome (both reference and variant data, if requested), performs Cutting Frequency Determination (CFD) analysis, and selects targets.

- [**complete-test**](#complete-test): Tests the full CRISPRme pipeline using a small input dataset, allowing you to verify the tool's functionality before running larger analyses.

- [**targets-integration**](#targets-integration): Integrates in-silico predicted targets with empirical data to generate a usable target panel.

- [**gnomAD-converter**](#gnomad-converter): Converts gnomAD v3.1 vcf.bgz files into a format compatible with CRISPRme.

- [**generate-personal-card**](#generate-personal-card): Generates a personalized card for a specific sample, extracting all private off-targets unique to that individual.

- [**web-interface**](#web-interface): Activates CRISPRme's web interface, allowing you to interact with the tool through a browser locally.

#### Complete-search

The `complete-search` functionality executes a comprehensive variant- and haplotype-aware off-target analysis on the provided reference genome and variant datasets. This function processes all steps in CRISPRme's pipeline, including off-target identification, annotation, and reporting. It generates detailed reports on population- and sample-specific off-targets, as well as graphical summaries of the findings. This section describes the input parameters, output data, and provides a usage example.

##### Input arguments

Below the list and description of the input arguments required by `complete-search` function:

- **--help**: Displays the help message and exits.

- **--genome**: Specifies the directory containing the reference genome in FASTA format. The genome should be divided into individual chromosome files (e.g., `chr1.fa`, `chr2.fa`, etc.).

- **--vcf**: ext file specifying the subdirectories containing the VCFs used to enrich the input reference genome. If not specified, CRISPRme searches off-targets only on the reference genome sequence. [OPTIONAL]

- **--guide**: Text file containing a list of guide sequences (one per row) to search within the input data. This option is mutually exclusive with the `--sequence` argument.

- **--sequence**: Text file listing guide sequences formatted as a FASTA file. This option cannot be used in conjunction with `--guide`.

- **--pam**: Text file containing the PAM sequence to be used in the search.

- **--be-window**: Specifies the base editor window, requiring start and stop coordinates (1-based, with respect to the target start), separated by a comma. This option is used to search for susceptibility to base editing. [OPTIONAL]

- **--be-base**: Specifies the base editor nucleotide(s) to check for the chosen editor. [OPTIONAL]

- **--annotation**: BED file containing annotation data, such as a list of genetic regions with associated functions (e.g., DNase hypersensitive sites, enhancers, etc.).

- **--samplesID**: Text file containing a list of sample IDs (one per line) corresponding to the VCF datasets used. These files must be located within the sampleID directory. This option is mandatory if `--vcf` is specified. [OPTIONAL]

- **--gene_annotation**: BED file containing gene data, such as Gencode annotations, used to compute gene proximity for the identified targets.

- **--mm**: Maximum number of mismatches allowed during the target search.

- **--bDNA**: Maximum bulge size allowed on DNA.

- **--bRNA**: Maximum bulge size allowed on RNA.

- **--merge**: Specifies the window size (in base pairs) used to merge candidate off-targets. When merging, the pivot targets are selected based on the highest score, particularly when CFD and CRISTA scores are computed, using the criteria defined by the `--sorting-criteria-scoring` argument. If CFD and CRISTA scores are not available, the merging follows the criteria specified in `--sorting-criteria`. The default window size for merging is set to 3 base pairs. [OPTIONAL]

- **--sorting-criteria-scoring**: Defines the criteria for sorting and merging targets based on CFD/CRISTA scores, prioritizing the highest score. The criteria must be provided as a comma-separated list. Available options include `mm` (number of mismatches), `bulges` (bulge size), and `mm+bulges` (total number of mismatches and bulges). By default, `mm+bulges` is used as the sorting criterion. [OPTIONAL]

- **--sorting-criteria**: Specifies the criteria for sorting and merging targets based on the fewest mismatches and bulges, particularly when CFD and CRISTA scores cannot be computed. The criteria should be provided as a comma-separated list. Available options include `mm` (number of mismatches), `bulges` (bulge size), and `mm+bulges` (total mismatches and bulges). By default, the tool uses `mm+bulges` as the primary criterion and `mm` as a secondary criterion, with `mm+bulges` serving as the main pivot for merging. [OPTIONAL]

- **--output**: Name of the output directory, which will contain all the results of CRISPRme's analysis. This directory will be located within the `Results` directory.

- **--thread**: Number of threads to use during the computation. By default, CRISPRme uses 4 threads. [OPTIONAL]

- **debug**: Runs the tool in debug mode. [OPTIONAL]

##### Output data

The `complete-search` functionality generates various reports detailing the targets identified and prioritized by CRISPRme, along with statistical summaries of these targets categorized by sample or input guide. Additionally, CRISPRme provides graphical summaries of the search results and statistics in image format.

Below is a description of the contents of each output file:

- **`*.bestMerge.txt`**: This file contains the top targets selected based on the input criteria. Each row represents the best target for a specific genomic location, chosen according to the highest CFD score (if computed), highest CRISTA score (if computed), and the fewest combined mismatches and bulges.

- **`*.altMerge.txt`**: This file includes all alternative alignments that were not selected for the `*.bestMerge.txt` file. It lists all additional targets identified by CRISPRme at each candidate off-target position. Each genomic position may have multiple reported targets in this file.

- **`*.integrated_results.tsv`**: This file details the best targets according to the input criteria, along with annotation data from the input files. It includes information on the gene proximity of each candidate off-target and its location within regulatory elements.

- **`*.all_results_with_alternative_alignments.tsv`**: This file lists all targets identified by CRISPRme, including alternative alignments at potential off-target sites. Each reported target is annotated with gene proximity and its location within functional genomic regulatory elements, based on the input annotation files.

- **`*.summary_by_guide.<guide-sequence>_CFD.txt`**: This file summarizes the number of candidate off-targets identified by CRISPRme for each guide, using the CFD score as the primary sorting criterion. It provides counts of targets based on bulge type, mismatch number, and bulge size, as well as how many targets were found in the reference sequence, variant genome, and how many were a result of PAM creation due to a variant. These counts are derived from the targets listed in the `*.bestMerge.txt` file.

- **`*.summary_by_guide.<guide-sequence>_CRISTA.txt`**: This file provides a summary similar to the `*_CFD.txt` file, but uses the CRISTA score as the primary sorting criterion. It details the number of candidate off-targets per guide, categorized by bulge type, mismatch number, bulge size, and the presence of targets in the reference sequence, variant genome, and those resulting from PAM creation due to a variant.

- **`*.summary_by_guide.<guide-sequence>_fewest.txt`**: This file summarizes the number of candidate off-targets identified by CRISPRme for each guide, using the fewest combined mismatches and bulge size as the primary sorting criterion. It provides counts based on bulge type, mismatch number, and bulge size, and notes how many targets were found in the reference sequence, variant genome, and those due to PAM creation by a variant.

- **`*.summary_by_samples.<guide-sequence>_CFD.txt`**: This file counts the candidate off-targets identified by CRISPRme that are private to each sample, using the CFD score as the primary sorting criterion. It also reports how many targets were found in the sample’s population and superpopulation, and how many involved PAM creation induced by a variant.

- **`*.summary_by_samples.<guide-sequence>_CRISTA.txt`**: Similar to the `*_CFD.txt` file, this file counts private candidate off-targets for each sample, using the CRISTA score as the primary sorting criterion. It also provides details on targets found in the population, superpopulation, and those involving PAM creation due to a variant.

- **`*.summary_by_samples.<guide-sequence>_fewest.txt`**: This file summarizes the candidate off-targets identified by CRISPRme that are private to each sample, sorted by the fewest combined mismatches and bulge size. It includes counts of targets found in the sample’s population and superpopulation, as well as those involving PAM creation due to a variant.

- **`imgs` directory**: This directory contains images that illustrate the impact of genetic variants on the top 1000 targets, selected according to each of the previously mentioned criteria (CFD score, CRISTA score, and the combined number of mismatches and bulge size). The images depict changes in these metrics due to genetic variations. Additionally, the directory includes bar plots that show the distribution of targets across different populations, categorized by the combined number of mismatches and bulge size identified by CRISPRme.

##### Usage example

Usage example of `complete-search` function:

- Via `conda`/`mamba`:
  ```
  crisprme.py complete-search --genome Genomes/hg38/ --vcf vcf_list.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/dhs+gencode+encode.hg38.bed --samplesID samplesID_list.txt --be-window 4,8 --be-base A --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --sorting-criteria-scoring mm+bulges --sorting-criteria mm+bulges,mm --output sg1617_NGG --thread 4
  ```

- Via `Docker`:
  ```
  docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crisprme crisprme.py complete-search --genome Genomes/hg38/ --vcf vcf_list.txt/ --guide sg1617.txt --pam ./PAMs/20bp-NGG-SpCas9.txt --annotation Annotations/dhs+encode+gencode.hg38.bed --samplesID samplesID_list.txt --be-window 4,8 --be-base A --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --sorting-criteria-scoring mm+bulges --sorting-criteria mm+bulges,mm --output sg1617_NGG --thread 4
  ```

#### Complete-test

The `complete-test` feature offers a fully automated pipeline for testing CRISPRme's installation. It handles the creation of the CRISPRme directory structure and downloads the essential files required to run the tool. `Complete-test` provides several testing options, allowing the user to perform quick tests on a single chromosome or conduct comprehensive tests on the entire genome. Furthermore, the feature supports testing CRISPRme using either the 1000 Genomes Phase 3 dataset or the Human Genome Diversity Project (HGDP) dataset. The testing environment can be configured by the user via command line parameters.

The necessary data for these tests can be accessed [here](https://github.com/pinellolab/CRISPRme/tree/gnomad-4.1-converter/test/data), except for the FASTA and VCF files, which are downloaded from the UCSC and 1000 Genomes/HGDP databases due to their large file sizes.

##### Input arguments

Below the list and description of the input arguments required by `complete-test` function:

- **--help**: Displays the help message and exits.

- **--chrom**: Specifies the chromosome for CRISPRme's test. The chromosome number must be prefixed with chr (e.g., chr22). By default, the test is conducted on the entire genome. [OPTIONAL]

- **--vcf_dataset**: Defines the variant dataset used for testing CRISPRme. The available options are `1000G` for the 1000 Genomes dataset and `HGDP` for the Human Genome Diversity Project dataset. The default dataset is 1000 Genomes. [OPTIONAL]

- **--debug**: Runs the tool in debug mode. [OPTIONAL]

##### Output data


##### Usage example

Usage example of `complete-test` function:

- Via `conda`/`mamba`:
```
# test only on chr 22 and 1000 Genomes
crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G

# test on the entire genome and HGDP 
crisprme.py complete-test --vcf_dataset HGDP

# test on the entire genome and 1000G+HGDP
crisprme.py complete-test --vcf_dataset 1000G+HGDP
```

- Via `Docker`:
```
# test only on chr 22 and 1000 Genomes
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crisprme crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G

# test on the entire genome and HGDP 
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crisprme crisprme.py complete-test --vcf_dataset HGDP

# test on the entire genome and 1000G+HGDP
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crisprme
crisprme.py complete-test --vcf_dataset 1000G+HGDP
```

#### Targets-integration

The `targets-integration` feature combines the targets identified by CRISPRme with empirically validated targets, such as those obtained through GUIDE-seq or CIRCLE-seq. The empirically validated targets should be provided in a `BED` file. After integrating the two datasets, this feature generates an `integrated_result` file within the specified output directory.

##### Input arguments

Below the list and description of the input arguments required by `targets-integration` function:

- **--help**: Displays the help message and exits.

- **--targets**: File containing the targets identified and processed from a CRISPRme search.

- **--empirical_data**: `BED` file with empirically verified off-targets, such as those obtained from GUIDE-seq, CIRCLE-seq, or other sequencing protocols.

- **--output**: NName of the output directory where the integrated targets file will be saved.

- **--debug**: Runs the tool in debug mode. [OPTIONAL]

##### Output data

##### Usage example

Usage example of `targets-integration` function:

- Via `conda`/`mamba`:
  ```
  crisprme.py targets-integration --targets results.integrated_results.tsv --empirical_data empirical_data.bed --output integrated_targets_dir
  ```

- Via `Docker`:
  ```
  docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py targets-integration --targets results.integrated_results.tsv --empirical_data empirical_data.bed --output integrated_targets_dir
  ```

#### gnomAD-converter 

The `gnomAD-converter` function provides a utility to convert gnomAD VCF files into a format suitable for use as input in CRISPRme searches. The converter currently supports gnomAD VCF files for the following versions:

- v3.1

- v4.0

- v4.1 (including joint VCFs)

For a successful conversion, the converter requires a sample ID file. Since individual samples are not provided in gnomAD VCFs, CRISPRme treats populations as single samples. This approach does not impact population-based statistics. However, if you require sample-specific statistics, we do not recommend using gnomAD, as the data are population-based rather than sample-based. Sample ID files for gnomAD [v3.1 and v4.0](https://github.com/pinellolab/CRISPRme/blob/gnomad-4.1-converter/test/data/hg38_gnomAD.samplesID.txt) and [v4.1](https://github.com/pinellolab/CRISPRme/blob/gnomad-4.1-converter/test/data/hg38_gnomAD.v4.1.samplesID.txt) are available on our GitHub page.

##### Input arguments

Below the list and description of the input arguments required by `targets-integration` function:

- **--help**: Displays the help message and exits.

- **--gnomAD_VCFdir**: Directory containing the gnomAD VCF files to process and convert into a format suitable for CRISPRme.

- **--samplesID**: Text file containing the sample IDs used during the gnomAD VCF conversion. In this file, gnomAD populations are treated as individual samples to construct VCFs suitable for CRISPRme.

- **--joint**: Flag specifying whether the gnomAD VCFs to convert are joint VCFs. By default the converter assumes the input VCFs are not joint VCFs [OPTIONAL]

- **--keep**: Flag specifying whether to retain all variants, regardless of their value in the FILTER field of the VCF. By default, the converter discards variants without a PASS value in the FILTER field. [OPTIONAL]

- **--multiallelic**: Flag indicating whether to merge variants mapped at the same position into a single line, creating multiallelic entries. By default, the converter keeps all variants biallelic and does not merge them. [OPTIONAL]

- **--thread**: Specifies the number of threads to use during the computation. By default, CRISPRme uses 4 threads. [OPTIONAL]

- **--debug**: Runs the tool in debug mode. [OPTIONAL]

##### Output data

##### Usage example

Usage example of `targets-integration` function:

- Via `conda`/`mamba`:
  ```
  # keep all variants and merge variants in multiallelic variants
  crisprme.py gnomAD-converter --gnomAD_VCFdir gnomad_vcf_dir --samplesID hg38_gnomAD.samplesID.txt --keep --multiallelic -thread 4
  ```

- Via `Docker`:
  ```
  # keep all variants and merge variants in multiallelic variants
  docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py gnomAD-converter --gnomAD_VCFdir gnomad_vcf_dir --samplesID hg38_gnomAD.samplesID.txt --keep --multiallelic -thread 4
  ```

#### Generate-personal-card

The `generate-personal-card` functionality generates a sample-specific report. This report, called personal card, contains all targets found by CRISPRme on the sample-specific genomic sequence, accounting for its private variants. This functionality is useful when searching private potential off-target sequences for a certain input guide.

##### Input arguments 

Below the list and description of the input arguments required by `generate-personal-card` function:

- **--help**: Displays the help message and exits.

- **--result_dir**: Directory containing the CRISPRme search results. Targets are extracted from the targets reports available in the input directory.

- **--guide_seq**: sequence of the guide of interest to use to extract the sample-specific targets from the input data

- **--sample_id**: ID of the sample to use to generate the personal card

- **--debug**: Runs the tool in debug mode.

##### Output data
- Set of plots generated with personal and private targets containing the variant CFD score and the reference CFD score
- Filtered file with private targets of the sample directly extracted from integrated file

##### Usage example

Usage example of `generate-personal-card` function:

- Via `conda`/`mamba`:
  ```
  crisprme.py generate-personal-card --result_dir Results/sg1617.6.2.2 --guide_seq CTAACAGTTGCTTTTATCACNNN --sample_id NA21129
  ```

- Via `Docker`:
  ```
  docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py generate-personal-card --result_dir Results/sg1617.6.2.2/ --guide_seq CTAACAGTTGCTTTTATCACNNN --sample_id NA21129
  ```

#### Web-interface

The `web-interface` deploys locally a graphical user interface similar to crisprme's website. This local server can be accessed using the major web-browser used (successfully tested on Chrome, Safari, and Firefox). The graphical interface provides an easy-to-use interface to run CRISPRme and explore the search results interactively. 

Note that this functionality is currently available only for CRISPRme distributions installed via `conda`/`mamba`.

##### Input arguments

Below the list and description of the input arguments required by `web-interface` function:

- **--debug**: Deploys the local server in debug mode.

##### Output data

This function does not produce output data.

##### Usage example

Usage example of `web-interface` function:

- Via `conda`/`mamba`:
  ```
  crisprme.py web-interface  # starts local server
  ```

## Test

To quickly verify the installation of CRISPRme, open a new terminal window and run the following command:

- Via `Conda`/`Mamba`:
```
crisprme.py --version  # v2.1.6
```

- Via `Docker`:
```
docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py --version  # v2.1.6
```

If the command returns the correct software version, CRISPRme has been successfully installed and is ready for use.

For a more detailed test, testing the main CRISPRme's functionality (`complete-search`), CRISPRme provides a dedicated [functionality](#complete-test). The test functionality provides two main test options, one quicker the other more detailed. The quicker test mode, performs an off-target search on a single chromosome sequence enriched with variants from either 1000 Genomes or Human Genome Diversity Project datasets. The more detailed test, instead, test crisprme searching potential off-targets across the entire genome. The test on the full genome replicates the analysis presented in our [paper](#citation), using the sg1617 single uguide RNA, `NGG` PAM, and the gencode and encode annotations.

To run the test on a single chromosome (chromosome 22) and using 1000 Genomes variants, open a new terminal window and run the following command:

- Via `Conda`/`Mamba`
```
crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G
```

- Via `Docker`:
```
docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py complete-test --chrom 22 --vcf_dataset 1000G
```

To run the test on the entire genome and using 1000 Genomes variants, open a new terminal window and run the following command:

- Via `Conda`/`Mamba`
```
crisprme.py complete-test --vcf_dataset 1000G
```

- Via `Docker`:
```
docker run -v ${PWD}:/DATA -w /DATA -i i pinellolab/crisprme crisprme.py complete-test --vcf_dataset 1000G
```

## Citation
If you use CRISPRme in your research, please cite our [paper](https://rdcu.be/c1GYQ):

Cancellieri S, Zeng J, Lin LY, Tognon M, Nguyen MA, Lin J, Bombieri N, Maitland SA, Ciuculescu MF, Katta V, Tsai SQ, Armant M, Wolfe SA, Giugno R, Bauer DE, Pinello L. Human genetic diversity alters off-target outcomes of therapeutic gene editing. Nat Genet. 2023 Jan;55(1):34-43. doi: [10.1038/s41588-022-01257-y](https://doi.org/10.1038/s41588-022-01257-y). Epub 2022 Dec 15. PMID: 36522432; PMCID: PMC10272994.

## Contacts

Luca Pinello: lpinello@mgh.harvard.edu<br>
Rosalba Giugno: rosalba.giugno@univr.it<br>
Daniel Bauer: bauer@bloodgroup.tch.harvard.edu

## License

CRISPRme is licensed under AGPL-3.0, allowing use for academic research only.

For-profit institutions must purchase a license before using CRISPRme. For more information, please contact lpinello@mgh.harvard.edu.
