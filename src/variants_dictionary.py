"""
This module provides functions for managing and processing variant and indel
dictionaries in genomic data. It includes utilities for searching allele frequencies,
formatting SNP and indel entries, and saving the results to various file formats.

Key functionalities include:
- Extracting and formatting allele frequencies and SNP identifiers.
- Creating and saving FASTA and log files for indels.
- Updating variant dictionaries with information from SNPs and indels.
- Handling errors related to file operations and data formatting.
"""

from typing import List, Tuple, Dict, Union

import json
import os

PADSIZE = 26
LOGINDELSHEADER = ["CHR", "SAMPLES", "rsID", "AF", "indel", "FAKEPOS", "refseq"]


def search_af(info: str, snpid: str) -> List[str]:
    """
    Searches for the allele frequency (AF) in the provided INFO field of a variant.
    This function extracts the AF values associated with the specified SNP ID from
    the INFO string.

    Args:
        info (str): The INFO field containing variant information.
        snpid (str): The identifier for the SNP being queried.

    Returns:
        List[str]: A list of allele frequency values extracted from the INFO field.

    Raises:
        ValueError: If the INFO field is empty or if no AF values are found for the
            specified SNP ID.
    """

    if info == ".":
        raise ValueError(f"INFO field empty for variant {snpid}. Cannot recover AF")
    # split info field -> ';' separator char
    fields = info.strip().split(";")
    idx = next((i for i, e in enumerate(fields) if e.startswith("AF")), -1)
    if idx == -1:  # is AF not present?
        raise ValueError(f"No AF found for variant {snpid}")
    return fields[idx].replace("AF=", "").split(",")


def retrieve_af(afs: List[str], allele_id: int) -> float:
    """
    Retrieves the allele frequency from a list of allele frequencies using the
    specified index. This function converts the value at the given index to a
    float and handles potential errors related to non-numeric values.

    Args:
        afs (List[str]): A list of allele frequency values as strings.
        allele_id (int): The index of the allele frequency to retrieve.

    Returns:
        float: The allele frequency as a float.

    Raises:
        ValueError: If the value at the specified index is not numeric.
    """
    try:
        return float(afs[allele_id])  # recover allele frequency
    except TypeError as e:
        raise ValueError(
            f"The Allele Frequency ({afs[allele_id]}) is not numeric"
        ) from e


def format_snpkey(chrom: str, pos: int) -> str:
    """
    Formats a SNP identifier for use as a key in a dictionary. This function combines
    the chromosome and position into a single string, separated by a comma.

    Args:
        chrom (str): The chromosome where the SNP is located.
        pos (int): The position of the SNP on the chromosome.

    Returns:
        str: The formatted SNP identifier in the form "chrom,pos".
    """

    # format snp identifier for dict (variant key)
    return f"{chrom},{pos}"


def format_snpentry(
    samples: List[str], ref: str, alt: str, snpid: str, af: float
) -> str:
    """
    Creates a formatted entry for a SNP based on the provided sample information.
    The entry includes the samples carrying the variant, reference and alternative
    alleles, SNP ID, and allele frequency, formatted appropriately based on whether
    any samples carry the variant.

    Args:
        samples (List[str]): A list of sample identifiers that carry the variant.
        ref (str): The reference allele.
        alt (str): The alternative allele.
        snpid (str): The identifier for the SNP.
        af (float): The allele frequency of the SNP.

    Returns:
        str: The formatted SNP entry as a string, structured based on the presence
            of samples carrying the variant.
    """

    # create the netry pointed by the corresponding snp key
    # the entry will look in two different ways depending on whether any sample
    # carries the variant or not
    #
    # If samples carry the variants the entry will look like:
    # chrX,100 -> sample1,sample2,sample3;A,T;rs123;0.01
    #
    # if no sample found the dictionary entry will look like:
    # chrX,100 -> ;A,T;rs123;0.01
    return ";".join([",".join(samples), f"{ref},{alt}", snpid, str(af)])


def retrieve_alternative_alleles(alt: str) -> Tuple[List[Tuple[str, int]], int]:
    """
    Extracts alternative alleles from a given string representing alternative alleles
    at a variant site. This function returns a list of tuples containing each
    alternative allele and its corresponding identifier, along with the total number
    of alternative alleles.

    Args:
        alt (str): A comma-separated string of alternative alleles.

    Returns:
        Tuple[List[Tuple[str, int]], int]: A tuple containing a list of tuples
            with each alternative allele and its identifier, and the total number
            of alternative alleles.
    """

    # split alt to retrieve all alternative alleles for multiallelic sites (snvs + indels)
    alleles_alt = alt.split(",")
    # recover the original number of alternative alleles
    aa_num = len(alleles_alt)  # required to match ids and allele frequencies
    # record allele number to recover the original genotype, e.g., 1,2,3,...,N
    return [(aa, i + 1) for i, aa in enumerate(alleles_alt)], aa_num


def retrieve_snpids(snpid: str, alleles_alt_num: int) -> List[str]:
    """
    Retrieves SNP identifiers for each alternative allele from a given SNP ID string.
    This function handles both single and multiple SNP IDs, ensuring that the correct
    number of identifiers is returned based on the number of alternative alleles.

    Args:
        snpid (str): A comma-separated string of SNP identifiers.
        alleles_alt_num (int): The expected number of alternative alleles.

    Returns:
        List[str]: A list of SNP identifiers corresponding to the alternative alleles.

    Raises:
        ValueError: If the SNP ID format is invalid for the specified number of
            alternative alleles.
    """

    # recover snp id for each alternative allele
    snpids = snpid.split(",")  # include snps and indels
    if len(snpids) == 1:  # single id?
        if snpids[0] == ".":  # missing value in VCF -> .
            return ["." for _ in range(alleles_alt_num)]
        else:
            raise ValueError(
                f"Forbidden snp id format ({snpid}) to denote {alleles_alt_num} alternative alleles"
            )
    return snpids


def retrieve_samples(
    samples: List[str], genotypes: List[str], aanum: int
) -> Dict[str, List[str]]:
    """
    Constructs a dictionary that links alternative alleles to their corresponding
    samples and genotypes. This function creates an efficient mapping of allele
    identifiers to lists of samples, each annotated with their respective genotype.

    Args:
        samples (List[str]): A list of sample identifiers.
        genotypes (List[str]): A list of genotype strings corresponding to the
            samples.
        aanum (int): The number of alternative alleles.

    Returns:
        Dict[str, List[str]]: A dictionary where each key is an allele identifier
            and each value is a list of samples with their genotypes.
    """

    # dictionary to efficiently link alleles and samples
    samples_dict = {
        f"{aid + 1}": [] for aid in range(aanum)
    }  # avoid initialize values with None for typing requirements
    for aid in samples_dict:  # consider snp and indels
        samples_dict[aid] = [
            f"{samples[i]}:{gt}" for i, gt in enumerate(genotypes) if aid in gt
        ]  # append genotype to each sample
    return samples_dict


def process_snp(
    samples_allele: List[str], ref: str, allele: str, snpid: str, allele_freq: float
) -> str:
    """
    Formats a SNP entry for inclusion in the variants dictionary. This function
    utilizes the provided sample allele information, reference allele, alternative
    allele, SNP ID, and allele frequency to create a properly formatted SNP entry.

    Args:
        samples_allele (List[str]): A list of samples carrying the allele, annotated
            with their genotypes.
        ref (str): The reference allele.
        allele (str): The alternative allele.
        snpid (str): The identifier for the SNP.
        allele_freq (float): The allele frequency of the SNP.

    Returns:
        str: The formatted SNP entry as a string.
    """

    # format the snp entry to insert in the variants dictionary
    # see format_snpentry() for info regarding the snp formatting style
    return format_snpentry(samples_allele, ref, allele, snpid, allele_freq)


def pad_genome_seq(pos: int, ref: str, n: int) -> Tuple[int, int]:
    """
    Computes the padded start and stop positions for a genomic sequence surrounding
    an indel. This function adjusts the positions by a specified number of base pairs
    upstream and downstream of the given position, taking into account the length of
    the reference allele.

    Args:
        pos (int): The position of the indel in the genome.
        ref (str): The reference allele sequence.
        n (int): The number of base pairs to pad upstream and downstream.

    Returns:
        Tuple[int, int]: A tuple containing the padded start and stop positions.
    """

    # compute start and stop position adding 26 bp up and downstream wrt indel pos
    return pos - n, pos + n + len(ref)


def retrieve_genome_seq(fasta: List[str], start: int, stop: int) -> List[str]:
    """
    Retrieves a subsequence from a given FASTA sequence based on specified start
    and stop positions.
    This function assumes that the provided FASTA list corresponds to the desired
    contig and returns the sequence within the specified range.

    Args:
        fasta (List[str]): The FASTA sequence represented as a list of characters.
        start (int): The starting index of the subsequence (0-based).
        stop (int): The ending index of the subsequence (exclusive).

    Returns:
        List[str]: A list containing the characters of the subsequence from the
            FASTA sequence.
    """

    # recover reference sequence at indel positon
    # assume fasta points to the desired contig
    return fasta[start:stop]


def compute_indel_sequence(genseq: List[str], ref: str, indel: str) -> Tuple[str, int]:
    """
    Computes the resulting sequence after inserting or deleting nucleotides based
    on the provided indel string.
    This function constructs a new sequence by combining parts of the original
    genome sequence with the specified indel.

    Args:
        genseq (List[str]): The original genome sequence represented as a list of
            characters.
        ref (str): The reference sequence that is being modified.
        indel (str): The insertion or deletion sequence to be applied.

    Returns:
        Tuple[str, int]: A tuple containing the modified sequence and the length
            of the indel.
    """

    # explode indel string in a list of str -> use faster list concatenation
    # instead of string concatenation
    indelseq = "".join(genseq[:25] + list(indel) + genseq[(25 + len(ref)) :])
    return indelseq, len(indelseq)


def format_indelkey(chrom: str, pos: int, ref: str, indel: str) -> str:
    """
    Creates a formatted identifier for an indel entry in a dictionary. This function
    combines the chromosome, position, reference allele, and indel sequence into a
    single string for easy indexing.

    Args:
        chrom (str): The chromosome where the indel is located.
        pos (int): The position of the indel on the chromosome.
        ref (str): The reference allele sequence.
        indel (str): The sequence of the indel.

    Returns:
        str: The formatted indel identifier in the form "chrom_pos_ref_indel".
    """

    # format indel identifier for dict (indel dict)
    return f"{chrom}_{pos}_{ref}_{indel}"


def format_indelentry(
    indelseq: str,
    chrom: str,
    start: int,
    stop: int,
    indelid: int,
    indel_samples: List[str],
    snpid: str,
    af: float,
    indelkey: str,
    indelstart: int,
    indelstop: int,
    refseq: List[str],
) -> Tuple[str, List[str]]:
    """
    Formats an indel entry for a genomic variant, creating a structured representation
    of the indel information.
    This function generates a tuple containing the indel sequence and a list of
    relevant details, including chromosome, position, samples, and reference sequence.

    Args:
        indelseq (str): The sequence of the indel.
        chrom (str): The chromosome where the indel is located.
        start (int): The start position of the indel.
        stop (int): The stop position of the indel.
        indelid (int): A unique identifier for the indel.
        indel_samples (List[str]): A list of samples associated with the indel.
        snpid (str): The SNP identifier related to the indel.
        af (float): The allele frequency of the indel.
        indelkey (str): A key representing the indel.
        indelstart (int): The start position of the indel in the reference sequence.
        indelstop (int): The stop position of the indel in the reference sequence.
        refseq (List[str]): The reference sequence represented as a list of
            characters.

    Returns:
        Tuple[str, List[str]]: A tuple containing the indel sequence and a list
            of formatted indel information.
    """

    # create the entry pointed by the corresponding indel key
    # the entry will look like:
    # chrX_100_A_GGGG -> (
    #   CCTTGGGGTTTT,
    #   [chrX_96-104_1, sample1,sample2,sample3, rs123, 0.01, chrX_100_A_GGGG, 0,12, CCTTATTTTGCA]
    # )
    indelinfos = [
        f"{chrom}_{start}-{stop}_{indelid}",
        ",".join([s.split(":")[0] for s in indel_samples]),
        snpid,
        str(af),
        indelkey,
        f"{indelstart},{indelstop}",
        "".join(refseq),
    ]
    return indelseq, indelinfos


def process_indel(
    fasta: List[str],
    chrom: str,
    pos: int,
    ref: str,
    indel: str,
    snpid: str,
    allele_freq: float,
    samples: List[str],
    indelstart: int,
    indelid: int,
) -> Tuple[str, Tuple[str, List[str]], int, int]:
    """
    Processes an indel by integrating it into the genomic sequence and creating an
    entry for it in the indel dictionary. This function retrieves the relevant genomic
    sequence, computes the new indel sequence, and formats the necessary information
    for storage.

    Args:
        fasta (List[str]): The genome sequence in FASTA format.
        chrom (str): The chromosome where the indel is located.
        pos (int): The position of the indel on the chromosome.
        ref (str): The reference allele sequence.
        indel (str): The sequence of the indel to be inserted.
        snpid (str): The identifier for the SNP associated with the indel.
        allele_freq (float): The allele frequency of the indel.
        samples (List[str]): A list of samples associated with the indel.
        indelstart (int): The starting position for the indel entry.
        indelid (int): The identifier for the indel.

    Returns:
        Tuple[str, Tuple[str, List[str]], int, int]: A tuple containing the indel
            key, the formatted indel entry, the next available stop position, and
            the next available indel identifier.
    """

    # initialize start and stop positions for the indel range
    start, stop = pad_genome_seq(pos, ref, PADSIZE)  # pad 26 bp up- and down-stream)
    # retrieve genome sequence at indel position
    # the sequence is extended by 26bp upstream and downstream
    genseq = retrieve_genome_seq(fasta, start, stop)
    # insert indel in the reference sequence, and return the new sequence size
    indelseq, indelsize = compute_indel_sequence(genseq, ref, indel)
    # retrieve ref sequence overlapping indel
    refseq = retrieve_genome_seq(fasta, start, start + indelsize)
    # compute the key for the current indel entry to insert it in the indels dict
    indelkey = format_indelkey(chrom, pos, ref, indel)
    # store start and stop positions of indel allele and create the entry for
    # the current indel
    indelstop = indelstart + indelsize
    indelentry = format_indelentry(
        indelseq,
        chrom,
        start,
        stop,
        indelid,
        samples,
        snpid,
        allele_freq,
        indelkey,
        indelstart,
        indelstop,
        refseq,
    )
    return (
        indelkey,
        indelentry,
        indelstop + 1,
        indelid + 1,
    )  # update indel start position and id


def update_variant_dictionaries(
    variant_dictionary: Dict[str, str],
    indels_dictionary: Dict[str, Tuple[str, List[str]]],
    fasta: List[str],
    chrom: str,
    snppos: int,
    snpid: str,
    ref: str,
    alt: str,
    info: str,
    genotypes: List[str],
    samples: List[str],
    indelstart: int,
    indelid: int,
) -> Tuple[Dict[str, str], Dict[str, Tuple[str, List[str]]], int, int]:
    """
    Updates the variant and indel dictionaries with information from a given SNP
    and its associated alleles.
    This function processes alternative alleles, retrieves relevant sample data,
    and constructs entries for both SNPs and indels, returning the updated dictionaries.

    Args:
        variant_dictionary (Dict[str, str]): A dictionary mapping SNP keys to their
            corresponding entries.
        indels_dictionary (Dict[str, Tuple[str, List[str]]]): A dictionary mapping
            indel keys to their corresponding entries.
        fasta (List[str]): The reference genome sequence represented as a list of
            characters.
        chrom (str): The chromosome where the variant is located.
        snppos (int): The position of the SNP in the genome.
        snpid (str): The SNP identifier.
        ref (str): The reference allele.
        alt (str): The alternate allele(s).
        info (str): Additional information about the variant.
        genotypes (List[str]): A list of genotype strings for the samples.
        samples (List[str]): A list of sample identifiers.
        indelstart (int): The start position of the indel.
        indelid (int): A unique identifier for the indel.

    Returns:
        Tuple[Dict[str, str], Dict[str, Tuple[str, List[str]]], int, int]: A
            tuple containing the updated variant dictionary, the updated indels
            dictionary, and the updated indel identifier and start positions.
    """

    # recover alternative alleles
    alleles_alt, alleles_num = retrieve_alternative_alleles(alt)
    # recover snp ids
    snpids = retrieve_snpids(snpid, alleles_num)
    # recover samples carrying each alternative allele
    samples_dict = retrieve_samples(samples, genotypes, alleles_num)
    # initialize snp key and entry list
    snpkey = format_snpkey(chrom, snppos)  # compute snp key
    snpentry = []  # snp entry -> separate different snps at same site with $
    # NOTE multiallelic sites with snps store each allele in the same entry
    # instead, indels are treated separately, for each indel allele a different
    # entry is created
    for allele, allele_id in alleles_alt:
        idx = allele_id - 1  # recover 0-based allele id -> required for list access
        samples_allele = samples_dict[f"{allele_id}"]  # samples carrying the allele
        allele_snpid = snpids[idx]  # current allele id
        # retrieve allele frequency for the current allele
        allele_freq = retrieve_af(search_af(info, allele_snpid), idx)
        if len(allele) == 1 and len(ref) == 1:  # snp
            snpentry.append(
                process_snp(samples_allele, ref, allele, allele_snpid, allele_freq)
            )
        else:  # indel
            # ignore unsolved alleles with format <*>
            # see chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://samtools.github.io/hts-specs/VCFv4.2.pdf
            # for further details on the allowed fields values
            if "<" in ref or "<" in alt:
                continue
            # compute indel key and entry for the current allele and update
            # indel's start position and id
            indelkey, indelentry, indelstart, indelid = process_indel(
                fasta,
                chrom,
                snppos,
                ref,
                allele,
                allele_snpid,
                allele_freq,
                samples_allele,
                indelstart,
                indelid,
            )
            indels_dictionary[indelkey] = indelentry
    variant_dictionary[snpkey] = "$".join(snpentry)  # update variant dictionary
    return variant_dictionary, indels_dictionary, indelid, indelstart


def save_variant_dictionary(
    variant_dictionary: Dict[str, str], contig: str, outdir: str
) -> None:
    """
    Saves the variant dictionary to a JSON file named according to the specified
    contig. This function handles serialization of the dictionary and manages
    potential errors that may arise during the file writing process.

    Args:
        variant_dictionary (Dict[str, str]): The dictionary containing variant
            information to be saved.
        contig (str): The name of the contig used to create the filename for the
            JSON file.
        outdir (str): The base directory where the SNPs dictionary will be saved.

    Raises:
        TypeError: If the data in the variant dictionary is not serializable to
            JSON.
        OSError: If there is an error writing to the file, such as permission issues
            or insufficient disk space.
        Exception: For any other unexpected errors that occur during the serialization
            process.
    """

    # variants dictionary is stored in a JSON file
    fname_json = os.path.join(outdir, f"my_dict_{contig}.json")
    try:
        with open(fname_json, mode="w") as f:
            json.dump(
                variant_dictionary, f
            )  # create JSON file from variants dictionary
    except TypeError as e:  # data is not serializable
        raise TypeError(f"Failed to serialize data to JSON: {e}") from e
    except OSError as e:  # permission error or disk full
        raise OSError(f"Failed to write to {fname_json}: {e}") from e
    except Exception as e:  # catch other exceptions
        raise Exception(f"Unexpected error occurred while serializing JSON: {e}") from e


def create_indels_fasta(
    indelfasta: str, indels: Union[List[str], None], contig: str
) -> None:
    """
    Creates a FASTA file containing indel sequences for a specified contig. This
    function formats the indels with a header and separates each sequence with a
    newline character followed by 'N'.

    Args:
        indelfasta (str): The path to the output FASTA file where the indels will
            be written.
        indels (List[str]): A list of indel sequences to be included in the FASTA
            file.
        contig (str): The name of the contig used in the FASTA header.

    Raises:
        OSError: If there is an error writing to the specified FASTA file.
    """

    # write indels fasta file
    indels = "" if indels is None else "\nN\n".join(indels) + "\nN\n"  # type: ignore -> suppress type warning
    try:
        with open(indelfasta, mode="w") as outfile:
            # write indels fasta header and content
            # each indel is separated by a N character placed below each sequence
            outfile.write(f">fake{contig}\n{indels}")
    except IOError as e:
        raise OSError(f"Failed writing indels to {indelfasta}") from e


def create_indels_log(logindels: str, indelsinfo: Union[List[str], None]) -> None:
    """
    Creates a log file containing information about indels in a tab-separated format.
    This function writes a header followed by the indel information to the specified
    log file.

    Args:
        logindels (str): The path to the output log file where the indel information
            will be written.
        indelsinfo (List[str]): A list of lists containing indel information to be
            logged.

    Raises:
        OSError: If there is an error writing to the specified log file.
    """

    # write log file with indels' info (tab-separated fields)
    indelsinfo = (
        "" if indelsinfo is None else "\n".join(["\t".join(e) for e in indelsinfo])
    )  # type: ignore -> suppress type warning
    try:
        with open(logindels, mode="w") as outfile:
            outfile.write("\t".join(LOGINDELSHEADER) + "\n")
            outfile.write(indelsinfo)  # type: ignore -> suppress type warning
    except IOError as e:
        raise OSError(f"Failed writing indels log to {logindels}") from e


def save_indels_dictionary(
    indels_dictionary: Dict[str, Tuple[str, List[str]]],
    contig: str,
    vcfdir: str,
    variants_genomedir: str,
    snps_genomedir: str,
) -> None:
    """
    Saves the indels dictionary to a specified output directory by creating a FASTA
    file and a log file. This function organizes the indel data and ensures that the
    necessary directories and files are created for storage.

    Args:
        indels_dictionary (Dict[str, Tuple[str, List[str]]]): A dictionary containing
            indel sequences and their associated information.
        contig (str): The name of the contig used in the filenames.
        outdir (str): The base directory where the output files will be saved.

    Raises:
        OSError: If there is an error creating the output directory or writing the
            files.
    """

    # unpack dictionary content
    indels, indelsinfo = None, None
    if indels_dictionary:
        indels, indelsinfo = zip(*indels_dictionary.values())
        indels, indelsinfo = list(indels), list(indelsinfo)
    # save indels in fasta file
    indelsdir = os.path.join(
        variants_genomedir, f"fake_{os.path.basename(vcfdir)}_{contig}"
    )
    if not os.path.isdir(indelsdir):
        os.mkdir(indelsdir)
    # create indels fasta and log
    indelfasta = os.path.join(indelsdir, f"fake{contig}.fa")
    create_indels_fasta(indelfasta, indels, contig)
    logindels = os.path.join(snps_genomedir, f"log{contig}.txt")
    create_indels_log(logindels, indelsinfo)
