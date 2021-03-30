Example of a possible call:

./automated_search.py --genome hg19_ref/ --vcf hg19_vcfs/ --guide guides.txt --pam pam.txt --annotation hg19.annotations.bed --samplesID samplesID.txt --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --output output/


INPUT explanation:

1) ref_folder = folder containing the reference genome divided in fasta files
2) [OPTIONAL] vcf_folder = folder containing the VCFs (with both SNP and INDELS). WARNING: the names of VCFs must contain the name of the chromosome in this fashion "something.chrNumber.somethingelse.anotherthing..." (the name of the chromosome must be in the second "field" delimited by the character "." ) 
3) guide_file = file contaning the guide(s)
4) pam_file = file containing the PAM used. The file must be structured as follows: repetitions of N (corresponding to the length of the hypotetical guide), followed by the PAM (or viceversa if the PAM is at the beginning of the guide), followed by a whitespace, followed by a number representing the length of the PAM
5) annotation_file = file containing the annotations
6) [OPTIONAL] sampleID = tab separated file containig information about the samples (population, superpopulation...) 
7) bMax = number of bulges used for indexing the genomes (both reference and variant)
8) mm = maximum number of mismatches allowed in the search phase
9) bDNA = number of DNA bulges permitted 
10) bRNA = number of RNA bulges permitted 
11) output_folder = folder used to store the results. Attention: this is also used to check whether some steps are already performed

