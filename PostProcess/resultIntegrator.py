#!/usr/bin/env python

from intervaltree import IntervalTree
import sys
import time
import glob
import subprocess


def rev_comp(a):
    # perform reverse complement on a string
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'


def createBedforMultiAlternative(variantList, samples):
    # 5096 conta degli alleli di 1000G
    # bcftools view -H -r 1:222531770 -s HG00001 originalVCFs/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
    # list containing snp position to search in the vcf
    toSearchList = list()
    # bed file with vcf data extracted
    vcfBedFile = outputDir + originFileName + '.temp.bed'
    # cycle over variant in list to extract the correct ones from vcf
    for elem in variantList:
        split = elem.strip().split('_')
        chrom = split[0]
        pos = split[1]
        vcfFile = ''
        for vcf in vcfList:
            if chrom+'.' in vcf:
                vcfFile = vcf
        vcfChr = chrom.replace('chr', '')
        toSearchList.append(vcfChr+':'+pos)
    # string with all the snps to search in the vcf file
    toSearchString = ','.join(toSearchList).strip()
    # subprocces run to call bcftools and extract the desidered snps
    subprocess.run(['bcftools', 'view', '-H', '-r',
                    toSearchString, '-s', samples, vcfFile, '-o', str(vcfBedFile)])
    haplotypeDict = dict()
    allele1 = list()
    allele2 = list()
    haplotype = 0
    # read line from vcf exctraction to obtain the alleles
    for line in open(vcfBedFile, 'r'):
        if 'INDEL' in line:
            continue
        split = line.strip().split('\t')
        tempHaploList = list()
        for elem in split:
            if '|' in elem:
                tempHaploList.append(elem)
        haplotypeDict[split[1]] = tempHaploList
        allele1 = [1]*len(tempHaploList)
        allele2 = [1]*len(tempHaploList)
    # count the correct haplotype frequence by and with sum operation
    # e.g. 0|1 and 1|1 will return 0|1 that will be counted as 1 for that sample haplotype
    for snp in haplotypeDict:
        for count, sample in enumerate(haplotypeDict[snp]):
            allele1_snp = sample.strip().split('|')[0]
            allele2_snp = sample.strip().split('|')[1]
            allele1[count] = int(allele1[count]) and int(allele1_snp)
            allele2[count] = int(allele2[count]) and int(allele2_snp)
    # final sum to generete haplotype frequency of the specific target
    for count, allele in enumerate(allele1):
        haplotype += int(allele1[count]) + int(allele2[count])
    # return to save the haplo freq divided by the number of alleles in 1000G (i.e. 2548samples*2alleles)
    return float(haplotype)/5096


print('READING INPUT FILES')

# READ INPUT FILES
crispritzResultFile = sys.argv[1]  # file with CRISPRitz results
empiricalResults = sys.argv[2]  # file with empirical results
annotationFile = sys.argv[3]  # annotation file used to find nearest gene
guideFile = sys.argv[4]  # real guide used in the search
outputDir = sys.argv[5]  # output file name
check = sys.argv[6].upper()  # output file name
genomeRelease = str(sys.argv[7]).strip()  # genome used in the search phase
# directory of vcf to perform haplotype count with more than one SVs in single target
vcfFileDirectory = sys.argv[8]

# OPEN INPUT FILES AND PREPARE OUTPUT FILE
# crispritz results file open
inCrispritzResults = open(crispritzResultFile, 'r')
# empirical-seq data open
inEmpiricalResults = open(empiricalResults, 'r').readlines()
# annotation data open
inAnnotationFile = open(annotationFile, 'r')
# guide file
# realguide = open(guideFile, 'r').readlines()
# open outputDir to write results
originFileName = crispritzResultFile.split('/')
originFileName = originFileName[len(originFileName)-1]
outFile = open(outputDir + originFileName +
               '.integrated_results.tsv', 'w')
# list file in vcf directory
checkVCF = False
if 'vuota' not in vcfFileDirectory:
    checkVCF = True
    vcfList = glob.glob(vcfFileDirectory+'/*.gz')
    redirectFile = open(outputDir + originFileName + '.redirectFile.out', 'w')
    for vcfFile in vcfList:
        # create tabix index if not already done
        subprocess.run(['tabix', str(vcfFile)], stderr=redirectFile)
    redirectFile.close()


empiricalTree = IntervalTree()
empiricalList = []
genomeDict = {}
empiricalDict = {}
valueDict = {}
# saveDict = {"real_guide": 'n', "genome": 'n', "chr": 'n', "prim_pos": 'n', "strand": 'n', "highest_CFD_guide_alignment": 'n', "highest_CFD_alignment(ref)": 'n',
#             "highest_CFD_alignment(alt)": 'n', "ref_seq_length": 'n', "ref_pos_alt(aligned_strand)": 'n', "pam": 'n', "annotation": 'n', "highest_CFD_score": 'n',
#             "highest_CFD_score(ref)": 'n', "highest_CFD_score(alt)": 'n', "risk_score": 'n', "absolute_risk_score": 'n', "highest_CFD_mismatch": 'n',
#             "highest_CFD_bulge": 'n', "highest_CFD_mismatch+bulge": 'n', "fewest_mm+bulge_guide_alignment": 'n', "fewest_mm+bulge_alignment(ref)": 'n',
#             "fewest_mm+bulge_alignment(alt)": 'n', "fewest_mm+bulge_CFD_score(ref)": 'n', "fewest_mm+bulge_CFD_score(alt)": 'n', "fewest_mismatch": 'n',
#             "fewest_bulge": 'n', "fewest_mismatch+bulge": 'n', "alt_haplotypes": 'n', "prim_origin": 'n', "prim_AF": 'n', "prim_samples": 'n',
#             "prim_SNP_ID(positive_strand)": 'n', "gene_name": 'n', "gene_ID": 'n', "gene_annotation": 'n', "gene_distance(kb)": 'n', "lowest_empirical": 'n',
#             "Nature2019": 'n', "Nature2019_mm+bul": 'n', "CHANGEseq": 'n', "CHANGEseq_mm+bul": 'n', "CIRCLEseq": 'n', "CIRCLEseq_mm+bul": 'n', "ONEseq": 'n',
#             "ONEseq_mm+bul": 'n', "GUIDEseq_293": 'n', "GUIDEseq_293_mm+bul": 'n', "GUIDEseq_CD34": 'n', "GUIDEseq_CD34_mm+bul": 'n', "GUIDEseq": 'n',
#             "GUIDEseq_mm+bul": 'n'}

saveDict = {"real_guide": 'n', "genome": 'n', "chr": 'n', "prim_pos": 'n', "strand": 'n', "highest_CFD_guide_alignment": 'n', "highest_CFD_alignment(ref)": 'n',
            "highest_CFD_alignment(alt)": 'n', "ref_seq_length": 'n', "ref_pos_alt(aligned_strand)": 'n', "pam": 'n', "annotation": 'n', "highest_CFD_score": 'n',
            "highest_CFD_score(ref)": 'n', "highest_CFD_score(alt)": 'n', "risk_score": 'n', "absolute_risk_score": 'n', "highest_CFD_mismatch": 'n',
            "highest_CFD_bulge": 'n', "highest_CFD_mismatch+bulge": 'n', "fewest_mm+bulge_guide_alignment": 'n', "fewest_mm+bulge_alignment(ref)": 'n',
            "fewest_mm+bulge_alignment(alt)": 'n', "fewest_mm+bulge_CFD_score(ref)": 'n', "fewest_mm+bulge_CFD_score(alt)": 'n', "fewest_mismatch": 'n',
            "fewest_bulge": 'n', "fewest_mismatch+bulge": 'n', "alt_haplotypes": 'n', "prim_origin": 'n', "prim_AF": 'n', "prim_samples": 'n',
            "prim_SNP_ID(positive_strand)": 'n', "gene_name": 'n', "gene_ID": 'n', "gene_annotation": 'n', "gene_distance(kb)": 'n', "lowest_empirical": 'n'}

start_time = time.time()

print('CREATING INTERVAL TREES')

for count, line in enumerate(inEmpiricalResults):
    empList = line.strip().split('\t')
    empList = [elem.strip() for elem in empList]
    empList.append(count)
    # adding empirical data to the tree
    empiricalTree[int(empList[1]):int(empList[2])] = empList
    # generate temp empirical dict to save empirical information per row
    empiricalDict[str(empList[4])] = 50
    # value of empirical row, keeping info about mm+bul for each empirical origin (seq data, in-silico, vitro, ecc)
    valueDict[str(empList[4])] = 'n'
    # update save dict with user-defined names from empirical data
    saveDict[str(empList[4])] = 'n'
    newkey = str(empList[4])+'_mm+bul'
    saveDict[newkey] = 'n'

# writing header in file
save = '#'
save += '\t'.join(list(saveDict.keys()))
# for key in saveDict:
#     save += str(key)+'\t'
save += '\n'
outFile.write(save)

print('INTEGRATING RESULTS')

if '#' in inCrispritzResults.readline():
    print('SKIP HEADER')
else:
    inCrispritzResults.seek(0)

for nline, line in enumerate(inCrispritzResults):
    x = line.strip().split('\t')
    annotationLine = inAnnotationFile.readline().strip().split('\t')
    lowestEmpirical = 100

    for key in saveDict:
        saveDict[key] = 'n'

    for key in empiricalDict:
        empiricalDict[key] = 50
        valueDict[key] = 'n'

    if 'NA' not in annotationLine and check == 'TRUE':
        for elem in annotationLine:
            if 'gene_id' in elem:
                temp = elem.strip().split(';')
                for name in temp:
                    if 'gene_id' in name:
                        saveDict['gene_ID'] = name.strip().split('=')[1]
                    if 'gene_name' in name:
                        saveDict['gene_name'] = name.strip().split('=')[1]
                    if 'gene_type' in name:
                        tipo = str(annotationLine[11])
                        if 'codon' in tipo or 'exon' in tipo:
                            tipo = 'CDS'
                        elif 'gene' in tipo or 'transcript' in tipo:
                            tipo = 'intron'
                        elif 'five_prime_UTR' in tipo:
                            tipo = "5'UTR"
                        elif 'three_prime_UTR' in tipo:
                            tipo = "3'UTR"
                    # saveDict['gene_annotation'] = str(annotationLine[11])
                        saveDict['gene_annotation'] = tipo
        saveDict['gene_distance(kb)'] = str(
            float(annotationLine[len(annotationLine)-1])/1000)
        if float(annotationLine[len(annotationLine)-1]) != 0:
            saveDict['gene_annotation'] = 'intergenic'

    origin = ''
    if 'n' in str(x[13]):
        origin = 'ref'
    else:
        origin = 'alt'

    variantList = ['n']
    if 'alt' in origin and str(x[20]) != str(x[21]):
        variantList = str(x[18]).strip().split(',')
        if len(variantList) > 1 and checkVCF:
            samples = x[13]
            x[17] = createBedforMultiAlternative(variantList, samples)
            x[17] = str(x[17])[:7]
        var_pos = []
        # generate variant position corrected to be in the positive strand
        if '+' in str(x[7]):
            refseq = str(x[3])
            altseq = str(x[2])
            for pos, nucleotide in enumerate(refseq):
                if altseq[pos].lower() != nucleotide.lower():
                    var_pos.append(pos+1)
        else:
            refseq = str(x[3])
            altseq = str(x[2])
            for pos in range(len(refseq)-1, -1, -1):
                if refseq[pos].lower() != altseq[pos].lower():
                    var_pos.append(pos+1)
        # generate variant nucleotide reverse complemented always in positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split('_')
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if split_one_len != split_second_len:
                correction = split_second_len
            if '+' in str(x[7]):
                variantList[count] = str(
                    split[2])+str(int(var_pos[count])+correction)+str(split[3])
            else:
                firstcomp = ''
                secondcomp = ''
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = ''.join(
                    reversed(firstcomp))+str(int(var_pos[count])+correction)+''.join(reversed(secondcomp))

    # saveDict['real_guide'] = str(realguide[0]).strip()
    saveDict['real_guide'] = str(x[15]).strip()
    saveDict['genome'] = genomeRelease
    saveDict['chr'] = str(x[4])
    saveDict['prim_pos'] = str(x[6])
    saveDict['strand'] = str(x[7])
    saveDict['highest_CFD_guide_alignment'] = str(x[1])
    saveDict['highest_CFD_alignment(alt)'] = str(x[2])
    saveDict['highest_CFD_alignment(ref)'] = str(x[3])
    saveDict['ref_seq_length'] = str(len(str(x[3])))
    saveDict['ref_pos_alt(aligned_strand)'] = ','.join(variantList)
    saveDict['pam'] = str(x[2])[-3:]
    saveDict['annotation'] = str(x[14])
    saveDict['highest_CFD_score(alt)'] = str(x[20])
    saveDict['highest_CFD_score(ref)'] = str(x[21])
    saveDict['highest_CFD_score'] = str(x[20]) if float(
        x[20]) > float(x[21]) else str(x[21])
    if 'ref' in origin:
        saveDict['highest_CFD_alignment(alt)'] = 'n'
        saveDict['highest_CFD_alignment(ref)'] = str(x[2])
        saveDict['highest_CFD_score(alt)'] = 'n'

    if str(x[8]) != str(x[30]) or str(x[9]) != str(x[31]):
        saveDict['fewest_mismatch'] = str(x[30])
        saveDict['fewest_bulge'] = str(x[31])
        saveDict['fewest_mismatch+bulge'] = str(x[32])
        saveDict['fewest_mm+bulge_guide_alignment'] = str(x[23])
        saveDict['fewest_mm+bulge_alignment(alt)'] = str(x[24])
        saveDict['fewest_mm+bulge_alignment(ref)'] = str(x[25])
        saveDict['fewest_mm+bulge_CFD_score(alt)'] = str(x[42])
        saveDict['fewest_mm+bulge_CFD_score(ref)'] = str(x[43])
        if 'ref' in origin:
            saveDict['fewest_mm+bulge_alignment(alt)'] = 'n'
            saveDict['fewest_mm+bulge_alignment(ref)'] = str(x[24])
            saveDict['fewest_mm+bulge_CFD_score(alt)'] = 'n'
    else:
        saveDict['fewest_mismatch'] = str(x[8])
        saveDict['fewest_bulge'] = str(x[9])
        saveDict['fewest_mismatch+bulge'] = str(x[10])
        saveDict['fewest_mm+bulge_guide_alignment'] = str(x[1])
        saveDict['fewest_mm+bulge_alignment(alt)'] = str(x[2])
        saveDict['fewest_mm+bulge_alignment(ref)'] = str(x[3])
        saveDict['fewest_mm+bulge_CFD_score(alt)'] = str(x[20])
        saveDict['fewest_mm+bulge_CFD_score(ref)'] = str(x[21])
        if 'ref' in origin:
            saveDict['fewest_mm+bulge_alignment(alt)'] = 'n'
            saveDict['fewest_mm+bulge_alignment(ref)'] = str(x[2])
            saveDict['fewest_mm+bulge_CFD_score(alt)'] = 'n'

    saveDict['risk_score'] = str(float(x[20])-float(x[21]))
    saveDict['absolute_risk_score'] = str(abs(float(x[20])-float(x[21])))
    saveDict['alt_haplotypes'] = str(x[19].strip().split('.')[0])
    saveDict['prim_origin'] = origin
    saveDict['prim_AF'] = str(x[17])
    saveDict['prim_samples'] = str(x[13])
    saveDict['prim_SNP_ID(positive_strand)'] = str(x[18])
    saveDict['highest_CFD_mismatch'] = str(x[8])
    saveDict['highest_CFD_bulge'] = str(x[9])
    saveDict['highest_CFD_mismatch+bulge'] = str(x[10])

    if saveDict['highest_CFD_score(ref)'] == saveDict['highest_CFD_score(alt)'] and origin == 'alt':
        saveDict['highest_CFD_alignment(alt)'] = 'n'
        saveDict['highest_CFD_score(alt)'] = 'n'
        saveDict['highest_CFD_score(ref)'] = str(x[21])
        saveDict['ref_seq_length'] = 1
        origin = 'ref'
        saveDict['prim_origin'] = origin
        saveDict['ref_pos_alt(aligned_strand)'] = 'n'
        saveDict['prim_AF'] = 'n'
        saveDict['prim_samples'] = 'n'
        saveDict['prim_SNP_ID(positive_strand)'] = 'n'
        saveDict['pam'] = str(x[3])[-3:]
        saveDict['fewest_mm+bulge_alignment(alt)'] = 'n'
        saveDict['fewest_mm+bulge_alignment(ref)'] = str(x[2])
        saveDict['fewest_mm+bulge_CFD_score(alt)'] = 'n'

    foundEmpirical = sorted(empiricalTree[int(x[6])-4:int(x[6])+4])

    for found in range(0, len(foundEmpirical)):
        empirical = foundEmpirical[found].data
        if str(saveDict['chr']) == str(empirical[0]):
            empiricalList.append(empirical[7])
            valueDict[str(empirical[4])] = empirical[5]
            empiricalDict[str(empirical[4])] = int(empirical[3])

    for key in empiricalDict:
        if int(empiricalDict[key]) < 50:
            saveDict[key] = str(valueDict[key])
            newkey = str(key)+'_mm+bul'
            saveDict[newkey] = empiricalDict[key]
            if int(empiricalDict[key]) < lowestEmpirical:
                saveDict['lowest_empirical'] = str(empiricalDict[key])

    save = ''
    for key in saveDict:
        save += str(saveDict[key])+'\t'
    save += '\n'

    outFile.write(save)

print('CHECKING MISSING RESULTS')

notFoundFile = open(outputDir + originFileName +
                    '.empirical_not_found.tsv', 'w')

for count, line in enumerate(inEmpiricalResults):
    if count not in empiricalList:
        notFoundFile.write(line)

print("INTEGRATION COMPLETED IN: %s seconds" % (time.time() - start_time))
