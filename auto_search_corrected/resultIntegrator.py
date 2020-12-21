#!/usr/bin/env python

from intervaltree import Interval, IntervalTree
import sys
import time
import concurrent.futures
import glob


def rev_comp(a):
    # perform reverse complement on a string
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'


print('READING INPUT FILES')

# READ INPUT FILES
crispritzResultFile = sys.argv[1]  # file with CRISPRitz results
empiricalResults = sys.argv[2]  # file with empirical results
annotationFile = sys.argv[3]  # annotation file used to find nearest gene
guideFile = sys.argv[4]  # real guide used in the search
outputDir = sys.argv[5]  # output file name
check = sys.argv[6].upper()  # output file name
genomeRelease = str(sys.argv[7])  # genome used in the search phase

# OPEN INPUT FILES AND PREPARE OUTPUT FILE
# crispritz results file open
inCrispritzResults = open(crispritzResultFile, 'r')
# empirical-seq data open
inEmpiricalResults = open(empiricalResults, 'r').readlines()
# annotation data open
inAnnotationFile = open(annotationFile, 'r')
# guide file
realguide = open(guideFile, 'r').readlines()

empiricalTree = IntervalTree()
empiricalList = []
genomeDict = {}
empiricalDict = {}
valueDict = {}
saveDict = {"real_guide": 'n', "genome": 'n', "chr": 'n', "prim_pos": 'n', "strand": 'n', "highest_CFD_guide_alignment": 'n', "highest_CFD_alignment(ref)": 'n',
            "highest_CFD_alignment(alt)": 'n', "ref_seq_length": 'n', "ref_pos_alt(aligned_strand)": 'n', "pam": 'n', "annotation": 'n', "highest_CFD_score": 'n',
            "highest_CFD_score(ref)": 'n', "highest_CFD_score(alt)": 'n', "risk_score": 'n', "absolute_risk_score": 'n', "highest_CFD_mismatch": 'n',
            "highest_CFD_bulge": 'n', "highest_CFD_mismatch+bulge": 'n', "fewest_mm+bulge_guide_alignment": 'n', "fewest_mm+bulge_alignment(ref)": 'n',
            "fewest_mm+bulge_alignment(alt)": 'n', "fewest_mm+bulge_CFD_score(ref)": 'n', "fewest_mm+bulge_CFD_score(alt)": 'n', "fewest_mismatch": 'n',
            "fewest_bulge": 'n', "fewest_mismatch+bulge": 'n', "alt_haplotypes": 'n', "prim_origin": 'n', "prim_AF": 'n', "prim_samples": 'n',
            "prim_SNP_ID(positive_strand)": 'n', "gene_name": 'n', "gene_ID": 'n', "gene_annotation": 'n', "gene_distance(kb)": 'n', "lowest_empirical": 'n',
            "Nature2019": 'n', "Nature2019_mm+bul": 'n', "CHANGEseq": 'n', "CHANGEseq_mm+bul": 'n', "CIRCLEseq": 'n', "CIRCLEseq_mm+bul": 'n', "ONEseq": 'n',
            "ONEseq_mm+bul": 'n', "GUIDEseq_293": 'n', "GUIDEseq_293_mm+bul": 'n', "GUIDEseq_CD34": 'n', "GUIDEseq_CD34_mm+bul": 'n', "GUIDEseq": 'n',
            "GUIDEseq_mm+bul": 'n'}

start_time = time.time()

print('CREATING INTERVAL TREES')

for count, line in enumerate(inEmpiricalResults):
    empList = line.strip().split('\t')
    empList = [elem.strip() for elem in empList]
    empList.append(count)
    # adding empirical data to the tree
    empiricalTree[int(empList[2]):int(empList[3])] = empList
    empiricalDict[str(empList[5])] = 50
    valueDict[str(empList[5])] = 'n'

# open outputDir to write results
originFileName = crispritzResultFile.split('/')
originFileName = originFileName[len(originFileName)-1]
outFile = open(outputDir + originFileName +
               '.integrated_results.tsv', 'w')

# writing header in file
save = '#'
for key in saveDict:
    save += str(key)+'\t'
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

    saveDict['real_guide'] = str(realguide[0]).strip()
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
            valueDict[str(empirical[5])] = empirical[6]
            empiricalDict[str(empirical[5])] = int(empirical[4])

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
