#!/usr/bin/env python

import gzip
import sys
import glob
import os
import multiprocessing

vcfDIR = sys.argv[1]
inPop = open(sys.argv[2], 'r').readlines()
threads = 0
try:
    threads = int(sys.argv[3])
except:
    threads = len(os.sched_getaffinity(0))-2

def convertVCF(inVCF):
    #open VCF file
    outVCFname = inVCF.replace('bgz', 'gz')
    outVCF = gzip.open(inVCF.replace('bgz', 'gz'), 'wt')
    inVCF = gzip.open(inVCF, 'rt')

    #read samplesID from file
    header = list()
    popDict = dict()
    for pop in inPop:
        if '#' in pop:
            continue
        popDict[pop.split()[0].strip()] = '0/0'

    #read header from original VCF
    for line in inVCF:
        if '##' in line:
            header.append(line.strip())
        else:
            popheader = '\t'.join(popDict.keys())
            header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Sample Collapsed Genotype">')
            header.append(line.strip()+'\tFORMAT'+'\t'+popheader+'\n')
            break
    #insert header into new VCF
    outVCF.write('\n'.join(header))

    #read each variant into the VCF
    for line in inVCF:
        split = line.strip().split('\t')
        if split[6] != 'PASS': #skip rows with no PASS in FILTER
            continue
        info = split[7].strip().split(';')
        for pop in popDict:
            popDict[pop] = '0/0'
            for index, data in enumerate(info): #read AC for each gnomAD population and insert a fake GT for each sample
                if 'AC_'+str(pop)+'=' in data:
                    ACvalue = int(data.strip().split('=')[1])
                    if ACvalue > 0:
                        popDict[pop] = '0/1'
        split[7] = info[2]
        #write each line passing the filtering into the new VCF
        outVCF.write('\t'.join(split[:8])+'\tGT\t'+'\t'.join(popDict.values())+'\n')

    inVCF.close()
    outVCF.close()
    return outVCFname

def bcftools_merging(outVCFname):    
    # finalOutVCF = outVCFname.replace('genomes','genomes.collapsed')
    tempName = outVCFname.strip().split('.')
    tempName[-3] = tempName[-3]+'.collapsed'
    finalOutVCF = '.'.join(tempName)
    os.system(f"bcftools norm -m+ -O z -o {finalOutVCF} {outVCFname}")
    
def full_process(inVCF):
    #function to process each VCF from initial conversion to collapsing multi-variant alleles with shared reference into one entry
    outVCFname = convertVCF(inVCF) #convert VCF from gnomADv3.1 to VCF4.2 compatible with CRISRPme
    bcftools_merging(outVCFname) #merge multi-variant alleles into single entry
    os.system(f"rm -f {outVCFname}") #removed unused temp vcf file


if __name__ == '__main__':
    pool = multiprocessing.Pool(threads)
    #call convert vcf for each vcf in dir
    listVCF = glob.glob(vcfDIR+'/*.vcf.bgz')
    for inVCF in listVCF:
        pool.apply_async(full_process, args=(inVCF,))
    # wait until all threads are completed than join
    pool.close()
    pool.join()
