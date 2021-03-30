#!/usr/bin/env python

import gzip
import sys
import re
inVCF = gzip.open(sys.argv[1], 'rt')
inPop = open(sys.argv[2], 'r').readlines()
outVCF = gzip.open(sys.argv[1].replace('bgz', 'gz'), 'wt')

header = list()
popDict = dict()
for pop in inPop:
    popDict[pop.strip()] = '0|0'

for line in inVCF:
    if '##' in line:
        header.append(line.strip())
    else:
        popheader = '\t'.join(popDict.keys())
        header.append(line.strip()+'\t'+popheader+'\n')
        break

# print('\n'.join(header))
outVCF.write('\n'.join(header))


for line in inVCF:
    split = line.strip().split('\t')
    info = split[7].strip().split(';')
    for pop in popDict:
        popDict[pop] = '0|0'
        for index, data in enumerate(info):
            if 'AC-'+str(pop)+'=' in data:
                ACvalue = int(data.strip().split('=')[1])
                if ACvalue > 0:
                    popDict[pop] = '0|1'
    outVCF.write('\t'.join(split[:8])+'\t'+'\t'.join(popDict.values())+'\n')

inVCF.close()
outVCF.close()
# print('\t'.join(split[:8]), '\t', '\t'.join(popDict.values()))
