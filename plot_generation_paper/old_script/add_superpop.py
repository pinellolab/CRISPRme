import sys


inSample = open(sys.argv[1],'r').readlines()
inTarget = open(sys.argv[2],'r')


sampleDict = dict()
for line in inSample:
    split = line.strip().split('\t')
    superpop = split[2]
    sampleDict[split[0]] = superpop
    
    
for line in inTarget:
    split = line.strip().split('\t')
    samples = split[13].strip().split(',')
    superpop = set()
    for sample in samples:
        superpop.add(sampleDict[sample])
    superpop = list(superpop)
    print(line.strip()+'\t'+','.join(superpop))
