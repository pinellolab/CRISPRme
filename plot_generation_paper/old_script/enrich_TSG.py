import sys

in_targets = open(sys.argv[1], 'r')
in_human_TSG = open(sys.argv[2], 'r')
header = in_targets.readline().strip()+"\tGene_description\tGene_type"

geneDict = dict()
for line in in_human_TSG:
    split = line.strip().split('\t')
    # inserisco description e function in dict con gene come key
    #genesymbol = [description,genetype]
    geneDict[split[1]] = [split[6], split[7]]

out_targets = open(sys.argv[1]+'_annotated.txt', 'w')
out_targets.write(header+'\n')
for line in in_targets:
    split = line.strip().split('\t')
    new_cols = ['NA', 'NA']
    try:
        new_cols = geneDict[split[47]]
    except:
        pass
    split.extend(new_cols)
    out_targets.write('\t'.join(split)+'\n')
