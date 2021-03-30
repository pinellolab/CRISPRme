'''
Script reads the top_1.sample file, for each cluster create a entry in dictionary
guide,chr,clusterpos -> union of all samples in this cluster
The opens the cluster targets file and for each top1 of the cluster assign the corresponding 
dictionary entry. For the other targets in the cluster writes 'n' -> TODO pensare a cosa scrivere
'''

#TODO migliorarlo per eventuale file top1_sample con dimensioni molto grandi
# sys1 is targets.cluster file
# sys2 is top1.samples file
# sys3 is job_id
import sys

sample_dict = dict()
current_pos = '0'
with open(sys.argv[2], 'r') as top_samples:
    for line in top_samples:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        try:
            sample_dict[line[1].replace('-','') + line[3] + line[5]] = sample_dict[line[1].replace('-','') + line[3] + line[5]].union(set(line[-2].split(',')))
        except:
            sample_dict[line[1].replace('-','') + line[3] + line[5]] = set(line[-2].split(','))
       
with open(sys.argv[1]) as targets, open(sys.argv[3] + '.final.txt', 'w+') as result:
    for line in targets:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        if current_pos != line[1].replace('-','') + line[3] + line[5]:
            try:
                line.append(','.join(sample_dict[line[1].replace('-','') + line[3] + line[5]]))
            except:
                line.append('n')
            current_pos = line[1].replace('-','') + line[3] + line[5]
        else:
            line.append('n')
        result.write('\t'.join(line) + '\t' + line[1].replace('-','') + '\n') #line[1].replace('-','') added to do better grep

