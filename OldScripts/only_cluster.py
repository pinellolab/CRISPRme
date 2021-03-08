#Script to cluster. Add columns total and cluster pos only if sys2 is false

#sys1 è semicommon file
#sys2 is True if total and cluster pos are already in the file, else is False
#Output column (not written): Bulge_type, Guide, Target, chr, pos, pos_cluster, direction, mms, bulge, total
import pandas as pd
import time
import sys

start = time.time()
total_targets = []
guides_dict = dict()
result_name = sys.argv[1][:sys.argv[1].rfind('.')] + '.cluster.txt'
cluster_pos = False
if sys.argv[2] == 'True':
    cluster_pos = True

with open (sys.argv[1]) as targets:
    for line in targets:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        if not cluster_pos:
            line.append(str(int(line[6]) + int(line[7])))
            if line[5] == '+':
                if line[0] == 'DNA':
                    # line.append(str(int(line[4]) + int(line[7])))
                    line.insert(5, str(int(line[4]) + int(line[7])))
                else:
                    # line.append(str(int(line[4]) - int(line[7])))
                    line.insert(5,str(int(line[4]) - int(line[7])))
            else:
                # line.append(line[4])
                line.insert(5, line[4])
        try:
            guides_dict[line[1].replace('-','')].append(line)
        except:
            guides_dict[line[1].replace('-','')] = [line]
        #total_targets.append(line)
        
print('Created set of targets for all guides:', time.time() - start)
start = time.time()
# total_targets.sort(key = lambda x: ( x[3] , int(x[-1]) ))
for k in guides_dict.keys():
    guides_dict[k].sort(key = lambda x: ( x[3] , int(x[5]) ))
#total_targets.sort(key = lambda x: ( x[3] , int(x[5]) ))

print('Targets sorted:', time.time() - start)

print('Start clustering')
start_time = time.time()

with open(result_name, 'w+') as result:
    total_targets = []
    for k in guides_dict.keys():
        total_targets += guides_dict[k]
    total_list = []

    first_line = total_targets[0]
    # current_chr_pos = first_line[3] + ' ' + first_line[9]
    current_chr_pos = first_line[3] + ' ' + first_line[5]

    total_list.append([first_line])

    for line in total_targets[1:]:
        #if line[3] + ' ' + line[9] != current_chr_pos:
        if line[3] + ' ' + line[5] != current_chr_pos:
            # total_list[-1].sort(key = lambda x: int(x[8]))
            total_list[-1].sort(key = lambda x: int(x[9]))
            total_list.append([line])
            # current_chr_pos = line[3] + ' ' + line[9]
            current_chr_pos = line[3] + ' ' + line[5]
        else:
            total_list[-1].append(line)     

    total_list[-1].sort(key = lambda x: int(x[9]))

    total_list.sort(key = lambda x: int(x[0][9]))
    for cluster in total_list:
        for target in cluster:
            result.write('\t'.join(target) + '\n')

print("Clustering runtime: %s seconds" % (time.time() - start_time))