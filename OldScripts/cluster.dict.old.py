#test1
#script per testare la lettura del file semicommon con pandas, creare una nuova colonna total, chr, pos corretta (ovver la pos senza i bulges),
# e vedere se è fattibile in termini di tempo e memoria -> non fattibile in termini di tempo

#test2
#Leggo il file riga per riga, aggiungo le colonne e salvo la lista in una lista, che poi ordino -> ok, 1gb input -> 1.5 min e 11 gb ram
# Input 1.2 Gb -> tempo 2 min; 11Gb RAM

#test 3
#Come test2, ma ogni lista è salvata in un dizionario con key = guide, per risolvere il problema di avere più guide in uno stesso cluster. 
#Ora se in una posizione ho due guide, creo due cluster separati
# Input 2.8 Gb -> 25 gb ram
# Problema con 5 gb input -> forse servono 50 Gb ram

#sys1 è target file
#sys2 is 'addGuide' or 'no' -> only for web server, only for search with only ref
#sys3 is True to keep column 5 (Pos cluster) and 9 (Total) and added guide, False to do clusterization but do not report the added columns
#sys4 is True if cluster only (no append of Total column or adding of Cluster position, because already present), False otherwise
#Output column (not written): Bulge_type, Guide, Target, chr, pos, pos_cluster (optional), direction, mms, bulge, total(optional), real guide(optional)

import time
import sys

start = time.time()
total_targets = []
guides_dict = dict()
addGuide = False                #Add real guide to last column for better grep
if 'addGuide' in sys.argv[:]:
    addGuide = True
if sys.argv[3] == 'True':
    keep_columns = True
else:
    keep_columns = False

result_name = sys.argv[1][:sys.argv[1].rfind('.')] + '.cluster.txt'
cluster_only = False
if sys.argv[4] == 'True':
    cluster_only = True

with open (sys.argv[1]) as targets:
    for line in targets:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        if not cluster_only:
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

if not cluster_only:
    print('Created \'Total\' and \'Position Cluster\' columns:', time.time() - start)
else:
    print('Loaded targets: ', time.time() - start)



start = time.time()
# total_targets.sort(key = lambda x: ( x[3] , int(x[-1]) ))
for k in guides_dict.keys():
    guides_dict[k].sort(key = lambda x: ( x[3] , int(x[5]) ))
#total_targets.sort(key = lambda x: ( x[3] , int(x[5]) ))

print('Targets sorted:', time.time() - start)

print('Start clustering')
start_time = time.time()

with open(result_name, 'w+') as result:
    if 'total' not in sys.argv[:]:      #TODO fix for offline release, praticamente se sto facendo cluster su total.txt metto l'header custom
        if addGuide:
            if keep_columns:
                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tReal Guide\n')
            else:
                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\tReal Guide\n')
        else:
            if keep_columns:
                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\n')
            else:
                result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge_Size\n')
    else:
        result.write('#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tPam_disr\tPAM_gen\tVar_uniq\n')

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
    if addGuide:
        if keep_columns:
            for cluster in total_list:
                for target in cluster:
                    result.write('\t'.join(target) + '\t' + target[1].replace('-','') + '\n')
        else:
            for cluster in total_list:
                for target in cluster:
                    result.write('\t'.join(target[0:5] + target[6:-1]) + '\t' + target[1].replace('-','') + '\n')

    else:
        if keep_columns:
            for cluster in total_list:
                for target in cluster:
                    result.write('\t'.join(target) + '\n')
        else:
            for cluster in total_list:
                for target in cluster:
                    result.write('\t'.join(target[0:5] + target[6:-1]) + '\n')


print("Clustering runtime: %s seconds" % (time.time() - start_time))