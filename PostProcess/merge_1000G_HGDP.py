#!/usr/bin/env python
import sys
import time
from collections import Counter

def init_summary_by_sample(path_samplesID):
    dict_samples = {}
    dict_pop = {}
    dict_superpop = {}
    with open(path_samplesID, "r") as sid:
        sid.readline()
        for line in sid:
            splitted = line.strip().split("\t")
            dict_samples[splitted[0]] = [splitted[3],
                                         splitted[1], splitted[2], 0, 0, 0, 0]
            if splitted[1] not in dict_pop.keys():
                dict_pop[splitted[1]] = 0
            if splitted[2] not in dict_superpop.keys():
                dict_superpop[splitted[2]] = 0

    return dict_samples, dict_pop, dict_superpop


def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    list_1000G = []
    list_HGDP = []
    
    for ele in cluster:
        if ele[-1] == 'HGDP':
            list_HGDP.append(ele)
        else:
            list_1000G.append(ele)

    list_HGDP.sort(key = lambda x : (float(x[cfd]), -int(x[total])), reverse=True) #taaaac
    list_1000G.sort(key = lambda x : (float(x[cfd]), -int(x[total])), reverse=True)

    best_HGDP = []
    best_1000G = []
    best_cfd_1000G_or_ref = 0
    if len(list_HGDP) > 0:
        best_HGDP = list_HGDP[0]
        best_cfd_1000G_or_ref = float(best_HGDP[cfd+1])
    if len(list_1000G) > 0:
        best_1000G = list_1000G[0]
        if best_cfd_1000G_or_ref < float(best_1000G[cfd]):
            best_cfd_1000G_or_ref = float(best_1000G[cfd])
        if best_cfd_1000G_or_ref < float(best_1000G[cfd+1]):
            best_cfd_1000G_or_ref = float(best_1000G[cfd+1])

    #append diff_cfd from best 1000G in cluster
    for ele in list_HGDP:
        ele.append(str(float(ele[cfd]) - best_cfd_1000G_or_ref))
    
    for ele in list_1000G:
        ele.append(str(float(ele[cfd]) - best_cfd_1000G_or_ref))


    #write best target in bestFile
    populations = []
    n_seq_in_cluster = str(len(list_1000G) + len(list_HGDP) - 1)
    if len(best_HGDP) > 0 and float(best_HGDP[cfd]) - best_cfd_1000G_or_ref > 0:
        best_HGDP[cfd-1] = n_seq_in_cluster
        best_HGDP[2*cfd+1] = n_seq_in_cluster
        
        for sample in best_HGDP[true_guide-2].split(','):
            populations.append(dict_samples_HGDP[sample][2])
        pops = Counter(populations)
        string_populations = []
        for pop in pops:
            string_populations.append(str(pops[pop])+"_"+str(pop))
        best_HGDP.append(','.join(string_populations))
        
        fileOut.write('\t'.join(best_HGDP)+"\n")
        list_HGDP.pop(0)
    elif len(best_1000G) > 0:
        best_1000G[cfd-1] = n_seq_in_cluster
        best_1000G[2*cfd+1] = n_seq_in_cluster

        for sample in best_1000G[true_guide-2].split(','):
            populations.append(dict_samples_1000G[sample][2])
        pops = Counter(populations)
        string_populations = []
        for pop in pops:
            string_populations.append(str(pops[pop])+"_"+str(pop))
        best_1000G.append(','.join(string_populations))

        fileOut.write('\t'.join(best_1000G)+"\n")
        list_1000G.pop(0)

    list_1000G.extend(list_HGDP)
    list_1000G.sort(key = lambda x : (float(x[cfd]), -int(x[total])), reverse=True)

    for ele in list_1000G:
        populations = []
        if ele[-2] == "HGDP":
            for sample in ele[true_guide-2].split(','):
                populations.append(dict_samples_HGDP[sample][2])
        else:
            for sample in ele[true_guide-2].split(','):
                populations.append(dict_samples_1000G[sample][2])
        pops = Counter(populations)
        string_populations = []
        for pop in pops:
            string_populations.append(str(pops[pop])+"_"+str(pop))
        ele[cfd-1] = n_seq_in_cluster
        ele[2*cfd+1] = n_seq_in_cluster
        ele.append(','.join(string_populations))
        fileOut_disc.write('\t'.join(ele)+"\n")


tau = int(sys.argv[3])
chrom = int(sys.argv[4])-1 
pos = int(sys.argv[5])-1 
total = int(sys.argv[6])-1
true_guide = int(sys.argv[7])-1
snp_info = int(sys.argv[8])-1 
cfd = int(sys.argv[9])-1
# -1 is to get the correct "python enumeration" from the bash script

samples1000G = sys.argv[10]
samplesHGDP = sys.argv[11]

dict_samples_1000G, dict_pop_1000G, dict_superpop_1000G = init_summary_by_sample(samples1000G)
dict_samples_HGDP, dict_pop_HGDP, dict_superpop_HGDP = init_summary_by_sample(samplesHGDP)

start = time.time()
with open(sys.argv[1], 'r') as fileIn:
    #header = fileIn.readline()
    with open(sys.argv[2], 'w') as fileOut:
        with open(sys.argv[2]+'.discarded_samples', 'w') as fileOut_disc:
            #fileOut.write(header)
            #fileOut_disc.write(header)
            prev_pos = -(tau+1)
            best_row = ""
            prev_guide = ""
            prev_chr = ""
            prev_snp = ""
            cluster = []
            for line in fileIn:
                splitted = line.strip().split("\t")
                if prev_guide != splitted[true_guide] or prev_chr != splitted[chrom] or int(splitted[pos]) - prev_pos > tau:
                    get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)
                    cluster = [splitted]
                else:
                    cluster.append(splitted)
                prev_guide = splitted[true_guide]
                prev_pos = int(splitted[pos])
                prev_chr = splitted[chrom]
                prev_snp = splitted[snp_info]
            
            get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Mergin done in: "+str(time.time()-start))
