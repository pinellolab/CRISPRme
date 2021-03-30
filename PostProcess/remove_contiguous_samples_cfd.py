#!/usr/bin/env python
"""
Created on Fri Aug 28 15:58:04 2020

@author: francesco
"""
import sys
import time

'''
def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    list_ref = []
    dict_var = dict()
    for ele in cluster:
        if ele[snp_info] == 'n':
            list_ref.append(ele)
        else:
            if ele[snp_info] in dict_var.keys():
                dict_var[ele[snp_info]].append(ele)
            else:
                dict_var[ele[snp_info]] = [ele]
    
    list_ref.sort(key = lambda x : x[total])
    if len(list_ref) > 1:
        best_ref = list_ref[0]
        for ele_ref in list_ref[1:]:
            if float(ele_ref[cfd]) > float(best_ref[cfd]):
                fileOut_disc.write("\t".join(best_ref))
                best_ref = ele_ref
            else:
                fileOut_disc.write("\t".join(ele_ref))
        fileOut.write("\t".join(best_ref))
    elif len(list_ref) == 1:
        fileOut.write("\t".join(list_ref[0]))
    
    for key in dict_var.keys():
        list_var = dict_var[key]
        list_var.sort(key = lambda x : x[total])
        best_var = list_var[0]
        if len(list_var) > 1:
            for ele_var in list_var[1:]:
                if float(ele_var[cfd]) > float(best_var[cfd]):
                    fileOut_disc.write("\t".join(best_var))
                    best_var = ele_var
                else:
                    fileOut_disc.write("\t".join(ele_var))
            fileOut.write("\t".join(best_var))
        else:
            fileOut.write("\t".join(best_var))
'''   
def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    final_list = []
    list_ref = []
    dict_var = dict()
    for ele in cluster:
        if ele[snp_info] == 'n':
            list_ref.append(ele)
        else:
            if ele[snp_info] in dict_var.keys():
                dict_var[ele[snp_info]].append(ele)
            else:
                dict_var[ele[snp_info]] = [ele]
    
    list_ref.sort(key = lambda x : int(x[total]))
    if len(list_ref) > 1:
        best_ref = list_ref[0]
        for ele_ref in list_ref[1:]:
            if float(ele_ref[cfd]) > float(best_ref[cfd]):
                best_ref = ele_ref
        final_list.append(best_ref)
    elif len(list_ref) == 1:
        final_list.append(list_ref[0])
    
    for key in dict_var.keys():
        list_var = dict_var[key]
        list_var.sort(key = lambda x : int(x[total]))
        best_var = list_var[0]
        if len(list_var) > 1:
            for ele_var in list_var[1:]:
                if float(ele_var[cfd]) > float(best_var[cfd]):
                    best_var = ele_var
            final_list.append(best_var)
        else:
            final_list.append(best_var)
    
    n_ele = len(final_list)
    if n_ele > 1:
        final_list.sort(key = lambda x : float(x[cfd]), reverse=True)
        if final_list[0][cfd] == final_list[1][cfd] and final_list[1][cfd-2] == 'n':
            final_list[1][cfd-1] = str(n_ele-1)
            final_list[1][2*cfd+1] = str(n_ele-1)
            fileOut.write("\t".join(final_list[1]))
            final_list.pop(1)
        else:
            final_list[0][cfd-1] = str(n_ele-1)
            final_list[0][2*cfd+1] = str(n_ele-1)
            fileOut.write("\t".join(final_list[0]))
            final_list.pop(0)
        for ele in final_list:
            ele[cfd-1] = str(n_ele-1)
            ele[2*cfd+1] = str(n_ele-1)
            fileOut_disc.write("\t".join(ele))
    elif n_ele == 1:
        fileOut.write("\t".join(final_list[0]))

tau = int(sys.argv[3])
chrom = int(sys.argv[4])-1 
pos = int(sys.argv[5])-1 
total = int(sys.argv[6])-1
true_guide = int(sys.argv[7])-1
snp_info = int(sys.argv[8])-1 
cfd = int(sys.argv[9])-1
# -1 is to get the correct "python enumeration" from the bash script

start = time.time()
with open(sys.argv[1], 'r') as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], 'w') as fileOut:
        with open(sys.argv[2]+'.discarded_samples', 'w') as fileOut_disc:
            fileOut.write(header)
            fileOut_disc.write(header)
            prev_pos = -(tau+1)
            best_row = ""
            prev_guide = ""
            prev_chr = ""
            prev_snp = ""
            cluster = []
            for line in fileIn:
                splitted = line.split("\t")
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
