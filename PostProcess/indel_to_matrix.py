#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 01:04:50 2020

@author: francesco
"""

import gzip
import sys
import json
import time
import os
import pandas as pd
import numpy as np
#from biclustlib.algorithms import BitPatternBiclusteringAlgorithm
from itertools import combinations

# argv1 is ALLchrx.gzip --> Create the json file containing chr, pos, ref, alt, list of samples (HG001,HG002...)
# argv2 is chr number (eg 1)


#vcf = "/home/whitebreeze/Tirocinio_Giugno/correct_vcf/CHR15/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.indels_only_v3_matrix"#sys.argv[1]
vcf = sys.argv[1] #"/home/whitebreeze/Tirocinio_Giugno/correct_vcf/CHR5/matrix.tsv"

#chr_dict = dict()
start_time = time.time()
#out = open(vcf+".indels_only", 'w')

df = pd.read_csv(vcf, sep='\t', dtype=np.bool)
#df3 = df.T
#df3.to_csv('/home/whitebreeze/Tirocinio_Giugno/correct_vcf/t_matrix_c.csv', index=False)
#print(df3.groupby(df3.columns, axis=1).first())
#clus = BitPatternBiclusteringAlgorithm(min_rows=2, min_cols=2)
#results = clus.run(df)
print(df.shape)
indexes = [i for i in range(df.shape[1])]
row_sum  = np.sum(df.values, axis=0)
max_index = np.argmax(row_sum)
print(max_index)
max_sum = max(row_sum)
indexes.remove(max_index)
cluster = df.iloc[:,max_index]#np.zeros([df.shape[0], 1])
group = [max_index]

with open("log_indels.txt", "w") as log_file:
    log_file.write("Starting with "+str(max_sum)+" and sample "+str(df.columns[max_index])+"("+str(max_index)+").\n")
    
    it = 1
    while True:
        #print(it, len(indexes), max_sum)
        if len(indexes) == 0:
            break
        tmp_cluster = None
        tmp_index = -1
        max_sum = -1
        for confronto in indexes:
            and_confronto = cluster & df.iloc[:,confronto]
            sum_confronto = sum(and_confronto)
            if sum_confronto > max_sum:
                max_sum = sum_confronto
                tmp_cluster = and_confronto
                tmp_index = confronto
                
        indexes.remove(tmp_index)
        cluster = tmp_cluster
        group.append(tmp_index)
        log_file.write("It: "+str(it)+". New cluster has "+str(max_sum)+" indels in it. We added sample "+str(df.columns[tmp_index])+"("+str(tmp_index)+").\n")
        it += 1
        if max_sum == 0:
            print("Breaking early, 0 indels overlapping!")
            log_file.write("Breaking early, 0 indels overlapping!")
            break

print("Done in "+str(time.time()-start_time))
with open('trace.txt', 'w') as trace:
    trace.write(','.join(group))

"""
mock = pd.DataFrame(np.ones([20,20]), dtype=np.bool)
mock.iloc[0,2] = 0
mock.iloc[1,6] = 0
mock.iloc[0,3] = 0
indexes_m = [i for i in range(mock.shape[1])]
dict_samples_m = {}
it = 1
while True:
    print(it, len(indexes_m))
    if len(indexes_m) == 0:
        break
    perno = indexes_m[0]
    cluster = [perno]
    if len(indexes_m) > 1:
        for confronto in indexes_m[1:]:
            if (mock.iloc[:,perno] == mock.iloc[:,confronto]).all():
                cluster.append(confronto)
    for ele in cluster:
        indexes_m.remove(ele)
    dict_samples_m[perno] = cluster
    it += 1
"""
"""
dict_indels = {}
indexes = [i for i in range(df.shape[0])]
it = 1
while True:
    print(it, len(indexes))
    if len(indexes) == 0:
        break
    perno = indexes[0]
    cluster = [perno]
    if len(indexes) > 1:
        for confronto in indexes[1:]:
            if (df.iloc[perno,:] == df.iloc[confronto,:]).all():
                cluster.append(confronto)
    for ele in cluster:
        indexes.remove(ele)
    dict_indels[perno] = cluster
    it += 1
"""
