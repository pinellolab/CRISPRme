#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 21:49:49 2020

@author: francesco
"""
import sys
import os
import time

fileCFD = sys.argv[1]

file_out = open(fileCFD+".tmp.txt","w")
start = time.time()
with open(fileCFD, "r") as f:
    file_out.write(f.readline())
    prev_key = "empty"
    cluster = []
    for line in f:
        splitted = line.rstrip().split("\t")
        key = splitted[3] + " " + splitted[5] + " " + splitted[6]
        if key == prev_key:
            cluster.append(splitted)
        else:
            if cluster:
                #cluster.sort(key = lambda x : float(x[-1]), reverse = True)
                cluster = ['\t'.join(x) for x in cluster]
                file_out.write('\t'.join(cluster)+"\n")
                prev_key = key
                cluster = [splitted]
            else:
                prev_key = key
                cluster = [splitted]
    
    #last cluster
    if len(cluster) > 1:
        #cluster.sort(key = lambda x : float(x[-1]), reverse = True)
        cluster = ['\t'.join(x) for x in cluster]
        file_out.write('\t'.join(cluster)+"\n")
    else:
        cluster = ['\t'.join(x) for x in cluster]
        file_out.write('\t'.join(cluster)+"\n")
    
    os.system("mv "+fileCFD+".tmp.txt "+fileCFD)

print('Adjusting results done in: '+str(time.time()-start))
        
            
