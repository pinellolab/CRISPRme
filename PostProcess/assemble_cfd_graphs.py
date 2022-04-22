#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:55:37 2020

@author: francesco
"""
import sys
import os
import pandas as pd
import numpy as np

output_dir = sys.argv[1]

files_fake = []
files_true = []
for f in os.listdir(output_dir):
    if 'CFDGraph' in f:
        if 'fake' in f:
            files_fake.append(f)
        else:
            files_true.append(f)

os.chdir(output_dir)
if files_fake:    
    empty_m = np.zeros([101,2])
    cumulative_graph_fake = pd.DataFrame(empty_m, columns=['ref','var'])
    for f in files_fake:
        df = pd.read_csv(f, sep='\t')
        cumulative_graph_fake = cumulative_graph_fake.add(df)
        os.system(f"rm {f}")
        
    cumulative_graph_fake = cumulative_graph_fake.astype('int64')
    cumulative_graph_fake.to_csv("indels.CFDGraph.txt", sep='\t', index=False)

if files_true:    
    empty_m = np.zeros([101,2])
    cumulative_graph_true = pd.DataFrame(empty_m, columns=['ref','var'])
    for f in files_true:
        df = pd.read_csv(f, sep='\t')
        cumulative_graph_true = cumulative_graph_true.add(df)
        os.system(f"rm {f}")
        
    cumulative_graph_true = cumulative_graph_true.astype('int64')
    cumulative_graph_true.to_csv("snps.CFDGraph.txt", sep='\t', index=False)