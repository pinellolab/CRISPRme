#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 12:32:36 2020

@author: francesco
"""

import pandas as pd
import sys
import numpy as np
import time

start_time = time.time()
df = pd.read_csv(sys.argv[1], sep='\t')
df_out = pd.DataFrame(df.values, columns = df.columns)
df_out = df_out.rename(columns={"Chromosome":"Chromosome_fake"})
n_cols = df.shape[1]
df_out["Chromosome"] = np.zeros([df.shape[0],1])

chroms = df_out["Chromosome_fake"].copy()
for line in range(df.shape[0]):
    splitted_dash = chroms[line].split('-')
    chrom_fake, start_position = splitted_dash[0].split('_')[0], int(splitted_dash[0].split('_')[1])
    #samples = chroms[line].split('_')[-1]
    df_out.iloc[line,4] = df.iloc[line,4]+start_position
    df_out.iloc[line,5] = df.iloc[line,5]+start_position
    df_out.iloc[line,n_cols] = chrom_fake
    #df_out.iloc[line,-1] = samples
    #df_out.iloc[line,3] = splitted_dash[0].split(':')[0]

cols = df_out.columns.tolist()
cols.remove("Chromosome_fake")
cols.remove("Chromosome")
df_cols = cols[0:3] + ["Chromosome"] + cols[3:] + ["Chromosome_fake"]
df_out = df_out[df_cols]
df_out.to_csv(sys.argv[1]+".corrected", index=False, sep='\t')
print("Adjusted position in ", time.time()-start_time)

