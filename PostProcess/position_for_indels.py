#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:31:45 2020

@author: francesco
"""
import pandas as pd
import numpy as np
import time
import sys

print("Constructing file for indel extraction from reference fasta file")

start_time = time.time()
indels_file = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(indels_file, sep="\t")
# df_out = pd.DataFrame(np.zeros([df.shape[0], 3]), columns=["CHR", "START", "END"], dtype=np.int)

with open(out + ".pos_indels", "w") as out_file:

    chro = str(df.iloc[0, 0])
    for line in range(df.shape[0]):
        true_start_pos, ref = df.iloc[line, 1:3]
        start_pos = true_start_pos - 26
        len_ref = len(ref)
        end_pos = true_start_pos + len_ref + 25
        # df_out.iloc[line,:] = ["chr"+chro, start_pos, end_pos]
        out_file.write(
            "chr" + chro + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n"
        )

# df_out.to_csv(out+".pos_indels", sep='\t', index=False, header=None)

print("Done in " + str(time.time() - start_time))
