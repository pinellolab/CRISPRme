#!/usr/bin/env python
"""
Created on Sun May 30 15:49:39 2021

@author: franc
"""

"""
Script used to replace n and .  with NA in bestMerge
"""

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd
import sys
import os

in_path = sys.argv[1]

chunksize_ = 50000
chunks = pd.read_csv(in_path, sep="\t", chunksize=chunksize_)


header = True
for chunk in chunks:

    chunk = chunk.replace("n", "NA")
    # chunk = chunk.replace(regex=['\*.,\*', '\*,.\*'], value='NA')
    chunk["rsID"] = chunk["rsID"].str.replace(".", "NA")

    chunk.to_csv(in_path + ".tmp", header=header, mode="a", sep="\t", index=False)

    header = False

os.system(f"mv {in_path}.tmp {in_path}")
