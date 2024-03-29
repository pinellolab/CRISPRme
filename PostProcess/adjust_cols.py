#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:40:01 2020

@author: francesco
"""

import pandas as pd 
import sys
import time
import os


# start = time.time()
# df = pd.read_csv(inFile, sep='\t')
# cols = df.columns.tolist()
# cols.remove('CFD')
# cols.remove('CFD_ref')
# cols.remove('Reference')
# try:
#     cols.remove('Chromosome_fake')
# except:
#     pass
# good_cols = cols[0:3] + ['Reference'] + cols[3:] + ['CFD', 'CFD_ref']
# df = df[good_cols]
# df.to_csv(inFile+'.tmp', index=False, sep='\t')
# os.system('mv '+inFile+'.tmp'+" "+inFile)

# print('Results adjusted in: '+str(time.time()-start))

#path to infile
inFile = sys.argv[1]
#start tracking time
start = time.time()

#read file to adjust in chunks
chunksize_ = 100000
chunks = pd.read_csv(inFile, sep = '\t', chunksize=chunksize_)

#write header the first time a chuck is processed
header = True
for chunk in chunks:
    #extract columns and then remove unwanted ones
    cols = chunk.columns.tolist()
    cols.remove('CFD')
    cols.remove('CFD_ref')
    cols.remove('Reference')
    try:
        cols.remove('Chromosome_fake')
    except:
        pass
    #reoreder cols in wanted order
    good_cols = cols[0:3] + ['Reference'] + cols[3:] + ['CFD', 'CFD_ref']
    chunk = chunk[good_cols]    
    chunk.to_csv(inFile+'.tmp', header = header, mode = 'a', sep = '\t', index=False) 
    #stop writing header
    header = False

#change name to tmp file
os.system('mv '+inFile+'.tmp'+" "+inFile)
print('Results adjusted in: '+str(time.time()-start))
