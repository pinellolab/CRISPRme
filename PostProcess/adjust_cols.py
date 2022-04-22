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

inFile = sys.argv[1]

start = time.time()
df = pd.read_csv(inFile, sep='\t')
cols = df.columns.tolist()
cols.remove('CFD')
cols.remove('CFD_ref')
cols.remove('Reference')
try:
    cols.remove('Chromosome_fake')
except:
    pass
good_cols = cols[0:3] + ['Reference'] + cols[3:] + ['CFD', 'CFD_ref']
df = df[good_cols]
df.to_csv(inFile+'.tmp', index=False, sep='\t')
os.system('mv '+inFile+'.tmp'+" "+inFile)

print('Results adjusted in: '+str(time.time()-start))