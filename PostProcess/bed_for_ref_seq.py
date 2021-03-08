#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 15:48:44 2020

@author: whitebreeze
"""
import sys

targets = sys.argv[1]
out = sys.argv[2]

with open(targets, 'r') as targets:
    with open(out+'.bed_for_ref', 'w') as out:
        targets.readline()
        for line in targets:
            splitted = line.strip().split('\t')
            end_pos = int(splitted[4]) + len(splitted[1])
            out.write(splitted[3] +'\t'+ splitted[4] +'\t'+ str(end_pos) +'\n')