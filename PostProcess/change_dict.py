#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 01:03:59 2020

@author: francesco
"""

import gzip
import sys
import json
import time
import bisect
import os


def findOccurrences(s, ch):
    """
    Finds the positions of a character ch in a string s.
    """
    return [i for i, letter in enumerate(s) if letter == ch]


def findVars(s):
    """
    Find all the variations in a dictionary entry s.
    Creates a dicitonary result with variations as keys and position in the entry s 
    as item associated to the key
    """
    posVar = findOccurrences(s, ";")
    length = len(posVar)
    result = {}
    if length > 1:  # if there are more than one variant for a certain chr position
        for i in range(0, length-2, 2):
            # get a single variation at a time
            result[s[posVar[i]+1:posVar[i+1]]] = posVar[i]
        result[s[posVar[length-1]+1:]] = posVar[length-1]  # get last variation
        return result
    else:
        result[s[posVar[0]+1:]] = posVar[0]  # get the only variation
        return result


def getSamples(entry, variations, target):
    """
    Extract sample list for the chosen variation and return it along 
    its starting and ending positions.
    """
    keys = list(variations.keys())
    if len(keys) > 1:  # if there are more than one variant for a certain chr position
        if target == keys[0]:
            # get first variation and its position
            return entry[:variations[target]], 0, variations[keys[0]]
        elif target == keys[len(keys)-1]:
            # get last variation and its position
            return entry[variations[keys[len(keys)-2]]+5:variations[target]], variations[keys[len(keys)-2]]+5, variations[target]
        else:
            # get in between variation and its position
            return entry[variations[keys[keys.index(target)-1]]+5:variations[target]], variations[keys[keys.index(target)-1]]+5, variations[target]
    else:
        # get the only variation and its position
        return entry[:variations[keys[0]]], 0, variations[keys[0]]


def putIntoEntrySorted(entry, toInsert):
    """
    Insert into the ordered entry the new sample maintaining the order (if not already present)
    """
    index = bisect.bisect(entry, toInsert)
    inserted = False
    if entry[index-1] != toInsert:
        entry.insert(index, toInsert)
        inserted = True
    return entry, inserted


def putIntoEntry(entry, toInsert):
    """
    Insert into the unordered entry the new sample (if not already present)
    """
    inserted = False
    if not(toInsert in entry):
        entry.append(toInsert)
        inserted = True
    return entry, inserted


def updateDictionary(oldDictionaryFile, newVCFFile):
    '''
    oldDictionaryFile : dictionary .json that will be updated
    newVCFFile: .vcf.gz file containing the new entries to be added to the old .json
    '''
    isSorted = True
    verbose = False
    oldDict = None
    oldDict = json.load(open(oldDictionaryFile))

    chr_dict = dict()
    start_time = time.time()
    oldEntry = 0
    newEntry = 0
    with gzip.open(newVCFFile, 'rb') as targets:
        if verbose:
            logFile = open("logUpdateDictionaries/log__"+os.path.basename(
                oldDictionaryFile)+"__"+os.path.basename(newVCFFile)+".log", 'w')
        # Skip vcf header
        for line in targets:
            line = line.decode('ascii')
            if ('#CHROM') in line:
                # Save this header for retrieving sample id
                column_vcf = line.strip().split('\t')
                break
        # Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        for line in targets:
            line = line.decode('ascii').strip().split('\t')
            list_samples = []
            list_chars = []
            # if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
            for pos, i in enumerate(line[9:]):
                if ('1' in i):
                    list_samples.append(column_vcf[pos + 9])
            if "chr" not in line[0]:
                line[0] = "chr"+str(line[0])
            chr_pos_string = line[0] + ',' + line[1]
            # Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
            # If no sample was found, the dict is chrX,100 -> ;A,T
            list_chars.append(line[3])
            list_chars.append(line[4])

            if chr_pos_string in oldDict:  # entry already present
                oldEntry += 1
                variations = findVars(oldDict[chr_pos_string])
                varToFind = list_chars[0]+","+list_chars[1]
                if varToFind in variations:  # variant already present
                    entry, start, end = getSamples(
                        oldDict[chr_pos_string], variations, varToFind)
                    entry = entry.split(",")
                    goodNewSamples = 0
                    for sample in list_samples:
                        if isSorted:
                            entry, inserted = putIntoEntrySorted(entry, sample)
                        else:
                            entry, inserted = putIntoEntry(entry, sample)
                        if inserted:
                            goodNewSamples += 1
                    if goodNewSamples > 0:
                        if verbose:
                            print("New Sample(s) in "+chr_pos_string +
                                  ": #"+str(goodNewSamples)+" new", file=logFile)
                        oldDict[chr_pos_string] = oldDict[chr_pos_string][:start] + \
                            ','.join(entry) + oldDict[chr_pos_string][end:]
                else:  # variant not present
                    if verbose:
                        print("New Variation in "+chr_pos_string+": "+varToFind +
                              " with #"+str(len(list_samples))+" samples", file=logFile)
                    oldDict[chr_pos_string] = oldDict[chr_pos_string] + \
                        "/" + ','.join(list_samples) + ';' + \
                        ','.join(list_chars)
            else:  # entry not present
                newEntry += 1
                if verbose:
                    print("New Entry "+chr_pos_string+" with " +
                          str(len(list_samples))+" samples", file=logFile)
                try:
                    oldDict[chr_pos_string] = ','.join(
                        list_samples) + ';' + ','.join(list_chars)
                except:
                    oldDict[chr_pos_string] = ';' + ','.join(list_chars)
        if verbose:
            logFile.close()

    with open(os.path.dirname(oldDictionaryFile) + "/my_dict_" + str(line[0]) + ".json", 'w') as f:
        json.dump(oldDict, f)
    print('Updated ' + oldDictionaryFile + ' in', time.time() - start_time)
    print("Old entries updated: "+str(oldEntry),
          "New entries inserted: "+str(newEntry))
