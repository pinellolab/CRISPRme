#!/usr/bin/env python

import time
from intervaltree import Interval, IntervalTree
import sys
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

file_final_results = sys.argv[1]
inAnnotationFile = sys.argv[2]
file_annotated = sys.argv[3]

print("Starting annotation")
start_time = time.time()
annotationsTree = IntervalTree()
annotationsSet = set()
# guidesSet = set()       #NOTE/BUG if guide finds 0 targets, it will not be annotated

with open(inAnnotationFile, 'r') as annotations:
    for line in annotations:
        x = line.split('\t')
        annotations_list = str(x[3]).strip()
        annotationsTree[int(x[1]):int(x[2])] = str(x[0])+'\t'+annotations_list
        annotationsSet.add(annotations_list)


with open(file_final_results, 'r') as f_in:
    with open(file_annotated, 'w') as f_out:
        header = f_in.readline()
        f_out.write(header)
        for line in f_in:
            splitted = line.rstrip().split('\t')
            guide_no_bulge = splitted[1].replace('-', '')
            # Calcolo annotazioni
            foundAnnotations = sorted(
                annotationsTree[int(splitted[5]):(int(splitted[5])+int(len(guide_no_bulge))+1)])
            string_annotation = []
            found_bool = False
            for found in range(0, len(foundAnnotations)):
                guide = foundAnnotations[found].data
                guideSplit = guide.split('\t')
                if(str(guideSplit[0]) == str(splitted[4])):
                    found_bool = True
                    string_annotation.append(str(guideSplit[1]))
            if not found_bool:
                last_annotation = 'n'
            else:
                last_annotation = ','.join(list(set(string_annotation)))
            splitted[14] = last_annotation  # bestCFD
            # splitted[36] = last_annotation #fewestMM_BUL
            # splitted[58] = last_annotation #bestCRISTA

            f_out.write('\t'.join(splitted)+'\n')

print(f"Annotation done in {time.time()-start_time}")
