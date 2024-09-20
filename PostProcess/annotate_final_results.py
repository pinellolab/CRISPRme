#!/usr/bin/env python

import time
from typing import Dict, List, Set, Tuple
from intervaltree import IntervalTree
import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# create the IntervalTree only of the chromosome to be queried later
def extractGenomesFromResults(files_final_result: List[str]) -> Set[str]:
    genome_to_process = set(())
    for f in files_final_result:
        with open(f, 'r') as f_in:
            for line in f_in:
                # add the genome in the result file
                splitted = line.rstrip().split('\t')
                genome_to_process.add(str(splitted[4]))
    return genome_to_process

def loadAnnotationsFiltered(inAnnotationFile: str, genome_to_process: Set[str]) -> Dict[str, IntervalTree]:
    annotationDict = dict()
    with open(inAnnotationFile, 'r') as annotations:
        for line in annotations:
            x = line.split('\t')
            if 'vuoto.txt' in inAnnotationFile:
                break
            
            genome = str(x[0])
            if genome not in genome_to_process:
                continue
            if genome not in annotationDict.keys():
                annotationDict[genome]=IntervalTree()
            
            annotations_list = str(x[3]).strip()
            annotationDict[genome][int(x[1]):int(x[2])] = f"{str(x[0])}\t{annotations_list}"
    return annotationDict

def writeAnnotated(file_result: str, file_annotated: str, annotationDict: Dict[str, IntervalTree]):
    with open(file_result, 'r') as f_in:
        with open(file_annotated, 'w') as f_out:
            header = f_in.readline()
            f_out.write(header)
            for line in f_in:
                splitted = line.rstrip().split('\t')
                guide_no_bulge = splitted[1].replace('-', '')
                # Calcolo annotazioni
                genome = str(splitted[4])
                foundAnnotations = list()
                if genome in annotationDict.keys():
                    some_idx = int(splitted[5])
                    foundAnnotations = list(annotationDict[genome][some_idx:(some_idx+int(len(guide_no_bulge))+1)])
                # TODO if genome is not found, should it be an error?

                # sorted(annotationsTree[int(splitted[5]):(int(splitted[5])+int(len(guide_no_bulge))+1)])
                string_annotation = list()
                for found in range(0, len(foundAnnotations)):
                    guide = foundAnnotations[found].data
                    guideSplit = guide.split('\t')
                    string_annotation.append(str(guideSplit[1]))
                if len(string_annotation) == 0:
                    last_annotation = 'n'
                else:
                    last_annotation = ','.join(list(set(string_annotation)))
                splitted[14] = last_annotation  # bestCFD
                # splitted[36] = last_annotation #fewestMM_BUL
                # splitted[58] = last_annotation #bestCRISTA

                tab = '\t'
                f_out.write(f"{tab.join(splitted)}\n")


def annotate_results(inAnnotationFile: str, files: List[Tuple[str, str]]):
    genome_to_process = extractGenomesFromResults([f[0] for f in files])
    annotationDict = loadAnnotationsFiltered(inAnnotationFile, genome_to_process)
    for (file_result, file_annotated) in files:
        writeAnnotated(file_result, file_annotated, annotationDict)


if __name__ == '__main__':
    
    # logic to activate new version: parameters starts with a dash or the number of arguments is wrong
    # in the old logic
    new_version = len(sys.argv) != 4 or any([a.startswith('-') for a in sys.argv[1:]])

    inAnnotationFile = ''
    file_final_results = []
    file_annotated = []
    if new_version:
        import argparse
        parser=argparse.ArgumentParser(prog='Annotate Results',
                    description='This file combine the output of the genomes with the enrichment',
                    # epilog='Text at the bottom of help'
                    )
        parser.add_argument('annotation_file', help='The annotation file to load the annotations.')
        parser.add_argument('-p','--process',action='append',nargs=2, required=True,
            metavar=('input','output'),help='The pair that specifies input file and output (enriched)')
        args = parser.parse_args()

        inAnnotationFile = args.annotation_file
        file_final_results = [f[0] for f in args.process]
        file_annotated = [f[1] for f in args.process]
    else:
        file_final_results = [sys.argv[1]]
        inAnnotationFile = sys.argv[2]
        file_annotated = [sys.argv[3]]

    print("Starting annotation")
    start_time = time.time()
    annotate_results(inAnnotationFile, [(result, output) for result, output in zip(file_final_results, file_annotated)])
    print(f"Annotation done in {time.time()-start_time}")
