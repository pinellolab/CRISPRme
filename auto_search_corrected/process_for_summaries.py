#!/usr/bin/env python

import sys
import os
import pandas as pd


def init_dict_samples(path_samples):
    dict_samples = {}
    dict_population = {}
    dict_superpopulation = {}
    with open(path_samples, 'r') as s:
        header = s.readline()
        for line in s:
            splitted = line.strip().split('\t')
            dict_samples[splitted[0]] = [splitted[3],
                                         splitted[1], splitted[2], 0, 0, 0, 0]
            if splitted[1] not in dict_population.keys():
                dict_population[splitted[1]] = 0
            if splitted[2] not in dict_superpopulation.keys():
                dict_superpopulation[splitted[2]] = 0

    return dict_samples, dict_population, dict_superpopulation


path_final_results_best = sys.argv[1]
path_final_results_alt = sys.argv[2]
path_guides = sys.argv[3]
path_samples = sys.argv[4]

os.system(
    f"cat {path_final_results_best} > {path_final_results_best}.cumulative")
os.system(
    f"tail -n +2 {path_final_results_alt} >> {path_final_results_best}.cumulative")

guides = []
with open(path_guides, 'r') as g:
    for line in g:
        guides.append(line)


for guide in guides:
    os.system(
        f"grep {guide} {path_final_results_best}.cumulative > {path_final_results_best}.cumulative.{guide}")
    os.system(
        f"sort -T {os.path.dirname(sys.argv[1])} -k5,5 -k6,6n {path_final_results_best}.cumulative.{guide} -o {path_final_results_best}.cumulative.{guide}.sorted")
    dict_guide_summary = {}
    dict_samples_summary, dict_population, dict_superpopulation = init_dict_samples(
        path_samples)
    with open(f"{path_final_results_best}.cumulative.sorted.{guide}", 'r') as f:
        header = f.readline()
        for line in f:
            splitted = line.strip().split('\t')
            var = False
            pam_creation = False
            if splitted[13] != 'n':  # check if we have at least one sample
                var = True
            if splitted[11] != 'n':  # check pam_creation
                pam_creation = True
            # BulgeType\tMismatches\tBulges
            class_target = splitted[0] + "\t" + \
                splitted[8] + '\t' + splitted[9]
            if class_target in dict_guide_summary.keys():
                if var:
                    dict_guide_summary[class_target][1] += 1
                else:
                    dict_guide_summary[class_target][0] += 1
                dict_guide_summary[class_target][2] += 1
                if pam_creation:
                    dict_guide_summary[class_target][3] += 1
            else:
                if var:
                    if pam_creation:
                        dict_guide_summary[class_target] = [0, 1, 1, 1]
                    else:
                        dict_guide_summary[class_target] = [0, 1, 1, 0]
                else:
                    if pam_creation:
                        dict_guide_summary[class_target] = [1, 0, 1, 1]
                    else:
                        dict_guide_summary[class_target] = [1, 0, 1, 0]

            if var:
                samples = splitted[13].split(',')
                for sample in samples:
                    if sample != 'NO_SAMPLES':
                        dict_samples_summary[sample][3] += 1
                        if pam_creation:
                            dict_samples_summary[sample][6] += 1
                        dict_population[dict_samples_summary[sample][1]] += 1
                        dict_superpopulation[dict_samples_summary[sample][2]] += 1

        for key in dict_samples_summary.keys():
            dict_samples_summary[key][4] = dict_population[dict_samples_summary[key][1]]
            dict_samples_summary[key][5] = dict_superpopulation[dict_samples_summary[key][2]]

    print(dict_guide_summary)
    print(dict_samples_summary)
