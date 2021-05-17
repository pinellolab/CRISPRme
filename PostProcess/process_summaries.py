#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys
import os


path_best = sys.argv[1]
#path_alt = sys.argv[2]
path_guides = sys.argv[2]
path_samplesID = sys.argv[3]
mms = int(sys.argv[4])
bulge = int(sys.argv[5])
path_output = sys.argv[6]
genome_type = sys.argv[7]


def init_summary_by_sample(path_samplesID):
    dict_samples = {}
    dict_pop = {}
    dict_superpop = {}
    with open(f"{path_samplesID}", "r") as sid:
        sid.readline()
        for line in sid:
            splitted = line.strip().split("\t")
            dict_samples[splitted[0]] = [splitted[3],
                                         splitted[1], splitted[2], 0, 0, 0, 0]
            if splitted[1] not in dict_pop.keys():
                dict_pop[splitted[1]] = 0
            if splitted[2] not in dict_superpop.keys():
                dict_superpop[splitted[2]] = 0

    return dict_samples, dict_pop, dict_superpop

#os.system(f"cp {path_best} {path_best}.cumulative")
#os.system(f"tail -n +2 {path_alt} >> {path_best}.cumulative")


guides = []
general_table = {}
with open(f"{path_guides}", "r") as guide_file:
    for line in guide_file:
        if line.strip():
            stripped = line.strip()
            guides.append(stripped)

            general_table[stripped] = dict()
            general_table[stripped]['ref'] = np.zeros((bulge+1, mms+1), dtype=int)
            general_table[stripped]['var'] = np.zeros((bulge+1, mms+1), dtype=int)

add_to_general_table = {}
count_superpop = {}
# with open(path_output + '.general_target_count.txt', 'w+') as general_count:
count_for = '(' + ' - '.join([str(tot)
                                for tot in range(1, mms + bulge + 1)]) + ' Mismatches + Bulges)'
# general_count.write('#Guide\tOn-Targets (Reference - Enriched)\tOff-Targets Reference ' +
#                     count_for + '\tOff-Targets Enriched ' + count_for + '\n')

for guide in guides:
    os.system(f"LC_ALL=C fgrep {guide} {path_best} > {path_best}.{guide}")
    #os.system(f"LC_ALL=C fgrep {guide} {path_alt} >> {path_best}.{guide}")

    sum_cfds = 0
    on_targets = []

    dict_summary_by_guide = {}
    add_to_general_table[guide] = {}
    add_to_general_table[guide]['distributions'] = [
        [0] * (bulge + 1) for x in range(10)]
    if genome_type == "var":
        dict_samples, dict_pop, dict_superpop = init_summary_by_sample(
            path_samplesID)
        count_superpop[guide] = {}
        for superpop in dict_superpop.keys():
            count_superpop[guide][superpop] = {}
            count_superpop[guide][superpop]['distributions'] = [
                [0] * (bulge + 1) for x in range(10)]

    with open(f"{path_best}.{guide}", "r") as f:
        for line in f:
            splitted = line.strip().split("\t")
            var = True
            pam_creation = splitted[11] != "n"

            if splitted[13] == "n":
                var = False
                general_table[guide]['ref'][int(
                    splitted[9]), int(splitted[8])] += 1
            else:
                general_table[guide]['var'][int(
                    splitted[9]), int(splitted[8])] += 1

            bulge_type = f"{splitted[0]}\t{splitted[8]}\t{splitted[9]}"
            if bulge_type in dict_summary_by_guide.keys():
                if var:
                    dict_summary_by_guide[bulge_type][1] += 1
                else:
                    dict_summary_by_guide[bulge_type][0] += 1
                dict_summary_by_guide[bulge_type][2] += 1
            else:
                if var:
                    dict_summary_by_guide[bulge_type] = [0, 1, 1, 0]
                else:
                    dict_summary_by_guide[bulge_type] = [1, 0, 1, 0]
            if pam_creation:
                dict_summary_by_guide[bulge_type][3] += 1

            if var:
                samples = set(splitted[13].split(","))
                seen_superpop = set()
                seen_pop = set()
                for sample in samples:
                    if sample != "NO_SAMPLES" and sample != '':
                        dict_samples[sample][3] += 1
                        seen_pop.add(dict_samples[sample][1])
                        seen_superpop.add(dict_samples[sample][2])
                        if pam_creation:
                            dict_samples[sample][6] += 1
                for superpop in seen_superpop:
                    dict_superpop[superpop] += 1
                    count_superpop[guide][superpop]['distributions'][int(
                        splitted[8]) + int(splitted[9])][int(splitted[9])] += 1
                for pop in seen_pop:
                    dict_pop[pop] += 1
            else:
                add_to_general_table[guide]['distributions'][int(
                    splitted[8]) + int(splitted[9])][int(splitted[9])] += 1

            sum_cfds += float(splitted[20])
            if int(splitted[10]) == 0:
                on_targets.append(line)

    if genome_type == "var":
        for sample in dict_samples.keys():
            dict_samples[sample][4] = dict_pop[dict_samples[sample][1]]
            dict_samples[sample][5] = dict_superpop[dict_samples[sample][2]]

    with open(f"{path_output}.summary_by_guide.{guide}.txt", "w") as summary_by_guide:
        summary_by_guide.write(
            "Bulge Type\tMismatches\tBulge Size\tReference\tEnriched\tCombined\tPAM Creation\n")
        for key in dict_summary_by_guide.keys():
            entry = dict_summary_by_guide[key]
            summary_by_guide.write(
                f"{key}\t{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")

    if genome_type == "var":
        with open(f"{path_output}.summary_by_samples.{guide}.txt", "w") as summary_by_samples:
            summary_by_samples.write(guide+"\n")
            for key in dict_samples.keys():
                entry = dict_samples[key]
                summary_by_samples.write(
                    f"{key}\t{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\t{entry[6]}\n")

    with open(f"{path_output}.acfd.txt", "a") as acfd:
        if sum_cfds == 0:
            acfd.write(guide+"\t0\tNA\tNA\n")
        else:
            acfd.write(guide+"\t"+str(100/(100+sum_cfds))+"\tNA\tNA\n")

    df_general_count = pd.DataFrame(general_table[guide]['ref'])
    df_general_count = df_general_count.append(
        pd.DataFrame(general_table[guide]['var']))
    df_general_count.to_csv(
        path_output+".general_target_count."+guide+".txt", sep='\t', index=False)
    # general_count.write(guide + '\t' + str(general_table[guide]['ref'][0] + general_table[guide]['var'][0]) + ' (' + str(general_table[guide]['ref'][0]) + ' - ' + str(general_table[guide]['var'][0]) + ')\t' +
    #                    str(sum(general_table[guide]['ref'][1:])) + ' (' + ' - '.join([ str(x) for x in general_table[guide]['ref'][1:]]) + ')\t' +
    #					str(sum(general_table[guide]['var'][1:])) + ' (' + ' - '.join([ str(x) for x in general_table[guide]['var'][1:]]) + ')\n'
    #					)

    os.system(
        f"sort -k2,2n -k3,3n {path_output}.summary_by_guide.{guide}.txt -o {path_best}.summary_by_guide.{guide}.sorted")
    os.system(
        f"mv {path_best}.summary_by_guide.{guide}.sorted {path_output}.summary_by_guide.{guide}.txt")
    os.system(f"rm {path_best}.{guide}")


with open(path_output + '.PopulationDistribution.txt', 'w+') as pop_distribution:
    for g in guides:
        pop_distribution.write('-Summary_' + g + '\n')
        pop_distribution.write('REF' + '\t' + '\t'.join([','.join(
            str(v) for v in t) for t in add_to_general_table[g]['distributions']]) + '\n')
        if genome_type == "var":
            for superpop in dict_superpop.keys():
                try:
                    pop_distribution.write(superpop + '\t' + '\t'.join([','.join(
                        str(v) for v in t) for t in count_superpop[g][superpop]['distributions']]) + '\n')
                except:
                    pop_distribution.write(superpop + '\t' + '\t'.join([','.join(
                        str(v) for v in t) for t in [[0] * (bulge + 1) for x in range(10)]]) + '\n')
