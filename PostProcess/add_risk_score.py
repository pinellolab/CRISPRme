#!/usr/bin/env python

import sys


file_in = sys.argv[1]
file_out = sys.argv[2]
alt = sys.argv[3]
if alt == 'True':
    alt = True
else:
    alt = False

with open(file_in, 'r') as fin:
    with open(file_out, 'w') as fout:
        header = fin.readline().strip().split('\t')
        # header.insert(22, 'Highest_CFD_Risk_Score')
        # header.insert(23, 'Highest_CFD_Absolute_Risk_Score')
        # header.append('MMBLG_CFD_Risk_Score')
        # header.append('MMBLG_CFD_Absolute_Risk_Score')
        header.append('Highest_CFD_Risk_Score')
        header.append('Highest_CFD_Absolute_Risk_Score')
        if alt:
            header.append('CLUSTER_ID')
        fout.write('\t'.join(header)+'\n')
        for line in fin:
            splitted = line.strip().split('\t')
            cfd_diff = float(splitted[20]) - float(splitted[21])
            abs_diff = abs(cfd_diff)
            # mmblg_cfd_diff = float(splitted[42]) - float(splitted[43])
            # mmblg_abs_diff = abs(mmblg_cfd_diff)
            fout.write('\t'.join(splitted)+'\t' +
                       str(cfd_diff)+'\t'+str(abs_diff)+'\n')
            # if alt:
            #     fout.write('\t'.join(splitted[:22])+'\t'+"{:.3f}".format(cfd_diff)+'\t'+"{:.3f}".format(abs_diff)+"\t"+"\t".join(
            #         splitted[22:-1])+'\t'+"{:.3f}".format(mmblg_cfd_diff)+'\t'+"{:.3f}".format(mmblg_abs_diff)+"\t"+splitted[-1]+'\n')
            # else:
            #     fout.write('\t'.join(splitted[:22])+'\t'+"{:.3f}".format(cfd_diff)+'\t'+"{:.3f}".format(abs_diff)+"\t"+"\t".join(
            #         splitted[22:])+'\t'+"{:.3f}".format(mmblg_cfd_diff)+'\t'+"{:.3f}".format(mmblg_abs_diff)+'\n')
