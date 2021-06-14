#!/usr/bin/env python
"""
Created on Sat May 29 18:02:45 2021

@author: franc
"""

'''
Script used to convert from old bestMerge format to new alt_results format
'''
import pandas as pd
import numpy as np
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

in_path = sys.argv[1]
out_path = sys.argv[2]

chunksize_ = 5000000000000

# if not isFileLocked(in_path):
chunks = pd.read_csv(in_path, sep='\t', chunksize=chunksize_, na_filter=False)

new_order = ['Real_Guide', 'Chromosome', 'Position', 'Direction', 'Cluster_Position', 'crRNA', 'Reference', 'DNA', 'Mismatches', 'Bulge_Size', 'Total',
             '#Bulge_type', 'PAM_gen', 'CFD', 'CFD_ref', 'Highest_CFD_Risk_Score', 'Highest_CFD_Absolute_Risk_Score', 'SNP', 'AF', 'rsID',
             'Samples', 'Var_uniq', '#Seq_in_cluster', 'MMBLG_Real_Guide', 'MMBLG_Chromosome', 'MMBLG_Position', 'MMBLG_Cluster_Position', 'MMBLG_Direction', 'MMBLG_crRNA',
             'MMBLG_Reference', 'MMBLG_DNA', 'MMBLG_Mismatches', 'MMBLG_Bulge_Size', 'MMBLG_Total', 'MMBLG_#Bulge_type', 'MMBLG_PAM_gen', 'MMBLG_CFD', 'MMBLG_CFD_ref',
             'MMBLG_CFD_Risk_Score', 'MMBLG_CFD_Absolute_Risk_Score', 'MMBLG_Var_uniq', 'MMBLG_SNP', 'MMBLG_AF', 'MMBLG_rsID', 'MMBLG_Samples', 'MMBLG_#Seq_in_cluster',
             'MMBLG_Annotation_Type', 'Annotation_Type']
# print(len(new_order))

to_remove = ['Real_Guide', 'Cluster_Position', 'Highest_CFD_Absolute_Risk_Score', 'MMBLG_Real_Guide', 'MMBLG_Chromosome', 'MMBLG_Cluster_Position',
             'MMBLG_CFD_Absolute_Risk_Score', 'MMBLG_Var_uniq', 'MMBLG_#Seq_in_cluster', 'MMBLG_Annotation_Type'
             ]
# print(len(to_remove))

new_names = ['Chromosome', 'Start_coordinate_(highest_CFD)', 'Strand_(highest_CFD)', 'Aligned_spacer+PAM_(highest_CFD)', 'Aligned_protospacer+PAM_REF_(highest_CFD)', 'Aligned_protospacer+PAM_ALT_(highest_CFD)',
             'Mismatches_(highest_CFD)', 'Bulges_(highest_CFD)', 'Mismatches+bulges_(highest_CFD)', 'Bulge_type_(highest_CFD)', 'PAM_creation_(highest_CFD)', 'CFD_score_(highest_CFD)', 'CFD_score_REF_(highest_CFD)',
             'CFD_risk_score_(highest_CFD)', 'Variant_info_genome_(highest_CFD)', 'Variant_MAF_(highest_CFD)', 'Variant_rsID_(highest_CFD)', 'Variant_samples_(highest_CFD)', 'Not_found_in_REF', 'Other_motifs', 
             'Start_coordinate_(fewest_mm+b)', 'Strand_(fewest_mm+b)', 'Aligned_spacer+PAM_(fewest_mm+b)', 'Aligned_protospacer+PAM_REF_(fewest_mm+b)', 'Aligned_protospacer+PAM_ALT_(fewest_mm+b)', 
             'Mismatches_(fewest_mm+b)', 'Bulges_(fewest_mm+b)', 'Mismatches+bulges_(fewest_mm+b)', 'Bulge_type_(fewest_mm+b)', 'PAM_creation_(fewest_mm+b)', 'CFD_score_(fewest_mm+b)', 
             'CFD_score_REF_(fewest_mm+b)', 'CFD_risk_score_(fewest_mm+b)', 'Variant_info_genome_(fewest_mm+b)', 'Variant_MAF_(fewest_mm+b)', 'Variant_rsID_(fewest_mm+b)', 
             'Variant_samples_(fewest_mm+b)', 'Annotation']


header = True
for chunk in chunks:

    chunk = chunk[new_order]
    chunk = chunk.drop(to_remove, axis=1)
    chunk.columns = new_names
    #chunk[chunk.columns.difference(['Not_found_in_REF'])] = chunk[chunk.columns.difference(['Not_found_in_REF'])].replace('n', 'NA')
    chunk = chunk.replace('n', 'NA')
    #chunk = chunk.replace(regex=['\*.,\*', '\*,.\*'], value='NA')
    chunk['Variant_rsID_(highest_CFD)'] = chunk['Variant_rsID_(highest_CFD)'].str.replace(
        '.', 'NA')
    mask = chunk['Aligned_protospacer+PAM_REF_(highest_CFD)'] == 'NA'
    chunk['Aligned_protospacer+PAM_REF_corrected_(highest_CFD)'] = np.where(
        mask, chunk['Aligned_protospacer+PAM_ALT_(highest_CFD)'], chunk['Aligned_protospacer+PAM_REF_(highest_CFD)'])
    chunk['Aligned_protospacer+PAM_ALT_(highest_CFD)'] = np.where(
        mask, chunk['Aligned_protospacer+PAM_REF_(highest_CFD)'], chunk['Aligned_protospacer+PAM_ALT_(highest_CFD)'])

    chunk['Aligned_protospacer+PAM_REF_(highest_CFD)'] = chunk['Aligned_protospacer+PAM_REF_corrected_(highest_CFD)']
    chunk.drop('Aligned_protospacer+PAM_REF_corrected_(highest_CFD)',
               axis=1, inplace=True)

    mask = chunk['Aligned_protospacer+PAM_REF_(fewest_mm+b)'] == 'NA'
    chunk['Aligned_protospacer+PAM_REF_corrected_(fewest_mm+b)'] = np.where(
        mask, chunk['Aligned_protospacer+PAM_ALT_(fewest_mm+b)'], chunk['Aligned_protospacer+PAM_REF_(fewest_mm+b)'])
    chunk['Aligned_protospacer+PAM_ALT_(fewest_mm+b)'] = np.where(
        mask, chunk['Aligned_protospacer+PAM_REF_(fewest_mm+b)'], chunk['Aligned_protospacer+PAM_ALT_(fewest_mm+b)'])

    chunk['Aligned_protospacer+PAM_REF_(fewest_mm+b)'] = chunk['Aligned_protospacer+PAM_REF_corrected_(fewest_mm+b)']
    chunk.drop('Aligned_protospacer+PAM_REF_corrected_(fewest_mm+b)',
               axis=1, inplace=True)

    chunk.to_csv(out_path, header=header, mode='w',
                 sep='\t', index=False, na_rep='NA')

    header = False
