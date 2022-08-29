import sys
import pandas as pd

print('start processing')
original_df = pd.read_csv(sys.argv[1], sep="\t", index_col=False,
                          na_values=['n'])
data_frames_list = list()

for guide in original_df['Spacer+PAM'].unique():

    guide_df = original_df.loc[(original_df['Spacer+PAM'] == guide)]

    alt_target_count_CFD = guide_df.apply(lambda x: True
                                          if x['REF/ALT_origin_(highest_CFD)'] == "alt" else False, axis=1)
    alt_target_count_MMvBUL = guide_df.apply(lambda x: True
                                             if x['REF/ALT_origin_(fewest_mm+b)'] == "alt" else False, axis=1)

    # Count number of True in the series
    guide_df['alt_target_cfd'] = len(
        alt_target_count_CFD[alt_target_count_CFD == True].index)

    guide_df['alt_target_mmvbul'] = len(
        alt_target_count_MMvBUL[alt_target_count_MMvBUL == True].index)

    guide_df['total_target'] = len(guide_df.index)

    on_target_df = guide_df.loc[(
        guide_df['Mismatches+bulges_(fewest_mm+b)'] <= 1)]

    data_frames_list.append(on_target_df[['Spacer+PAM', 'Chromosome',
                                          'Start_coordinate_(fewest_mm+b)',
                                          'Aligned_protospacer+PAM_REF_(fewest_mm+b)',
                                          'Aligned_protospacer+PAM_ALT_(fewest_mm+b)',
                                          'Mismatches+bulges_(fewest_mm+b)',
                                          'Annotation_closest_gene_name', 'total_target',
                                          'alt_target_cfd', 'alt_target_mmvbul']])

final_df = pd.concat(data_frames_list)

final_df.to_csv(sys.argv[2]+str(sys.argv[1]).split('/')[-1]+'_extracted_target.tsv',
                sep='\t', na_rep='NA', index=False)
