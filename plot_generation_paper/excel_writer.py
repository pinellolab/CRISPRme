from operator import truediv
import sys
import pandas as pd


def count_ref_mm(aligned_ref_seq):
    count_mm = 0
    for nt in str(aligned_ref_seq):
        if nt.islower():
            count_mm += 1
    return count_mm


print('start processing')

if len(sys.argv[:]) < 4:
    print('some input is missing, please provide input')
    print('integrated.tsv out_dir sort_criteria(CFD/CRISTA/fewest)')
    exit(1)

# df with targets
original_df = pd.read_csv(sys.argv[1], sep="\t", index_col=False,
                          na_values=['n'],nrows=200000)
# excel writer
writer = pd.ExcelWriter(sys.argv[2]+str(sys.argv[2].strip().split('/')[1])+'_guide_sheets.xlsx')
# user sort criteria
sort_criteria = sys.argv[3]
if sort_criteria == 'CFD':
    sort_criteria = 'highest_CFD'
elif sort_criteria == 'CRISTA':
    sort_criteria = 'highest_CRISTA'
elif sort_criteria == 'fewest':
    sort_criteria = 'fewest_mm+b'

for guide in sorted(original_df['Spacer+PAM'].unique().tolist()):
    # filter df for guide
    guide_df = original_df.loc[(original_df['Spacer+PAM'] == guide)]

    drop_criteria = list()
    # sort df using user criteria
    if 'CFD' in sort_criteria:
        guide_df.sort_values('CFD_score_(highest_CFD)',
                             ascending=False, inplace=True)
        drop_criteria.append('fewest_mm+b')
        drop_criteria.append('highest_CRISTA')
    elif 'CRISTA' in sort_criteria:
        guide_df.sort_values('CRISTA_score_(highest_CRISTA)',
                             ascending=False, inplace=True)
        drop_criteria.append('fewest_mm+b')
        drop_criteria.append('highest_CFD')
    elif 'fewest' in sort_criteria:
        guide_df.sort_values('Mismatches+bulges_(fewest_mm+b)',
                             ascending=True, inplace=True)
        drop_criteria.append('highest_CFD')
        drop_criteria.append('highest_CRISTA')

    columns_to_drop = list()
    for column in list(guide_df.columns):
        if any(criteria in column for criteria in drop_criteria):
            columns_to_drop.append(column)
    guide_df.drop(columns_to_drop, axis=1, inplace=True)

    # extract top 1000 rows for each guide
    guide_df = guide_df.head(1000)
    # guide_df.insert(guide_df.columns.get_loc('Mismatches+bulges_(fewest_mm+b)'),
    #                 'Mismatches+bulges_(fewest_mm+b)_REF', 0)
    # guide_df.rename(columns={
    #                 "Mismatches+bulges_(fewest_mm+b)": "Mismatches+bulges_(fewest_mm+b)_ALT"}, inplace=True)
    # guide_df['Mismatches+bulges_(fewest_mm+b)_REF'] = guide_df['Aligned_protospacer+PAM_REF_(fewest_mm+b)'].apply(
    #     lambda x: count_ref_mm(x))
    # guide_df['Mismatches+bulges_(fewest_mm+b)_REF'] += guide_df['Bulges_(fewest_mm+b)']

    # generate excel sheets
    guide_df.to_excel(writer, sheet_name=str(guide), na_rep='NA', index=False)

# save the writer to excel file
writer.save()
writer.close()
