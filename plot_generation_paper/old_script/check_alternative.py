import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t", index_col=False,
                 na_values=['n'])


df.sort_values('CFD_score_(highest_CFD)', ascending=False, inplace=True)
dff_alt_CFD = df.head(100)
dff_alt_CFD.to_csv('top100_CFD',
                   sep='\t', na_rep='NA', index=False)
dff_alt_CFD = dff_alt_CFD.loc[(
    dff_alt_CFD['REF/ALT_origin_(highest_CFD)'] == 'alt')]
print('top100 row by CFD and filtered with alt', len(dff_alt_CFD))

df.sort_values("Mismatches+bulges_(fewest_mm+b)", ascending=True, inplace=True)
dff_alt_total = df.head(100)
dff_alt_total.to_csv('top100_total',
                     sep='\t', na_rep='NA', index=False)
dff_alt_total = dff_alt_total.loc[(
    dff_alt_total['REF/ALT_origin_(fewest_mm+b)'] == 'alt')]
print('top100 row by total and filtered with alt', len(dff_alt_total))
