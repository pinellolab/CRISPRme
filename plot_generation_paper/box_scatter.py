import sys
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import pandas as pd
import numpy as np
matplotlib.use('Agg')


print('creating dataframe')
out_folder=sys.argv[2] #outfolder for images
df = pd.read_csv(sys.argv[1], sep="\t", index_col=False, na_values=['n'])

# plt.figure(figsize=(15, 5))

ax = sns.stripplot(x="origin", y="value", data=df)
ax=sns.boxplot(x="origin", y="value", data=df)
for patch in ax.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b,0.2))

plt.ylabel('ALT_Targets/Top 100 targets')
plt.tight_layout()
plt.savefig(out_folder+"box_plot.png", bbox_inches="tight")
plt.clf()
plt.close('all')
