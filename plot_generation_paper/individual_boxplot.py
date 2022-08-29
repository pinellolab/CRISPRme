import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import warnings
from statistics import mean
import scipy

# SUPPRESS ALL WARNINGS
warnings.filterwarnings("ignore")
# do not use X11
matplotlib.use('Agg')
# set matplotlib for pdf editing
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_context("paper")

df_single_search = pd.read_csv(sys.argv[1], sep="\t", index_col=False, na_values=[
                               'n'], usecols=['Variant_samples_(highest_CFD)'])
df_double_search = pd.read_csv(sys.argv[2], sep="\t", index_col=False, na_values=[
                               'n'], usecols=['Variant_samples_(highest_CFD)'])
sample_file = open(sys.argv[3], 'r')
out_folder = sys.argv[4]
analyzed_set = sys.argv[5]+' Variants'

print('doing search', analyzed_set)


def count_ratio(boxplot_values, sample_dict: dict):
    for sample in sample_dict:
        # skip sample with 0 personal
        if sample_dict[sample][1] == 0:
            continue

        # calculate ratio with shared targets
        ratio = sample_dict[sample][0]/sample_dict[sample][1]
        boxplot_values[0].append(ratio)
        boxplot_values[1].append(sample_dict[sample][0])
        boxplot_values[2].append(sample_dict[sample][1])
        boxplot_values[3].append(str(sample_dict[sample][2]))


def count_personal_and_private(sample_string: str, sample_dict: dict):
    sample_list = sample_string.strip().split(',')
    for sample in sample_list:
        try:
            sample_dict[sample][1] += 1
            if len(sample_list) == 1:
                sample_dict[sample][0] += 1
        except:
            continue


sample_dict_single = dict()
sample_dict_double = dict()
for line in sample_file:
    if '#' in line:
        continue
    splitted = line.strip().split('\t')
    #sample_dict[individual] =[private,personal,superpop]
    sample_dict_single[splitted[0]] = [0, 0, splitted[2]]
    sample_dict_double[splitted[0]] = [0, 0, splitted[2]]

# analyze search
df_single_search['Variant_samples_(highest_CFD)'].apply(
    lambda x: count_personal_and_private(str(x), sample_dict_single))
df_double_search['Variant_samples_(highest_CFD)'].apply(
    lambda x: count_personal_and_private(str(x), sample_dict_double))

# list containing lists ratio for private_single_search/personal_single_search
boxplot_values_single_search = [[], [], [], []]
# private_double_search/personal_double_search
boxplot_values_double_search = [[], [], [], []]
count_ratio(boxplot_values_single_search, sample_dict_single)
count_ratio(boxplot_values_double_search, sample_dict_double)

print('mean value single search', mean(boxplot_values_single_search[0]))
print('mean value double search', mean(boxplot_values_double_search[0]))
print('t-test',scipy.stats.ttest_ind(boxplot_values_single_search[0],boxplot_values_double_search[0]))
df_complete = pd.DataFrame(
    {str(analyzed_set): boxplot_values_single_search[0], 'Private': boxplot_values_single_search[1], 'Personal': boxplot_values_single_search[2], '1000G+HGDP Variants': boxplot_values_double_search[0], 'Superpopulation': boxplot_values_single_search[3]})
# print(df_complete[str(analyzed_set)])
# print(df_complete['1000G+HGDP'])

# SCATTERPLOT
# plt.figure()
color_dict = {'AFR': 'tab:orange', 'AMR': 'tab:brown', 'CSA': 'tab:blue',
              'EAS': 'tab:pink', 'EUR': 'tab:red', 'MEA': 'tab:purple', 'OCE': 'tab:green', 'SAS': 'tab:cyan'}
sns.scatterplot(data=df_complete, x='Private', y='Personal',
                hue='Superpopulation', rasterized=True, palette=color_dict, alpha=0.5, linewidth=0)

# ax = sns.violinplot(x='Private', y='Personal', data=df_complete)
# for violin in ax.collections[::2]:
#     violin.set_alpha(0.2)
# ax = sns.stripplot(x='Private', y='Personal', data=df_complete)

plt.title(str(analyzed_set))
plt.tight_layout()
plt.savefig(
    out_folder+f"{analyzed_set}_scatterplot.pdf")
plt.clf()
plt.close('all')

# DISTPLOT
plt.figure()

# sns.displot(df_complete[[str(analyzed_set), '1000G+HGDP']])
# print(df_complete[[str(analyzed_set), '1000G+HGDP']])
ax = sns.violinplot(data=df_complete[[str(analyzed_set), '1000G+HGDP Variants']])
# for violin in ax.collections[::2]:
#     violin.set_alpha(0.2)
# ax = sns.stripplot(
#     data=df_complete[[str(analyzed_set), '1000G+HGDP']])
plt.title(str(analyzed_set).replace('Variants','Individuals'))
# plt.xlabel('Variant dataset')
plt.ylabel('Ratio of private/personal targets')
plt.tight_layout()
plt.savefig(
    out_folder+f"{analyzed_set}_distplot.pdf")
plt.clf()
plt.close('all')
