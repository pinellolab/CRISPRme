from textwrap import fill
from scipy import stats
import sys
import pandas as pd
import matplotlib
import seaborn as sns
from matplotlib import pyplot as plt
# from matplotlib_venn import venn2
import warnings
import numpy as np

# SUPPRESS ALL WARNINGS
warnings.filterwarnings("ignore")
# do not use X11
matplotlib.use('Agg')
# set matplotlib for pdf editing
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.style.use('seaborn-poster')
sns.set_context("paper")


# INPUT
# ARGV1 INTEGRATED FILE
# ARGV2 OUTPUT FOLDER


# color palette for hue
palette = {
    '0': 'tab:blue',
    '1': 'tab:green',
    '2': 'tab:orange'
}


def bulge_color(row: pd.Series):
    lower_bulge = min(int(row['Bulges_(highest_CFD)']),
                      int(row['Bulges_(highest_CRISTA)']))
    return str(lower_bulge)


def plot_correlation(original_df: pd.DataFrame, max_bulges: int):

    print('start plotting data')

    # df list with all guides dfs
    df_guide_list = list()
    count_list = list()

    for guide in original_df["Spacer+PAM"].unique():

        # filter the df to obtain single guide targets
        df_guide = original_df.loc[(original_df['Spacer+PAM'] == guide)]
        df_guide.sort_values(['Chromosome'], inplace=True)
        # rank scores in df
        df_guide['CFD_Rank'] = df_guide['CFD_score_(highest_CFD)'].rank(
            method='first', ascending=False)
        df_guide['CRISTA_Rank'] = df_guide['CRISTA_score_(highest_CRISTA)'].rank(
            method='first', ascending=False)
        df_guide['CFD_Rank'] = df_guide['CFD_Rank'].astype(int)
        df_guide['CRISTA_Rank'] = df_guide['CRISTA_Rank'].astype(int)

        df_CFD_top10000 = df_guide.loc[(
            df_guide['CFD_Rank'] <= 10000)]
        df_CFD_top10000.to_csv(
            sys.argv[2]+f'top10000_CFD_{guide}_{max_bulges}_bulges.tsv', sep='\t', na_rep='NA', index=False)

        df_CRISTA_top10000 = df_guide.loc[(
            df_guide['CRISTA_Rank'] <= 10000)]
        df_CRISTA_top10000.to_csv(
            sys.argv[2]+f'top10000_CRISTA_{guide}_{max_bulges}_bulges.tsv', sep='\t', na_rep='NA', index=False)

        plt.figure(figsize=(5, 5))
        sns.scatterplot(data=df_CFD_top10000, x="CFD_score_(highest_CFD)", y='CRISTA_score_(highest_CRISTA)',
                        hue='Bulge_count', rasterized=True, palette=palette, alpha=0.5, legend=False, linewidth=0)
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.tight_layout()
        plt.savefig(
            sys.argv[2]+f'corr_plot_top10000_CFD_{guide}_{max_bulges}_bulges.pdf', dpi=300)
        plt.clf()
        plt.close('all')

        plt.figure(figsize=(5, 5))
        sns.scatterplot(data=df_CRISTA_top10000, x="CFD_score_(highest_CFD)", y='CRISTA_score_(highest_CRISTA)',
                        hue='Bulge_count', rasterized=True, palette=palette, alpha=0.5, legend=False, linewidth=0)
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.tight_layout()
        plt.savefig(
            sys.argv[2]+f'corr_plot_top10000_CRISTA_{guide}_{max_bulges}_bulges.pdf', dpi=300)
        plt.clf()
        plt.close('all')

        # select only the top10000 for CFD and CRISTA ranks
        df_guide = df_guide.loc[(
            df_guide['CFD_Rank'] <= 10000) | (df_guide['CRISTA_Rank'] <= 10000)]
        # pair all ranks >10000 to 10000
        df_guide.loc[df_guide['CFD_Rank']
                     > 10000, 'CFD_Rank'] = 10000
        df_guide.loc[df_guide['CRISTA_Rank']
                     > 10000, 'CRISTA_Rank'] = 10000

        # save the guide df to file
        df_guide.to_csv(
            sys.argv[2]+f'original_{guide}_{max_bulges}_bulges.tsv', sep='\t', na_rep='NA', index=False)
        # append dfguide to df_guide list to generate the final df with all guides
        df_guide_list.append(df_guide)

        plt.figure(figsize=(5, 5))
        sns.scatterplot(data=df_guide, x="CFD_score_(highest_CFD)", y='CRISTA_score_(highest_CRISTA)',
                        hue='Bulge_count', rasterized=True, palette=palette, alpha=0.5, legend=False, linewidth=0)

        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.tight_layout()
        plt.savefig(
            sys.argv[2]+f'correlation_CFDvCRISTA_top10000_union_for_{guide}_with_{max_bulges}_bulges.pdf', dpi=300)
        plt.clf()
        plt.close('all')

        print('plot for guide top10000', guide)
        guide_list = list()
        guide_list.append(guide)
        # count total targets in union
        guide_list.append(len(df_guide.index))
        # CFD<100 & CRISTA<100
        guide_list.append(len(df_guide[(df_guide.CFD_Rank <= 100) & (
            df_guide.CRISTA_Rank <= 100)].index))
        # CFD<100 & CRISTA>100
        guide_list.append(len(df_guide[(df_guide.CFD_Rank <= 100) & (
            df_guide.CRISTA_Rank > 100)].index))
        # CFD>100 & CRISTA<100
        guide_list.append(len(df_guide[(df_guide.CFD_Rank > 100) & (
            df_guide.CRISTA_Rank <= 100)].index))
        # CFD>100 & CRISTA>100
        guide_list.append(len(df_guide[(df_guide.CFD_Rank > 100) & (
            df_guide.CRISTA_Rank > 100)].index))
        # append all counts to single list
        count_list.append(guide_list)

    # whole figure
    final_df = pd.concat(df_guide_list)
    final_df.to_csv(
        sys.argv[2] + f'final_df_ranks_with_{max_bulges}_bulges.tsv', sep='\t', na_rep='NA', index=False)
    count_df = pd.DataFrame(count_list, columns=[
                            'sgRNA', 'total_targets_union', 'CFD<=100&CRISTA<=100', 'CFD<=100&CRISTA>100', 'CFD>100&CRISTA<=100', 'CFD>100&CRISTA>100'])
    count_df.to_csv(sys.argv[2]+f'count_list_top10000_union_with_{max_bulges}_bulges.tsv',
                    sep='\t', na_rep='NA', index=False)

    plt.figure(figsize=(5, 5))
    # plot = sns.JointGrid(data=final_df, x='CFD_Rank',
    #                      y='CRISTA_Rank', hue='Bulge_count', marginal_ticks=True, palette=palette)
    # plot.plot_joint(sns.scatterplot, alpha=0.5, rasterized=True, legend=False)
    # plot.plot_marginals(sns.histplot)
    sns.scatterplot(data=final_df, x='CFD_Rank', y='CRISTA_Rank',
                    hue='Bulge_count', rasterized=True, palette=palette, alpha=0.5, legend=False, linewidth=0)
    plt.hlines(100, 0, 10000)
    plt.vlines(100, 0, 10000)
    plt.yscale('log')
    plt.xscale('log')
    # plot.ax_joint.axvline(x=100)
    # plot.ax_joint.axhline(y=100)
    # plot.ax_joint.set_xscale('log')
    # plot.ax_joint.set_yscale('log')
    # # plot.ax_joint.invert_yaxis()
    # # plot.ax_joint.invert_xaxis()
    # plot.set_axis_labels('CFD Rank', 'CRISTA Rank')

    plt.tight_layout()
    plt.savefig(
        sys.argv[2]+f'scatter_rank_CFDvCRISTA_top10000_union_with_{max_bulges}_bulges.pdf', dpi=300)
    plt.clf()
    plt.close('all')

    plt.figure(figsize=(5, 5))
    sns.scatterplot(data=final_df, x="CFD_score_(highest_CFD)", y='CRISTA_score_(highest_CRISTA)',
                    hue='Bulge_count', rasterized=True, palette=palette, alpha=0.5, legend=False, linewidth=0)

    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    plt.tight_layout()
    plt.savefig(
        sys.argv[2]+f'correlation_CFDvCRISTA_top10000_union_with_{max_bulges}_bulges.pdf', dpi=300)
    plt.clf()
    plt.close('all')

    # pearson corr
    corr = stats.pearsonr(final_df['CFD_score_(highest_CFD)'],
                          final_df['CRISTA_score_(highest_CRISTA)'])
    outfile = open(
        sys.argv[2]+f'pearson_corr_scores_with_{max_bulges}_bulges.txt', 'w')
    outfile.write(str(corr[0])+'\t'+str(corr[1]))
    outfile.close()

    corr = stats.pearsonr(final_df['CFD_Rank'],
                          final_df['CRISTA_Rank'])
    outfile = open(
        sys.argv[2]+f'pearson_corr_ranks_with_{max_bulges}_bulges.txt', 'w')
    outfile.write(str(corr[0])+'\t'+str(corr[1]))
    outfile.close()

    # sperman corr
    corr = stats.spearmanr(final_df['CFD_score_(highest_CFD)'],
                           final_df['CRISTA_score_(highest_CRISTA)'])

    outfile = open(
        sys.argv[2]+f'spearman_corr_scores_with_{max_bulges}_bulges.txt', 'w')
    outfile.write(str(corr[0])+'\t'+str(corr[1]))
    outfile.close()

    corr = stats.spearmanr(final_df['CFD_Rank'],
                           final_df['CRISTA_Rank'])
    outfile = open(
        sys.argv[2]+f'spearman_corr_ranks_with_{max_bulges}_bulges.txt', 'w')
    outfile.write(str(corr[0])+'\t'+str(corr[1]))
    outfile.close()


print('start processing')
original_df = pd.read_csv(sys.argv[1], sep="\t", index_col=False,
                          na_values=['n'], usecols=['Spacer+PAM', 'Chromosome', 'Bulges_(highest_CFD)', 'Bulges_(highest_CRISTA)', 'CFD_score_(highest_CFD)', 'CRISTA_score_(highest_CRISTA)'])

# select the max number of bulges allowed in targets
max_bulges = 0
print('bulge 0')
# filter out targets with bulges != max_bulges
filter_bulges = original_df.loc[(
    original_df['Bulges_(highest_CFD)'] == max_bulges) & (original_df['Bulges_(highest_CRISTA)'] == max_bulges)]
filter_bulges['Bulge_count'] = '0'
plot_correlation(filter_bulges, max_bulges)

# bulge <=1
max_bulges = 1
print('bulge 1')
# filter out targets with bulges > max_bulges
filter_bulges = original_df.loc[(
    original_df['Bulges_(highest_CFD)'] == max_bulges) | (original_df['Bulges_(highest_CRISTA)'] == max_bulges)]
filter_bulges = filter_bulges.loc[(filter_bulges['Bulges_(highest_CFD)'] != 2) & (
    filter_bulges['Bulges_(highest_CRISTA)'] != 2)]
filter_bulges['Bulge_count'] = filter_bulges.apply(bulge_color, axis=1)
plot_correlation(filter_bulges, max_bulges)

# bulge<=2
max_bulges = 2
print('bulge 2')
# filter out targets with bulges > max_bulges
filter_bulges = original_df.loc[(
    original_df['Bulges_(highest_CFD)'] == max_bulges) | (original_df['Bulges_(highest_CRISTA)'] == max_bulges)]
filter_bulges['Bulge_count'] = filter_bulges.apply(bulge_color, axis=1)
plot_correlation(filter_bulges, max_bulges)
