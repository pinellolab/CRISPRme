"""
Linda Lin
1/15/2021

Input: CRISPRme best file (top 1000 sites or more, without on-target)
Output: top OT plots for CRISPRme manuscript

(Note: will get runtime warnings but all are fine to ignore)
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
#import seaborn as sns
import sys
import matplotlib
import warnings

# ignore all warnings
warnings.filterwarnings("ignore")

# set matplotlib to not use X11 server
matplotlib.use('Agg')


def plot_with_MMvBUL(df, out_folder, guide):
    # Remove targets ref with mm+bul<=1 for on-targets and on-targets variant-induced
    df = df.loc[df["Mismatches+bulges_(fewest_mm+b)"] > 1]

    # new col to store the scoring value for non-SpCas9 targets
    df['Mismatches+bulges_REF_(fewest_mm+b)'] = 0
    df['Mismatches+bulges_ALT_(fewest_mm+b)'] = 0

    # if col is alt calculate score for ref and alt, if ref skip
    for index in df.index:
        if df.loc[index, 'REF/ALT_origin_(fewest_mm+b)'] == 'alt':
            refTarget = str(
                df.loc[index, 'Aligned_protospacer+PAM_REF_(fewest_mm+b)'])
            countMM = 0
            for nt in refTarget:
                if nt.islower():
                    countMM += 1
            df.loc[index, 'Mismatches+bulges_REF_(fewest_mm+b)'] = countMM + \
                int(df.loc[index, 'Bulges_(fewest_mm+b)'])
            df.loc[index, 'Mismatches+bulges_ALT_(fewest_mm+b)'] = df.loc[index,
                                                                          'Mismatches+bulges_(fewest_mm+b)']
        else:
            df.loc[index, 'Mismatches+bulges_REF_(fewest_mm+b)'] = df.loc[index,
                                                                          'Mismatches+bulges_(fewest_mm+b)']
            df.loc[index, 'Mismatches+bulges_ALT_(fewest_mm+b)'] = df.loc[index,
                                                                          'Mismatches+bulges_(fewest_mm+b)']

    # sort in order to have highest REF mm+bul on top
    df.sort_values('Mismatches+bulges_(fewest_mm+b)',
                   ascending=True, inplace=True)
    # top1000 targets
    df = df.head(100)
    # Make index column that numbers the OTs starting from 1
    df.reset_index(inplace=True)
    index_count = 1
    for index in df.index:
        df.loc[index, 'index'] = index_count
        index_count += 1
    # If prim_AF = 'n', then it's a ref-nominated site, so we enter a fake numerical AF
    # This will cause a warning of invalid sqrt later on, but that's fine to ignore
    # df["prim_AF"] = df["prim_AF"].fillna(-1)
    df["Variant_MAF_(fewest_mm+b)"] = df["Variant_MAF_(fewest_mm+b)"].fillna(-1)

    # If multiple AFs (haplotype with multiple SNPs), take min AF
    # Approximation until we have haplotype frequencies
    df["AF"] = df["Variant_MAF_(fewest_mm+b)"].astype(str).str.split(',')
    df["AF"] = df["AF"].apply(lambda x: min(x))
    df["AF"] = pd.to_numeric(df["AF"])

    # Adjustments for plotting purposes
    # so haplotypes that got rounded down to AF = 0 (min AF = 0.01) still appear in the plot
    df["plot_AF"] = df["AF"] + 0.001
    df["plot_AF"] *= 1000  # make points larger
    df["plot_AF"] = np.sqrt(df["plot_AF"])  # so size increase is linear

    # Calculate ref_AF as (1 – alt_AF)
    # Not precisely correct because there can be other non-ref haplotypes, but approximation should be accurate in most cases
    df["ref_AF"] = 1 - df["AF"]
    df["ref_AF"] *= 1000  # make points larger
    df["ref_AF"] = np.sqrt(df["ref_AF"])  # so size increase is linear

    # Transparent colors
    transparent_red = mcolors.colorConverter.to_rgba("red", alpha=0.5)
    transparent_blue = mcolors.colorConverter.to_rgba("blue", alpha=0.5)
    transparent_gray = mcolors.colorConverter.to_rgba("gray", alpha=0.5)

    # # Size legend
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')

    """
    Log, ref/alt, top 1000: for main text
    """
    # matplotlib plot settings
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams["figure.figsize"] = 7.5, 2.25
    plt.rcParams.update({'font.size': 7})
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Plot data
    ax = df.plot.scatter(x="index", y="Mismatches+bulges_REF_(fewest_mm+b)",
                         s="ref_AF", c=transparent_red, zorder=1)

    # ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)", s="ref_AF", c=transparent_red, zorder=1, ax=ax)
    df.plot.scatter(x="index", y="Mismatches+bulges_ALT_(fewest_mm+b)",
                    s="plot_AF", c=transparent_blue, zorder=2, ax=ax)

    ax.set_xscale("log")
    # plt.title("Top CRISPRme-identified sites for sgRNA 1617")
    plt.xlabel("Candidate off-target site")
    plt.ylabel("Mismatches+Bulges")

    # Boundaries
    plt.xlim(xmin=0.9, xmax=100)
    # plt.ylim(ymin=0, ymax=1)

    # Arrows
    for x, y, z in zip(df["index"], df["Mismatches+bulges_REF_(fewest_mm+b)"], df["Mismatches+bulges_ALT_(fewest_mm+b)"]-df["Mismatches+bulges_REF_(fewest_mm+b)"]):
        plt.arrow(x, y+0.02, 0, z-0.04, color='gray', head_width=(x*(10**0.005-10**(-0.005))),
                  head_length=0.02, length_includes_head=True, zorder=0, alpha=0.5)
        # +/- to avoid overlap of arrow w/ points, head_width calculated to remain constant despite log scale of x-axis

    # Size legend
    plt.gca().add_artist(plt.legend(
        handles=[s1, s01, s001], title="Allele frequency", ncol=3, loc=9))

    # Color legend
    red = mpatches.Patch(color=transparent_red, label="Reference")
    blue = mpatches.Patch(color=transparent_blue, label="Alternative")
    plt.legend(handles=[red, blue])

    # Save
    plt.tight_layout()
    plt.savefig(
        out_folder+f"CRISPRme_fewest_top_1000_log_for_main_text_{guide}.png")
    plt.clf()


def plot_with_CRISTA_score(df, out_folder, guide):
    # Remove targets with mm+bul<=1 since they are probably on-target introduced by variants
    df = df.loc[df["Mismatches+bulges_(highest_CRISTA)"] > 1]
    # sort values to have highest scored target on top
    df.sort_values('CRISTA_score_(highest_CRISTA)',
                   ascending=False, inplace=True)
    # keep top1000 targets
    df = df.head(100)
    # Make index column that numbers the OTs starting from 1
    df.reset_index(inplace=True)
    index_count = 1
    for index in df.index:
        df.loc[index, 'index'] = index_count
        index_count += 1

    # If prim_AF = 'n', then it's a ref-nominated site, so we enter a fake numerical AF
    # This will cause a warning of invalid sqrt later on, but that's fine to ignore
    df["Variant_MAF_(highest_CRISTA)"] = df["Variant_MAF_(highest_CRISTA)"].fillna(-1)

    # If multiple AFs (haplotype with multiple SNPs), take min AF
    # Approximation until we have haplotype frequencies
    df["AF"] = df["Variant_MAF_(highest_CRISTA)"].astype(str).str.split(',')
    df["AF"] = df["AF"].apply(lambda x: min(x))
    df["AF"] = pd.to_numeric(df["AF"])

    # Adjustments for plotting purposes
    # so haplotypes that got rounded down to AF = 0 (min AF = 0.01) still appear in the plot
    df["plot_AF"] = df["AF"] + 0.001
    df["plot_AF"] *= 1000  # make points larger
    df["plot_AF"] = np.sqrt(df["plot_AF"])  # so size increase is linear

    # Calculate ref_AF as (1 – alt_AF)
    # Not precisely correct because there can be other non-ref haplotypes, but approximation should be accurate in most cases
    df["ref_AF"] = 1 - df["AF"]
    df["ref_AF"] *= 1000  # make points larger
    df["ref_AF"] = np.sqrt(df["ref_AF"])  # so size increase is linear

    # Transparent colors
    transparent_red = mcolors.colorConverter.to_rgba("red", alpha=0.5)
    transparent_blue = mcolors.colorConverter.to_rgba("blue", alpha=0.5)
    transparent_gray = mcolors.colorConverter.to_rgba("gray", alpha=0.5)

    # # Size legend
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')

    """
    Log, ref/alt, top 1000: for main text
    """
    # matplotlib plot settings
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams["figure.figsize"] = 7.5, 2.25
    plt.rcParams.update({'font.size': 7})
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Plot data CFD SCORE
    ax = df.plot.scatter(x="index", y="CRISTA_score_REF_(highest_CRISTA)",
                         s="ref_AF", c=transparent_red, zorder=1)
    # ax = df.plot.scatter(x="index", y="highest_CRISTA_score(ref)", s="ref_AF", c=transparent_red, zorder=1, ax=ax)
    df.plot.scatter(x="index", y="CRISTA_score_ALT_(highest_CRISTA)",
                    s="plot_AF", c=transparent_blue, zorder=2, ax=ax)
    ax.set_xscale("log")

    plt.xlabel("Candidate off-target site")
    plt.ylabel("CRISTA score")

    # Boundaries
    plt.xlim(xmin=0.9, xmax=100)
    plt.ylim(ymin=0, ymax=1)

    # Arrows
    for x, y, z in zip(df["index"], df["CRISTA_score_REF_(highest_CRISTA)"], df["CRISTA_score_ALT_(highest_CRISTA)"]-df["CRISTA_score_REF_(highest_CRISTA)"]):
        plt.arrow(x, y+0.02, 0, z-0.04, color='gray', head_width=(x*(10**0.005-10**(-0.005))),
                  head_length=0.02, length_includes_head=True, zorder=0, alpha=0.5)
        # +/- to avoid overlap of arrow w/ points, head_width calculated to remain constant despite log scale of x-axis

    # Size legend
    plt.gca().add_artist(plt.legend(
        handles=[s1, s01, s001], title="Allele frequency", ncol=3, loc=9))

    # Color legend
    red = mpatches.Patch(color=transparent_red, label="Reference")
    blue = mpatches.Patch(color=transparent_blue, label="Alternative")
    plt.legend(handles=[red, blue])

    # Save
    plt.tight_layout()
    plt.savefig(
        out_folder+f"CRISPRme_CRISTA_top_1000_log_for_main_text_{guide}.png")
    plt.clf()


def plot_with_CFD_score(df, out_folder, guide):
    # Remove targets with mm+bul<=1 since they are probably on-target introduced by variants
    df = df.loc[df["Mismatches+bulges_(highest_CFD)"] > 1]

    # sort values to have highest scored target on top
    df.sort_values('CFD_score_(highest_CFD)',
                   ascending=False, inplace=True)
    # keep top1000 targets
    df = df.head(100)
    # Make index column that numbers the OTs starting from 1
    df.reset_index(inplace=True)
    index_count = 1
    for index in df.index:
        df.loc[index, 'index'] = index_count
        index_count += 1

    # If prim_AF = 'n', then it's a ref-nominated site, so we enter a fake numerical AF
    # This will cause a warning of invalid sqrt later on, but that's fine to ignore
    df["Variant_MAF_(highest_CFD)"] = df["Variant_MAF_(highest_CFD)"].fillna(-1)

    # If multiple AFs (haplotype with multiple SNPs), take min AF
    # Approximation until we have haplotype frequencies
    df["AF"] = df["Variant_MAF_(highest_CFD)"].astype(str).str.split(',')
    df["AF"] = df["AF"].apply(lambda x: min(x))
    df["AF"] = pd.to_numeric(df["AF"])

    # Adjustments for plotting purposes
    # so haplotypes that got rounded down to AF = 0 (min AF = 0.01) still appear in the plot
    df["plot_AF"] = df["AF"] + 0.001
    df["plot_AF"] *= 1000  # make points larger
    df["plot_AF"] = np.sqrt(df["plot_AF"])  # so size increase is linear

    # Calculate ref_AF as (1 – alt_AF)
    # Not precisely correct because there can be other non-ref haplotypes, but approximation should be accurate in most cases
    df["ref_AF"] = 1 - df["AF"]
    df["ref_AF"] *= 1000  # make points larger
    df["ref_AF"] = np.sqrt(df["ref_AF"])  # so size increase is linear

    # Transparent colors
    transparent_red = mcolors.colorConverter.to_rgba("red", alpha=0.5)
    transparent_blue = mcolors.colorConverter.to_rgba("blue", alpha=0.5)
    transparent_gray = mcolors.colorConverter.to_rgba("gray", alpha=0.5)

    """
    # Linear, annotated, top 100 (for supplement)
    # """

    # # Size legend
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')

    """
    Log, ref/alt, top 1000: for main text
    """
    # matplotlib plot settings
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams["figure.figsize"] = 7.5, 2.25
    plt.rcParams.update({'font.size': 7})
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Plot data CFD SCORE
    ax = df.plot.scatter(x="index", y="CFD_score_REF_(highest_CFD)",
                         s="ref_AF", c=transparent_red, zorder=1)
    # ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)", s="ref_AF", c=transparent_red, zorder=1, ax=ax)
    df.plot.scatter(x="index", y="CFD_score_ALT_(highest_CFD)",
                    s="plot_AF", c=transparent_blue, zorder=2, ax=ax)
    ax.set_xscale("log")

    plt.xlabel("Candidate off-target site")
    plt.ylabel("CFD score")

    # Boundaries
    plt.xlim(xmin=0.9, xmax=100)
    plt.ylim(ymin=0, ymax=1)

    # Arrows
    for x, y, z in zip(df["index"], df["CFD_score_REF_(highest_CFD)"], df["CFD_score_ALT_(highest_CFD)"]-df["CFD_score_REF_(highest_CFD)"]):
        plt.arrow(x, y+0.02, 0, z-0.04, color='gray', head_width=(x*(10**0.005-10**(-0.005))),
                  head_length=0.02, length_includes_head=True, zorder=0, alpha=0.5)
        # +/- to avoid overlap of arrow w/ points, head_width calculated to remain constant despite log scale of x-axis

    # Size legend
    plt.gca().add_artist(plt.legend(
        handles=[s1, s01, s001], title="Allele frequency", ncol=3, loc=9))

    # Color legend
    red = mpatches.Patch(color=transparent_red, label="Reference")
    blue = mpatches.Patch(color=transparent_blue, label="Alternative")
    plt.legend(handles=[red, blue])

    # Save
    plt.tight_layout()
    plt.savefig(
        out_folder+f"CRISPRme_CFD_top_1000_log_for_main_text_{guide}.png")
    plt.clf()


# Read file
df_guide = pd.read_csv(sys.argv[1], sep="\t",
                       index_col=False, na_values=['n'])
out_folder = sys.argv[2]
guide = sys.argv[3]

# for guide in df['Spacer+PAM'].unique():
# reset df temp after every process to avoid memory problems
# df_guide = df.loc[df["Spacer+PAM"] == guide]
# takes df in input and produces the correlated plot, guide is used to save with correct name
plot_with_CFD_score(df_guide, out_folder, guide)
plot_with_CRISTA_score(df_guide, out_folder, guide)
plot_with_MMvBUL(df_guide, out_folder, guide)
