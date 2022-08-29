# from adjustText import adjust_text
import warnings
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# from matplotlib.image import BboxImage
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys
import os
import pandas as pd
import numpy as np
import math
import seaborn as sns

# ignore all warnings
warnings.filterwarnings("ignore")
# set matplotlib to not use X11 server
matplotlib.use('Agg')

# set matplotlib to print in pdf editable format
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.style.use('seaborn-poster')
sns.set_context("poster", font_scale=0.8)


gene_target_dict = {
    "AAAGGCTGCTGATGACACCTNNN": "TTR",
    "CCCAGAAGGGGACAGTAAGANNN": "CCR5",
    "CCCGCACCTTGGCGCAGCGGNNN": "PCSK9",
    "CTAACAGTTGCTTTTATCACNNN": "BCL11A",
    "CTTGCCCCACAGGGCAGTAANNN": "HBB",
    "CTTGTCAAGGCTATTGGTCANNN": "HGB1",
    "GAGTCCGAGCAGAAGAAGAANNN": "EMX1",
    "GGAATCCCTTCTGCAGCACCNNN": "FANFC",
    "GGAGAATGACGAGTGGACCCNNN": "TRBC",
    "GGCGCCCTGGCCAGTCGTCTNNN": "PDCD1",
    "GGGTGGGAAAATAGACTAATNNN": "HBB",
    "TCACTATGCTGCCGCCCAGTNNN": "CCR5",
    "TGTGCTAGACATGAGGTCTANNN": "TRAC1_TRAC2",
    "TTTATCACAGGCTCCAGGAANNN": "BCL11A"
}

iupac_code_set = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "r": {"A", "G"},
    "y": {"C", "T"},
    "s": {"G", "C"},
    "w": {"A", "T"},
    "k": {"G", "T"},
    "m": {"A", "C"},
    "b": {"C", "G", "T"},
    "d": {"A", "G", "T"},
    "h": {"A", "C", "T"},
    "v": {"A", "C", "G"},
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
    "a": {"a"},
    "t": {"t"},
    "c": {"c"},
    "g": {"g"},
    'N': {'A', 'T', 'G', 'C'}
}


def crisprme_plot_MMvBUL(df, guide, out_folder, max_mm_bul_value, pam_first_nucleotide, pam_len):

    # new col to store the scoring value for non-SpCas9 targets
    df['Mismatches+bulges_(fewest_mm+b)_REF'] = 0
    df['Mismatches+bulges_(fewest_mm+b)_DELTA'] = 0
    df['Mismatches+bulges_(fewest_mm+b)_REF_NORM'] = 0
    df['Mismatches+bulges_(fewest_mm+b)_ALT_NORM'] = 0

    # if col is alt calculate score for ref and alt, if ref calculate only ref
    for index in df.index:
        if df.loc[index, 'REF/ALT_origin_(fewest_mm+b)'] == 'alt':
            countMM = 0  # ref MM+bul counter (including PAM)

            # count MM in ref target wrt alt target
            refTarget = str(
                df.loc[index, 'Aligned_protospacer+PAM_REF_(fewest_mm+b)'])
            altTarget = str(
                df.loc[index, 'Aligned_protospacer+PAM_ALT_(fewest_mm+b)'])

            for nt in refTarget:
                if nt.islower():
                    countMM += 1

            # print('ref', refTarget, 'alt', altTarget)
            # check if PAM creation added MM
            if pam_first_nucleotide == 0:  # pat at start
                # print('ref PAM', refTarget[:pam_len])
                # print('alt PAM', altTarget[:pam_len])
                for pos in range(pam_len):
                    if refTarget[:pam_len][pos] != altTarget[:pam_len][pos]:
                        countMM += 1
            else:  # pam at end
                # print('ref PAM', refTarget[-pam_len:])
                # print('alt PAM', altTarget[-pam_len:])
                for pos in range(pam_len):
                    if refTarget[-pam_len:][pos] != altTarget[-pam_len:][pos]:
                        countMM += 1

            # update mm+bul ref with new count
            df.loc[index, 'Mismatches+bulges_(fewest_mm+b)_REF'] = int(
                df.loc[index, 'Bulges_(fewest_mm+b)'])+countMM
        else:
            # if ref no distance between REF and ALT sequence (ref/alt)
            df.loc[index, 'Mismatches+bulges_(fewest_mm+b)_REF'] = df.loc[index,
                                                                          'Mismatches+bulges_(fewest_mm+b)']
        df.loc[index, 'Mismatches+bulges_(fewest_mm+b)_DELTA'] = int(df.loc[index,
                                                                            'Mismatches+bulges_(fewest_mm+b)_REF']) - int(df.loc[index, 'Mismatches+bulges_(fewest_mm+b)'])

    df['Mismatches+bulges_(fewest_mm+b)_REF_NORM'] = df['Mismatches+bulges_(fewest_mm+b)_REF'] / \
        df['Mismatches+bulges_(fewest_mm+b)_REF'].max()
    # df.loc[df['Mismatches+bulges_(fewest_mm+b)_REF_NORM']
    #        > 1, 'Mismatches+bulges_(fewest_mm+b)_REF_NORM'] = 1
    df['Mismatches+bulges_(fewest_mm+b)_ALT_NORM'] = df['Mismatches+bulges_(fewest_mm+b)'] / \
        df['Mismatches+bulges_(fewest_mm+b)_REF'].max()

    # sort the df
    df.sort_values(['Mismatches+bulges_(fewest_mm+b)', 'Mismatches+bulges_(fewest_mm+b)_DELTA', 'Variant_MAF_(fewest_mm+b)'],
                   inplace=True, ascending=True)

    # Make index column that numbers the OTs starting from 1 after sorting
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
    transparent_black = mcolors.colorConverter.to_rgba("black", alpha=0.5)

    # # Size legend
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')
    s0001 = mlines.Line2D([], [], marker='o', label='0.001', linestyle='None',
                          markersize=math.sqrt(math.sqrt((0.001+0.001)*1000)), color='black')

    """
    Log, ref/alt, top 1000: for main text
    """
    # matplotlib plot settings
    plt.rcParams["figure.dpi"] = 600
    # plt.rcParams["figure.figsize"] = 7.5, 2.25
    plt.rcParams.update({'font.size': 7})
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    plt.figure(figsize=(10, 4))
    ax1 = plt.subplot(7, 1, (1, 5))
    # ax1.set_position([0.125, 0.4, 0.9, 0.9])
    # print('ax1-position', ax1.get_position())
    # Plot data
    ax1.scatter(x=df['index'], y=df['Mismatches+bulges_(fewest_mm+b)_REF_NORM'],
                s=df['ref_AF'], c=transparent_red)
    ax1.scatter(x=df['index'], y=df['Mismatches+bulges_(fewest_mm+b)_ALT_NORM'],
                s=df['plot_AF'], c=transparent_blue)
    # Plot data
    # ax = df.plot.scatter(x="index", y="Mismatches+bulges_(fewest_mm+b)_REF",
    #                      s="ref_AF", c=transparent_red, zorder=1)

    # # ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)", s="ref_AF", c=transparent_red, zorder=1, ax=ax)
    # df.plot.scatter(x="index", y="Mismatches+bulges_(fewest_mm+b)",
    #                 s="plot_AF", c=transparent_blue, zorder=2, ax=ax)

    # texts = list()
    # annotate the first 20 alternative targets with gene name (skip TSG)
    df_coordinates_CDS = df.loc[(df['REF/ALT_origin_(fewest_mm+b)']
                                 != 'ref') & (df['Annotation_GENCODE']).str.contains('CDS')]

    # counter = 0
    # for index, row in df_coordinates.iterrows():
    #     texts.append([row["index"], row["Mismatches+bulges_(fewest_mm+b)"], str(
    #         row["Annotation_closest_gene_name"]), 'blue'])
    #     counter += 1
    #     if counter >= 20:
    #         break

    # annotate alternative target if the gene is associated with tumor suppression activity and CDS
    # df_coordinates_TSG = df.loc[(df['Gene_description'].notnull()) & (
    # df['REF/ALT_origin_(fewest_mm+b)'] != 'ref') & (df['Annotation_GENCODE']).str.contains('CDS')]

    # counter = 0
    # for index, row in df_coordinates.iterrows():
    #     texts.append([row["index"], row["Mismatches+bulges_(fewest_mm+b)"], str(
    #         row["Annotation_closest_gene_name"]), 'red'])
    #     counter += 1
    #     if counter >= 20:
    #         break

    # annotate encode targets
    df_coordinates_ENCODE = df.loc[(
        df['REF/ALT_origin_(fewest_mm+b)'] != 'ref') & (df['Annotation_ENCODE']).notnull()]

    # counter = 0
    # for index, row in df_coordinates.iterrows():
    #     texts.append([row["index"], row["Mismatches+bulges_(fewest_mm+b)"], str(
    #         row["Annotation_ENCODE"]), 'black'])
    #     counter += 1
    #     if counter >= 20:
    #         break

    # sort text to order by x
    # texts.sort(key=lambda x: x[0])

    # plot labels and then adjust them to avoid overlap
    # plotted_text = [plt.text(texts[i][0], texts[i][1], texts[i][2],
    #                          ha='center', va='center', color=texts[i][3], fontsize=5) for i in range(len(texts))]
    # adjust_text(plotted_text)

    # plt.xlabel("Candidate off-target site")
    plt.ylabel("Mismatches + Bulges Normalized")

    # Boundaries
    ax1.margins(0.05)
    plt.xticks([1, 20, 40, 60, 80, 100])
    # plt.xlim(1,100)
    # plt.ylim(1, df['Mismatches+bulges_(fewest_mm+b)_REF'].max())

    # Arrows
    for x, y, z in zip(df["index"], df["Mismatches+bulges_(fewest_mm+b)_ALT_NORM"], df["Mismatches+bulges_(fewest_mm+b)_REF_NORM"]-df["Mismatches+bulges_(fewest_mm+b)_ALT_NORM"]):
        if z != 0:
            plt.arrow(x, y, 0, z-0.04, color='gray', zorder=0,
                      alpha=0.5, head_width=0.5, head_length=0.02)

    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')
    s0001 = mlines.Line2D([], [], marker='o', label='0.001', linestyle='None',
                          markersize=math.sqrt(math.sqrt((0.001+0.001)*1000)), color='black')
    s1_red = mlines.Line2D([], [], marker='o', label='Reference', linestyle='None',
                           markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_red)
    s1_blue = mlines.Line2D([], [], marker='o', label='Alternative', linestyle='None',
                            markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_blue)
    green_PLS = mpatches.Patch(color='tab:green', label="PLS")
    blue_pELS = mpatches.Patch(color='tab:blue', label="pELS")
    gray_dELS = mpatches.Patch(color='tab:gray', label="dELS")
    purple_CTCFonly = mpatches.Patch(color='tab:purple', label="CTCF-only")
    olive_DNaseH3K4me3 = mpatches.Patch(
        color='tab:olive', label="DNase-H3K4me3")
    cyan_CDS = mpatches.Patch(color='tab:cyan', label="CDS")
    red_TSG = mpatches.Patch(color='tab:red', label="TSG")

    # legend for allele frequency
    plt.gca().add_artist(plt.legend(handles=[s1, s01, s001, s0001], title='Allele frequency', bbox_to_anchor=(
        0.01, 1.03), loc='lower left', borderaxespad=0, ncol=4))
    # legend for target origin
    plt.gca().add_artist(plt.legend(handles=[s1_red, s1_blue], title='Target genome', bbox_to_anchor=(
        0.30, 1.03), loc='lower left', borderaxespad=0, ncol=3))
    # legend for encode annotation
    plt.gca().add_artist(plt.legend(handles=[green_PLS, cyan_CDS, blue_pELS, red_TSG, gray_dELS, purple_CTCFonly, olive_DNaseH3K4me3],
                                    loc='lower left', bbox_to_anchor=(0.54, 1.02), title='Annotations', borderaxespad=0, ncol=5))

    # ENCODE BARPLOT
    ax2 = plt.subplot(7, 1, 6, sharex=ax1)
    # print(list(df_coordinates_ENCODE['index']))
    bar_list = [0]*100
    bar_range = np.arange(1, 101, 1)
    bar_color = ['white']*100
    for index, row in df_coordinates_ENCODE.iterrows():
        bar_list[index] = 0.01
        if 'PLS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:green'
        elif 'pELS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:blue'
        elif 'dELS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:gray'
        elif 'CTCF-only' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:purple'
        elif 'DNase-H3K4me3' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:olive'
    # print(bar_list)
    # print(df['index'])
    # print(df_coordinates_ENCODE['index'])
    ax2.bar(bar_range, bar_list, width=0.3, align='center', color=bar_color)
    # plt.xticks([])
    plt.yticks([])
    # ax2.margins(0.05)
    # plt.xlim(1, 100)
    plt.ylabel('ALT cCRE')

    # CDS AND TSG BARPLOT
    ax3 = plt.subplot(7, 1, 7, sharex=ax1)
    bar_list = [0]*100
    bar_range = np.arange(1, 101, 1)
    bar_color = ['white']*100
    for index, row in df_coordinates_CDS.iterrows():
        bar_list[index] = 0.01
        bar_color[index] = 'tab:cyan'
        # print(str(row['Gene_description']))
        if str(row['Gene_description']) != 'nan':
            bar_color[index] = 'tab:red'
    ax3.bar(bar_range, bar_list, width=0.3, align='center', color=bar_color)
    # plt.xticks([])
    plt.yticks([])
    # ax3.margins(0.05)
    # plt.xlim(1, 100)
    plt.ylabel('ALT CDS')

    # columns to drop from df
    columns_to_drop = list()
    for column in list(df.columns):
        if 'highest_CFD' in column:
            columns_to_drop.append(column)
    columns_to_drop.append('index')
    df.drop(columns_to_drop, axis=1, inplace=True)

    # save df used in the analysis
    df.to_csv(out_folder+guide+'_processed_data.tsv',
              sep='\t', na_rep='NA', index=False)
    # Save
    plt.tight_layout()
    plt.savefig(
        out_folder+f"CRISPRme_top_1000_log_for_main_text_{guide}.pdf", transparent=True)
    plt.clf()
    plt.close('all')


def crisprme_plot_CFD(title, df, guide, out_folder):

    # Make index column that numbers the OTs starting from 0
    df.reset_index(inplace=True)
    index_count = 1
    for index in df.index:
        df.loc[index, 'index'] = index_count
        index_count += 1

    # If prim_AF = 'n', then it's a ref-nominated site, so we enter a fake numerical AF
    # This will cause a warning of invalid sqrt later on, but that's fine to ignore
    # df["prim_AF"] = df["prim_AF"].fillna(-1)
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
    transparent_black = mcolors.colorConverter.to_rgba("black", alpha=0.5)

    """
    Log, ref/alt, top 1000: for main text
    """
    # matplotlib plot settings
    # plt.rcParams["figure.dpi"] = 600
    # # plt.rcParams["figure.figsize"] = 7.5, 2.25
    # plt.rcParams.update({'font.size': 7})
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.rcParams['ps.fonttype'] = 42

    # plt.figure(figsize=(25, 22))
    plt.figure(figsize=(8.5, 10))
    # plt.figure()
    ax1 = plt.subplot(7, 1, (1, 5))
    # ax1.set_position([0.125, 0.4, 0.9, 0.9])
    # print('ax1-position', ax1.get_position())
    # Plot data
    ax1.scatter(x=df['index'], y=df['CFD_score_REF_(highest_CFD)'],
                s=df['ref_AF'], c=transparent_red)
    ax1.scatter(x=df['index'], y=df['CFD_score_ALT_(highest_CFD)'],
                s=df['plot_AF'], c=transparent_blue)

    # texts = list()
    # annotate alternative targets with CDS regions
    df_coordinates_CDS = df.loc[(df['REF/ALT_origin_(highest_CFD)']
                                 != 'ref') & (df['Annotation_GENCODE']).str.contains('CDS')]

    # annotate encode targets
    df_coordinates_ENCODE = df.loc[(
        df['REF/ALT_origin_(highest_CFD)'] != 'ref') & (df['Annotation_ENCODE']).notnull()]

    # plt.xlabel("Candidate off-target site")
    plt.title(title)
    plt.ylabel("CFD Score")

    # Boundaries
    ax1.margins(0.05)
    # plt.xlim(1,100)
    plt.xticks([1, 20, 40, 60, 80, 100])
    plt.ylim(0, 1)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    # Arrows
    # head_width=(x*(10**0.005-10**(-0.005)))
    #   plt.arrow(x, y+0.02, 0, z-0.04, color='gray', head_width=0.02, head_length=0.02,
    #           length_includes_head=True, zorder=0, alpha=0.5)
    for x, y, z in zip(df["index"], df["CFD_score_REF_(highest_CFD)"], df["CFD_score_ALT_(highest_CFD)"]-df["CFD_score_REF_(highest_CFD)"]):
        if z != 0:
            plt.arrow(x, y+0.02, 0, z-0.04, color='gray', zorder=0,
                      alpha=0.5, head_width=0.5, head_length=0.02)
        # +/- to avoid overlap of arrow w/ points, head_width calculated to remain constant despite log scale of x-axis

    # Size legend

    # plt.gca().add_artist(plt.legend(handles=[red_TSG, black_ENCODE, blue_CDS], title="Gene Labels", bbox_to_anchor=(
    #     0.35, 1.02), loc='lower left', borderaxespad=0, ncol=3))

    # Color legend
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')
    s0001 = mlines.Line2D([], [], marker='o', label='0.001', linestyle='None',
                          markersize=math.sqrt(math.sqrt((0.001+0.001)*1000)), color='black')
    s1_red = mlines.Line2D([], [], marker='o', label='Reference', linestyle='None',
                           markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_red)
    s1_blue = mlines.Line2D([], [], marker='o', label='Alternative', linestyle='None',
                            markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_blue)
    green_PLS = mpatches.Patch(color='tab:green', label="PLS")
    blue_pELS = mpatches.Patch(color='tab:blue', label="pELS")
    gray_dELS = mpatches.Patch(color='tab:gray', label="dELS")
    purple_CTCFonly = mpatches.Patch(color='tab:purple', label="CTCF-only")
    olive_DNaseH3K4me3 = mpatches.Patch(
        color='tab:olive', label="DNase-H3K4me3")
    cyan_CDS = mpatches.Patch(color='tab:cyan', label="CDS")
    red_TSG = mpatches.Patch(color='tab:red', label="TSG")

    # legend for allele frequency
    # plt.gca().add_artist(plt.legend(handles=[s1, s01, s001, s0001], title='Allele frequency', bbox_to_anchor=(
    #     -0.05, 1.03), loc='lower left', borderaxespad=0, ncol=4))
    # # legend for target origin
    # plt.gca().add_artist(plt.legend(handles=[s1_red, s1_blue], title='Target genome', bbox_to_anchor=(
    #     0.30, 1.03), loc='lower left', borderaxespad=0, ncol=3))
    # # legend for encode annotation
    # plt.gca().add_artist(plt.legend(handles=[green_PLS, cyan_CDS, blue_pELS, red_TSG, gray_dELS, purple_CTCFonly, olive_DNaseH3K4me3],
    #                                 loc='lower left', bbox_to_anchor=(0.54, 1.02), title='Annotations', borderaxespad=0, ncol=5))

    # ENCODE BARPLOT
    ax2 = plt.subplot(7, 1, 6, sharex=ax1)
    # print(list(df_coordinates_ENCODE['index']))
    bar_list = [0]*100
    bar_range = np.arange(1, 101, 1)
    bar_color = ['white']*100
    for index, row in df_coordinates_ENCODE.iterrows():
        bar_list[index] = 0.01
        if 'PLS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:green'
        elif 'pELS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:blue'
        elif 'dELS' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:gray'
        elif 'CTCF-only' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:purple'
        elif 'DNase-H3K4me3' in str(row['Annotation_ENCODE']):
            bar_color[index] = 'tab:olive'
    # print(bar_list)
    # print(df['index'])
    # print(df_coordinates_ENCODE['index'])
    ax2.bar(bar_range, bar_list, width=0.4, align='center', color=bar_color)
    # plt.xticks([])
    plt.yticks([])
    # ax2.margins(0.05)
    # plt.xlim(1, 100)
    plt.ylabel('ALT cCRE')

    # CDS AND TSG BARPLOT
    ax3 = plt.subplot(7, 1, 7, sharex=ax1)
    bar_list = [0]*100
    bar_range = np.arange(1, 101, 1)
    bar_color = ['white']*100
    for index, row in df_coordinates_CDS.iterrows():
        bar_list[index] = 0.01
        bar_color[index] = 'tab:cyan'
        # print(str(row['Gene_description']))
        if str(row['Gene_description']) != 'nan':
            bar_color[index] = 'tab:red'
    ax3.bar(bar_range, bar_list, width=0.4, align='center', color=bar_color)
    # plt.xticks([])
    plt.yticks([])
    # ax3.margins(0.05)
    # plt.xlim(1, 100)
    plt.ylabel('ALT CDS')

    # columns to drop
    columns_to_drop = list()
    for column in list(df.columns):
        if 'fewest_mm+b' in column:
            columns_to_drop.append(column)
    columns_to_drop.append('index')
    df.drop(columns_to_drop, axis=1, inplace=True)

    # save df used to create plots
    df.to_csv(out_folder+guide+'_processed_data.tsv',
              sep='\t', na_rep='NA', index=False)

    # Save
    # plt.tight_layout()
    # plt.subplots_adjust(hspace=0.8)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95,
                        bottom=0.05, hspace=0.8)

    plt.savefig(
        out_folder+f"CRISPRme_top_1000_log_for_main_text_{guide}.pdf", transparent=True, dpi=100)
    plt.clf()
    plt.close('all')


def plot_legend(guide: str, out_folder: str):

    fig = plt.figure(figsize=(8.5, 0.5))
    # Color legend
    transparent_red = mcolors.colorConverter.to_rgba("red", alpha=0.5)
    transparent_blue = mcolors.colorConverter.to_rgba("blue", alpha=0.5)
    s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                       markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
    s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                        markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
    s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                         markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')
    s0001 = mlines.Line2D([], [], marker='o', label='0.001', linestyle='None',
                          markersize=math.sqrt(math.sqrt((0.001+0.001)*1000)), color='black')
    s1_red = mlines.Line2D([], [], marker='o', label='Reference', linestyle='None',
                           markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_red)
    s1_blue = mlines.Line2D([], [], marker='o', label='Alternative', linestyle='None',
                            markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color=transparent_blue)
    green_PLS = mpatches.Patch(color='tab:green', label="PLS")
    blue_pELS = mpatches.Patch(color='tab:blue', label="pELS")
    gray_dELS = mpatches.Patch(color='tab:gray', label="dELS")
    purple_CTCFonly = mpatches.Patch(color='tab:purple', label="CTCF-only")
    olive_DNaseH3K4me3 = mpatches.Patch(
        color='tab:olive', label="DNase-H3K4me3")
    cyan_CDS = mpatches.Patch(color='tab:cyan', label="CDS")
    red_TSG = mpatches.Patch(color='tab:red', label="TSG")

    # legend for allele frequency
    plt.gca().add_artist(plt.legend(handles=[s1, s01, s001, s0001], title='Allele frequency', bbox_to_anchor=(
        0.01, 1.03), loc='lower left', borderaxespad=0, ncol=4))
    # legend for target origin
    plt.gca().add_artist(plt.legend(handles=[s1_red, s1_blue], title='Target genome', bbox_to_anchor=(
        0.30, 1.03), loc='lower left', borderaxespad=0, ncol=3))
    # legend for encode annotation
    plt.gca().add_artist(plt.legend(handles=[green_PLS, cyan_CDS, blue_pELS, red_TSG, gray_dELS, purple_CTCFonly, olive_DNaseH3K4me3],
                                    loc='lower left', bbox_to_anchor=(0.54, 1.02), title='Annotations', borderaxespad=0, ncol=5))

    plt.savefig(
        out_folder+f"Zlegend_{guide}.pdf", transparent=True)
    plt.clf()
    plt.close('all')


def plot_title_figure(guide: str, mm: int, bul: int, cas_protein: str, genome: str, out_folder: str):
    gene_target = gene_target_dict[guide]

    fig = plt.figure(figsize=(8.5, 0.2))
    fig.suptitle(guide+', '+gene_target+', '+cas_protein+', '+genome +
                 ', '+str(mm)+' mismatches'+' + '+str(bul)+' bulges', fontsize=12)
    plt.savefig(
        out_folder+f"Atitle_{guide}.pdf", transparent=True)
    plt.clf()
    plt.close('all')


def extraction_with_CFD(guide, df, out_dir, top_10_list, top_100_list, top_1000_list):
    # select the current analyzed guide and filter for CFD>=0.1
    df_single_guide = df.loc[(df["Spacer+PAM"] == guide)
                             & (df["CFD_score_(highest_CFD)"] >= 0.1)]

    # extract on_target_chr
    try:
        on_target_chr = df_single_guide.loc[(
            df_single_guide['Mismatches+bulges_(highest_CFD)'] == 0)]
        on_target_chr = on_target_chr.iloc[0]['Chromosome']
    except:
        on_target_chr = False

    # Remove targets with total=0 (any on-target)
    df_single_guide.drop(
        df_single_guide.loc[(df_single_guide["Mismatches+bulges_(highest_CFD)"] <= 1)].index, inplace=True)

    # sort values using score and annotation (NA goes last)
    df_single_guide.sort_values(
        ["CFD_score_(highest_CFD)", "Annotation_ENCODE"], inplace=True, ascending=False)

    # group by duplicates on REF sequence to count them
    df_duplicated = df_single_guide.groupby(
        ['Aligned_protospacer+PAM_REF_(highest_CFD)'], as_index=False, sort=False).size()

    # drop duplicate rows, keeping the first occurence (with non-null ENCODE annotation)
    df_single_guide.drop_duplicates(
        subset=['Aligned_protospacer+PAM_REF_(highest_CFD)'], keep='first', inplace=True)
    # rename column size to duplicate regions count
    df_duplicated.rename(
        columns={"size": "Duplicate_regions_count"}, inplace=True)
    # merge columns with identical ref seq to include the count
    df_single_guide = df_single_guide.merge(
        df_duplicated, how='inner', on='Aligned_protospacer+PAM_REF_(highest_CFD)')
    # reduce count of duplicate by one (avoid count the kept sequence)
    df_single_guide['Duplicate_regions_count'] -= 1

    # TOP10 ANALYSIS
    # count number of alternative targets in top10 ordered by CFD
    df_top10 = df_single_guide.head(
        10).loc[(df_single_guide['REF/ALT_origin_(highest_CFD)'] == 'alt')]
    df_top10.sort_values(
        "CFD_score_(highest_CFD)", inplace=True, ascending=False)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top10["AF"] = df_top10["Variant_MAF_(highest_CFD)"].astype(
        str).str.split(',')
    df_top10["AF"] = df_top10["AF"].apply(lambda x: min(x))
    df_top10["AF"] = pd.to_numeric(df_top10["AF"], errors='coerce')
    MAF_mean = df_top10['AF'].mean()
    MAF_std = df_top10['AF'].std()
    MAF_min = df_top10['AF'].min()
    MAF_max = df_top10['AF'].max()
    # save number of alt targets
    top_10_list.append(str(len(df_top10)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top10_temp = df_top10.loc[(
            df_top10['Chromosome'] == on_target_chr)]
        # save count
        top_10_list.append(str(len(df_top10_temp)))
    else:
        top_10_list.append('0')

    # count alt target with cCREs
    df_top10_temp = df_top10.loc[(df_top10['Annotation_ENCODE']).notnull()]
    top_10_list.append(str(len(df_top10_temp)))

    # count alt target with CDS
    df_top10_temp = df_top10.loc[(
        df_top10["Annotation_GENCODE"].str.contains("CDS"))]
    top_10_list.append(str(len(df_top10_temp)))
    # append mean, std and range of MAF for each guide
    top_10_list.append(str(MAF_mean))
    top_10_list.append(str(MAF_std))
    top_10_list.append(str(MAF_min)+'-'+str(MAF_max))

    # TOP100 ANALYSIS
    # count number of alternative targets in top100 ordered by CFD
    df_top100 = df_single_guide.head(
        100).loc[(df_single_guide['REF/ALT_origin_(highest_CFD)'] == 'alt')]
    df_top10.sort_values(
        "CFD_score_(highest_CFD)", inplace=True, ascending=False)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top100["AF"] = df_top100["Variant_MAF_(highest_CFD)"].astype(
        str).str.split(',')
    df_top100["AF"] = df_top100["AF"].apply(lambda x: min(x))
    df_top100["AF"] = pd.to_numeric(df_top100["AF"], errors='coerce')
    MAF_mean = df_top100['AF'].mean()
    MAF_std = df_top100['AF'].std()
    MAF_min = df_top100['AF'].min()
    MAF_max = df_top100['AF'].max()
    # save number of alt targets
    top_100_list.append(str(len(df_top100)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top100_temp = df_top100.loc[(
            df_top100['Chromosome'] == on_target_chr)]
        # save count
        top_100_list.append(str(len(df_top100_temp)))
    else:
        top_100_list.append('0')

    # count alt target with cCREs
    df_top100_temp = df_top100.loc[(df_top100['Annotation_ENCODE']).notnull()]
    top_100_list.append(str(len(df_top100_temp)))

    # count alt target with CDS
    df_top100_temp = df_top100.loc[(
        df_top100["Annotation_GENCODE"].str.contains("CDS"))]
    top_100_list.append(str(len(df_top100_temp)))

    # append mean, std and range of MAF for each guide
    top_100_list.append(str(MAF_mean))
    top_100_list.append(str(MAF_std))
    top_100_list.append(str(MAF_min)+'-'+str(MAF_max))

    # TOP1000 ANALYSIS
    # count number of alternative targets in top10 ordered by CFD
    df_top1000 = df_single_guide.head(
        1000).loc[(df_single_guide['REF/ALT_origin_(highest_CFD)'] == 'alt')]
    df_top1000.sort_values(
        "CFD_score_(highest_CFD)", inplace=True, ascending=True)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top1000["AF"] = df_top1000["Variant_MAF_(highest_CFD)"].astype(
        str).str.split(',')
    df_top1000["AF"] = df_top1000["AF"].apply(lambda x: min(x))
    df_top1000["AF"] = pd.to_numeric(df_top1000["AF"], errors='coerce')
    MAF_mean = df_top1000['AF'].mean()
    MAF_std = df_top1000['AF'].std()
    MAF_min = df_top1000['AF'].min()
    MAF_max = df_top1000['AF'].max()
    # save number of alt targets
    top_1000_list.append(str(len(df_top1000)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top1000_temp = df_top1000.loc[(
            df_top1000['Chromosome'] == on_target_chr)]
        # save count
        top_1000_list.append(str(len(df_top1000_temp)))
    else:
        top_1000_list.append('0')

    # count alt target with cCREs
    df_top1000_temp = df_top1000.loc[(
        df_top1000['Annotation_ENCODE']).notnull()]
    top_1000_list.append(str(len(df_top1000_temp)))

    # count alt target with CDS
    df_top1000_temp = df_top1000.loc[
        df_top1000["Annotation_GENCODE"].str.contains("CDS", na=False)]
    top_1000_list.append(str(len(df_top1000_temp)))

    # append MAF mean,std,range
    top_1000_list.append(str(MAF_mean))
    top_1000_list.append(str(MAF_std))
    top_1000_list.append(str(MAF_min)+'-'+str(MAF_max))

    # TITLE PLOT
    plot_title_figure(guide, 6, 2,
                      'SpCas9', 'hg38+1000G+HGDP', out_dir)

    # PLOT GENERATION
    # extract top100 targets with CFD>=0.1
    title = 'All OTs'
    crisprme_plot_CFD(title, df_single_guide.head(100),
                      guide+'_CFD01', out_dir)

    # extract targets with CDS in gencode annotation
    dff = df_single_guide.loc[df_single_guide["Annotation_GENCODE"].str.contains(
        "CDS", na=False)]
    title = 'CDS'
    crisprme_plot_CFD(title, dff.head(100), guide+'_ZCDS_CFD', out_dir)

    # extract targets with not null encode annotation
    dff = df_single_guide.loc[df_single_guide["Annotation_ENCODE"].notnull()]
    title = 'cCRE'
    crisprme_plot_CFD(title, dff.head(100), guide +
                      '_encode_annotated_CFD', out_dir)

    # extract targets with PAM creation
    dff = df_single_guide.loc[df_single_guide["PAM_creation_(highest_CFD)"].notnull(
    )]
    title = 'PAM creation'
    crisprme_plot_CFD(title, dff.head(100), guide+'_pam_creation_CFD', out_dir)

    # PLOT LEGEND PER PLOT
    # plot_legend(guide, out_dir)


def extraction_with_total(guide, df, out_dir, max_mm_bul_value, pam_first_nucleotide, pam_len, top_10_list, top_100_list, top_1000_list):
    # select the current analyzed guide
    df_single_guide = df.loc[df["Spacer+PAM"] == guide]

    # extract on_target_chr
    try:
        on_target_chr = df_single_guide.loc[(
            df_single_guide['Mismatches+bulges_(fewest_mm+b)'] == 0)]
        on_target_chr = on_target_chr.iloc[0]['Chromosome']
    except:
        on_target_chr = False

    # Remove targets with total=0 (any on-target)
    df_single_guide.drop(
        df_single_guide.loc[(df_single_guide["Mismatches+bulges_(fewest_mm+b)"] == 0)].index, inplace=True)

    # sort for total (ascending) and annotation (NA goes last)
    df_single_guide.sort_values(
        ["Mismatches+bulges_(fewest_mm+b)", "Annotation_ENCODE"], inplace=True, ascending=True)

    # group by duplicates on REF sequence to count them
    df_duplicated = df_single_guide.groupby(
        ['Aligned_protospacer+PAM_REF_(fewest_mm+b)'], as_index=False, sort=False).size()
    # drop duplicate rows, keeping the first occurence (with non-null ENCODE annotation)
    df_single_guide.drop_duplicates(
        subset=['Aligned_protospacer+PAM_REF_(fewest_mm+b)'], keep='first', inplace=True)
    # rename column size to duplicate regions count
    df_duplicated.rename(
        columns={"size": "Duplicate_regions_count"}, inplace=True)
    # merge columns with identical REF seq to include the count
    df_single_guide = df_single_guide.merge(
        df_duplicated, how='inner', on='Aligned_protospacer+PAM_REF_(fewest_mm+b)')
    # reduce count of duplicate by one (avoid count the kept sequence)
    df_single_guide['Duplicate_regions_count'] -= 1

    # TOP10 ANALYSIS
    # count number of alternative targets in top10 ordered by CFD
    df_top10 = df_single_guide.head(
        10).loc[(df_single_guide['REF/ALT_origin_(fewest_mm+b)'] == 'alt')]
    df_top10.sort_values(
        "Mismatches+bulges_(fewest_mm+b)", inplace=True, ascending=True)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top10["AF"] = df_top10["Variant_MAF_(fewest_mm+b)"].astype(
        str).str.split(',')
    df_top10["AF"] = df_top10["AF"].apply(lambda x: min(x))
    df_top10["AF"] = pd.to_numeric(df_top10["AF"], errors='coerce')
    MAF_mean = df_top10['AF'].mean()
    MAF_std = df_top10['AF'].std()
    MAF_min = df_top10['AF'].min()
    MAF_max = df_top10['AF'].max()
    # save number of alt targets
    top_10_list.append(str(len(df_top10)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top10_temp = df_top10.loc[(
            df_top10['Chromosome'] == on_target_chr)]
        # save count
        top_10_list.append(str(len(df_top10_temp)))
    else:
        top_10_list.append('0')

    # count alt target with cCREs
    df_top10_temp = df_top10.loc[(df_top10['Annotation_ENCODE']).notnull()]
    top_10_list.append(str(len(df_top10_temp)))

    # count alt target with CDS
    df_top10_temp = df_top10.loc[(
        df_top10["Annotation_GENCODE"].str.contains("CDS"))]
    top_10_list.append(str(len(df_top10_temp)))

    # append MAF mean,std,range
    top_10_list.append(str(MAF_mean))
    top_10_list.append(str(MAF_std))
    top_10_list.append(str(MAF_min)+'-'+str(MAF_max))

    # TOP100 ANALYSIS
    # count number of alternative targets in top100 ordered by CFD
    df_top100 = df_single_guide.head(
        100).loc[(df_single_guide['REF/ALT_origin_(fewest_mm+b)'] == 'alt')]
    df_top100.sort_values(
        "Mismatches+bulges_(fewest_mm+b)", inplace=True, ascending=True)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top100["AF"] = df_top100["Variant_MAF_(fewest_mm+b)"].astype(
        str).str.split(',')
    df_top100["AF"] = df_top100["AF"].apply(lambda x: min(x))
    df_top100["AF"] = pd.to_numeric(df_top100["AF"], errors='coerce')
    MAF_mean = df_top100['AF'].mean()
    MAF_std = df_top100['AF'].std()
    MAF_min = df_top100['AF'].min()
    MAF_max = df_top100['AF'].max()
    # save number of alt targets
    top_100_list.append(str(len(df_top100)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top100_temp = df_top100.loc[(
            df_top100['Chromosome'] == on_target_chr)]
        # save count
        top_100_list.append(str(len(df_top100_temp)))
    else:
        top_100_list.append('0')

    # count alt target with cCREs
    df_top100_temp = df_top100.loc[(df_top100['Annotation_ENCODE']).notnull()]
    top_100_list.append(str(len(df_top100_temp)))

    # count alt target with CDS
    df_top100_temp = df_top100.loc[(
        df_top100["Annotation_GENCODE"].str.contains("CDS"))]
    top_100_list.append(str(len(df_top100_temp)))

    # append MAF mean,std,range
    top_100_list.append(str(MAF_mean))
    top_100_list.append(str(MAF_std))
    top_100_list.append(str(MAF_min)+'-'+str(MAF_max))

    # TOP1000 ANALYSIS
    # count number of alternative targets in top10 ordered by CFD
    df_top1000 = df_single_guide.head(
        1000).loc[(df_single_guide['REF/ALT_origin_(fewest_mm+b)'] == 'alt')]
    df_top1000.sort_values(
        "Mismatches+bulges_(fewest_mm+b)", inplace=True, ascending=True)
    # create column to contain lower MAF in case of multiple value (2+ variants)
    df_top1000["AF"] = df_top1000["Variant_MAF_(fewest_mm+b)"].astype(
        str).str.split(',')
    df_top1000["AF"] = df_top1000["AF"].apply(lambda x: min(x))
    df_top1000["AF"] = pd.to_numeric(df_top1000["AF"], errors='coerce')
    MAF_mean = df_top1000['AF'].mean()
    MAF_std = df_top1000['AF'].std()
    MAF_min = df_top1000['AF'].min()
    MAF_max = df_top1000['AF'].max()
    # save number of alt targets
    top_1000_list.append(str(len(df_top1000)))

    if on_target_chr:
        # count number of alt sequence with same chr as on_target
        df_top1000_temp = df_top1000.loc[(
            df_top1000['Chromosome'] == on_target_chr)]
        # save count
        top_1000_list.append(str(len(df_top1000_temp)))
    else:
        top_1000_list.append('0')

    # count alt target with cCREs
    df_top1000_temp = df_top1000.loc[(
        df_top1000['Annotation_ENCODE']).notnull()]
    top_1000_list.append(str(len(df_top1000_temp)))

    # count alt target with CDS
    df_top1000_temp = df_top1000.loc[(
        df_top1000["Annotation_GENCODE"].str.contains("CDS"))]
    top_1000_list.append(str(len(df_top1000_temp)))

    # append MAF mean,std,range
    top_1000_list.append(str(MAF_mean))
    top_1000_list.append(str(MAF_std))
    top_1000_list.append(str(MAF_min)+'-'+str(MAF_max))

    # PLOT GENERATION
    crisprme_plot_MMvBUL(df_single_guide.head(
        100), guide+'_top100', out_dir, max_mm_bul_value, pam_first_nucleotide, pam_len)

    # extract targets with CDS in gencode annotation
    dff = df_single_guide.loc[df_single_guide["Annotation_GENCODE"].str.contains(
        "CDS", na=False)]
    crisprme_plot_MMvBUL(dff.head(100), guide+'_CDS_MMvBUL',
                         out_dir, max_mm_bul_value, pam_first_nucleotide, pam_len)

    # extract targets with not null encode annotation
    dff = df_single_guide.loc[df_single_guide["Annotation_ENCODE"].notnull()]
    crisprme_plot_MMvBUL(dff.head(100), guide +
                         '_encode_annotated_MMvBUL', out_dir, max_mm_bul_value, pam_first_nucleotide, pam_len)

    # extract targets with PAM creation
    dff = df_single_guide.loc[df_single_guide["PAM_creation_(fewest_mm+b)"].notnull(
    )]
    crisprme_plot_MMvBUL(dff.head(100), guide +
                         '_pam_creation_MMvBUL', out_dir, max_mm_bul_value, pam_first_nucleotide, pam_len)


in_targets_raw = sys.argv[1]
in_targets_raw_open = open(in_targets_raw, 'r')
in_humanTSGs = open(sys.argv[2], 'r')
outdir = sys.argv[3]

out_top10_counters = open(
    outdir+str(outdir.split('/')[1])+'_top10_counter.tsv', 'w')
out_top100_counters = open(
    outdir+str(outdir.split('/')[1])+'_top100_counter.tsv', 'w')
out_top1000_counters = open(
    outdir+str(outdir.split('/')[1])+'_top1000_counter.tsv', 'w')
out_top10_counters.write(
    'sgRNA\tALT-targets(fewest_mm+b)\tALT-same_chr(fewest_mm+b)\tALT-cCREs(fewest_mm+b)\tALT-CDS(fewest_mm+b)\tMAF_mean\tMAF_std\tMAF_range')
out_top10_counters.write(
    '\tALT-targets(highest_CFD)\tALT-same_chr(highest_CFD)\tALT-cCREs(highest_CFD)\tALT-CDS(highest_CFD)\tMAF_mean\tMAF_std\tMAF_range\n')
out_top100_counters.write(
    'sgRNA\tALT-targets(fewest_mm+b)\tALT-same_chr(fewest_mm+b)\tALT-cCREs(fewest_mm+b)\tALT-CDS(fewest_mm+b)\tMAF_mean\tMAF_std\tMAF_range')
out_top100_counters.write(
    '\tALT-targets(highest_CFD)\tALT-same_chr(highest_CFD)\tALT-cCREs(highest_CFD)\tALT-CDS(highest_CFD)\tMAF_mean\tMAF_std\tMAF_range\n')
out_top1000_counters.write(
    'sgRNA\tALT-targets(fewest_mm+b)\tALT-same_chr(fewest_mm+b)\tALT-cCREs(fewest_mm+b)\tALT-CDS(fewest_mm+b)\tMAF_mean\tMAF_std\tMAF_range')
out_top1000_counters.write(
    '\tALT-targets(highest_CFD)\tALT-same_chr(highest_CFD)\tALT-cCREs(highest_CFD)\tALT-CDS(highest_CFD)\tMAF_mean\tMAF_std\tMAF_range\n')

if os.path.isfile(in_targets_raw+'_annotated.txt') == False:
    header = in_targets_raw_open.readline().strip() + "\tGene_description\tGene_type"
    print('enriching with TSGs')
    geneDict = dict()
    for line in in_humanTSGs:
        split = line.strip().split('\t')
        # inserisco description e function in dict con gene come key
        geneDict[split[1]] = [split[6], split[7]]

    out_targets = open(in_targets_raw+'_annotated.txt', 'w')
    out_targets.write(header+'\n')
    for line in in_targets_raw_open:
        # annote gene if in TSGs
        split = line.strip().split('\t')
        new_cols = ['NA', 'NA']
        try:
            # check if gene name is in human TSG annotation
            # print(split[80])
            new_cols = geneDict[split[80]]
        except:
            pass
        split.extend(new_cols)
        out_targets.write('\t'.join(split)+'\n')
    out_targets.close()

print('creating dataframe')
df = pd.read_csv(in_targets_raw+'_annotated.txt', sep="\t", index_col=False,
                 na_values=['n'])
max_mm_bul_value = int(
    df["Mismatches+bulges_(fewest_mm+b)"].max())

for guide in df["Spacer+PAM"].unique():
    # identify position of PAM first nucleotide and PAM length
    pam_first_nucleotide = guide.find('N')
    pam_len = guide.count('N')
    # list to generate dataframe with counted data
    top_10_list = list()
    top_100_list = list()
    top_1000_list = list()
    top_10_list.append(str(guide))
    top_100_list.append(str(guide))
    top_1000_list.append(str(guide))

    # print('doing extraction and plots using mm+bul for guide', guide)
    # extraction_with_total(guide, df, outdir, max_mm_bul_value,
    #                       pam_first_nucleotide, pam_len, top_10_list, top_100_list, top_1000_list)

    print('doing extraction and plots using CFD for guide', guide)
    extraction_with_CFD(guide, df, outdir, top_10_list,
                        top_100_list, top_1000_list)

    # write top10 and top 100 data in file
    out_top10_counters.write('\t'.join(top_10_list)+'\n')
    out_top100_counters.write('\t'.join(top_100_list)+'\n')
    out_top1000_counters.write('\t'.join(top_1000_list)+'\n')

out_top10_counters.close()
out_top100_counters.close()
out_top1000_counters.close()
