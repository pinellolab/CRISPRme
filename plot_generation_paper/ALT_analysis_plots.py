#!/usr/bin/env python

import sys
import time
import random
from tokenize import Floatnumber
from unicodedata import decimal
import pandas as pd
from pandas.core.indexes.api import all_indexes_same
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# from upsetplot import generate_counts
from upsetplot import UpSet
# import matplotlib as mpl
from upsetplot import from_memberships
import warnings
import matplotlib
import math
# SUPPRESS ALL WARNINGS
warnings.filterwarnings("ignore")
# do not use X11
matplotlib.use('Agg')
# set matplotlib for pdf editing
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.style.use('seaborn-poster')


def annotation_analysis(row, on_target_dict):
    # empty category list
    categories_list = list()

    if 'nan' not in str(row['PAM_creation_(highest_CFD)']):
        categories_list.append('PAM creation')
    if 'CDS' in str(row['Annotation_GENCODE']):
        categories_list.append('CDS')
        if 'nan' not in str(row['Gene_description']):
            categories_list.append('TSG')
    if 'nan' not in str(row['Annotation_ENCODE']):
        categories_list.append('cCRE')
    if str(row['Chromosome']) == on_target_dict[str(row['Spacer+PAM'])]:
        categories_list.append('Same chr')
    if len(categories_list):  # if any annotation is found return it, else empty
        return (','.join(categories_list))
    else:
        return 'empty'


def keep_one_decimal(float_number):
    int_part = str(float_number).strip().split('.')[0]
    decimal_part = str(float_number).strip().split('.')[1][0]
    float_return = float(int_part+'.'+decimal_part)
    return float_return


def num_of_decimal_zeros(float_number):
    if float_number == 0:
        return math.pow(10, -5)  # 0.00001
    if float_number >= 0.1:
        return 1  # messo in categoria 0.1-1
    decimals = str(float_number).split('.')[1]
    count_zeros = 0
    for decimal in decimals:
        if decimal == '0':
            count_zeros += 1
        else:
            break
    # add 1 to respect the exp representation (10^-1 will have count_zeros=0, but should be 1)
    return math.pow(10, -(count_zeros+1))


def generate_distribution_plot_MMBUL(original_df):
    filtered_df = original_df.loc[(
        original_df["Mismatches+bulges_(fewest_mm+b)"] <= 4)]
    filtered_df["Variant_MAF_(fewest_mm+b)"] = filtered_df["Variant_MAF_(fewest_mm+b)"].fillna(-1)

    # If multiple AFs (haplotype with multiple SNPs), take min AF
    # Approximation until we have haplotype frequencies
    filtered_df["AF"] = filtered_df["Variant_MAF_(fewest_mm+b)"].astype(
        str).str.split(',')
    filtered_df["AF"] = filtered_df["AF"].apply(lambda x: min(x))
    filtered_df["AF"] = pd.to_numeric(filtered_df["AF"])
    filtered_df.sort_values(
        ['Mismatches+bulges_(fewest_mm+b)'], inplace=True, ascending=True)

    andamento_ALT_MAF005 = list()
    andamento_ALT_MAF05 = list()
    andamento_ALT_MAF0 = list()
    altTarget_MAF005 = 0
    altTarget_MAF05 = 0
    altTarget_MAF0 = 0

    for index, row in filtered_df.iterrows():
        if row['AF'] > 0.005:
            altTarget_MAF005 += 1
        if row['AF'] > 0.05:
            altTarget_MAF05 += 1
        if row['AF'] >= 0:
            altTarget_MAF0 += 1
        andamento_ALT_MAF005.append(altTarget_MAF005)
        andamento_ALT_MAF05.append(altTarget_MAF05)
        andamento_ALT_MAF0.append(altTarget_MAF0)

    plt.plot(andamento_ALT_MAF0, label='MAF>0')
    plt.plot(andamento_ALT_MAF005, label='MAF>0.005')
    plt.plot(andamento_ALT_MAF05, label='MAF>0.05')

    plt.ylabel('ALT Targets')
    plt.xlabel('Targets')
    plt.title(
        'Distribution of targets with different MAFs filtered with MM+BUL <= 4')
    plt.legend()

    plt.tight_layout()
    plt.savefig(out_folder+'distribution_plt_MMBUL.pdf')
    plt.clf()
    plt.close('all')


def generate_upset_plot_MMBUL(original_df):
    # MMBUL analysis
    df_alt = original_df.loc[(original_df['REF/ALT_origin_(fewest_mm+b)']
                              == 'alt') & (original_df['Mismatches+bulges_(fewest_mm+b)'] <= 4)]

    # create dict
    on_target_dict = dict()
    for guide in df_alt["Spacer+PAM"].unique():
        on_target_dict[str(guide)] = 'empty'

    # extract on-target chr
    try:
        on_target_chr = original_df.loc[(
            original_df['Mismatches+bulges_(fewest_mm+b)'] == 0)]
        on_target_guide = on_target_chr.iloc[0]['Spacer+PAM']
        on_target_chr = on_target_chr.iloc[0]['Chromosome']
        on_target_dict[str(on_target_guide)] = str(on_target_chr)
    except:
        print('on target not found for guide')
    # create empty categories col
    df_alt['Categories'] = 'empty'

    # process df to obtain categories of belonging for each target
    for index in df_alt.index:
        categories_list = list()
        if 'CDS' in str(df_alt.loc[index, 'Annotation_GENCODE']):
            categories_list.append('CDS')
        if 'nan' not in str(df_alt.loc[index, 'Annotation_ENCODE']):
            categories_list.append('ENCODE')
        if 'nan' not in str(df_alt.loc[index, 'Gene_description']) and 'CDS' in str(df_alt.loc[index, 'Annotation_GENCODE']):
            categories_list.append('TSG')
        if str(df_alt.loc[index, 'Chromosome']) == on_target_dict[str(df_alt.loc[index, 'Spacer+PAM'])]:
            categories_list.append('On-Target_Chromosome')
        if len(categories_list):
            df_alt.loc[index, 'Categories'] = ','.join(categories_list)

    # remove targets with empty categories
    df_alt = df_alt.loc[(df_alt['Categories'] != 'empty')]
    # collect categories per target
    categories_per_target = from_memberships(
        df_alt.Categories.str.split(','), data=df_alt)
    # print(categories_per_target)
    # create figure
    figu = plt.figure()
    upset_plot = UpSet(categories_per_target, show_counts=True,
                       sort_by='cardinality', sort_categories_by=None)
    upset_plot.plot(fig=figu)
    plt.title('ALT targets overlapping categories filtered with MM+BUL <= 4')
    # plt.tight_layout()
    plt.savefig(
        out_folder+'overlapping_alt_targets_categories_MMBUL.pdf')
    plt.clf()
    plt.close('all')


def generate_heatmap_CFD(original_df):
    # discard on-targets and quasi-on-targets(1 mm+bul)
    original_df = original_df.loc[(
        original_df["Mismatches+bulges_(highest_CFD)"] > 1)]
    # filter df to use only heatmap related columns
    df_heatmap = original_df[[
        'CFD_score_(highest_CFD)', 'Variant_MAF_(highest_CFD)']]
    df_heatmap = df_heatmap.loc[(
        df_heatmap["Variant_MAF_(highest_CFD)"].notnull())]

    # MAF conversion and filtering
    # df_heatmap["Variant_MAF_(highest_CFD)"] = df_heatmap["Variant_MAF_(highest_CFD)"].fillna(-1)
    df_heatmap["Variant_MAF_(highest_CFD)"] = df_heatmap["Variant_MAF_(highest_CFD)"].astype(
        str).str.split(',')
    df_heatmap["Variant_MAF_(highest_CFD)"] = df_heatmap["Variant_MAF_(highest_CFD)"].apply(
        lambda x: min(x))
    df_heatmap["Variant_MAF_(highest_CFD)"] = pd.to_numeric(
        df_heatmap["Variant_MAF_(highest_CFD)"], downcast="float")
    # df_heatmap = df_heatmap.loc[(df_heatmap['Variant_MAF_(highest_CFD)']) >= 0]
    # conversion to count of decimal zeros
    df_heatmap["Variant_MAF_(highest_CFD)"] = df_heatmap["Variant_MAF_(highest_CFD)"].apply(
        lambda x: num_of_decimal_zeros(x))

    # df_fake_target = pd.DataFrame([[0, 0.00001], [0, 0.0001], [0, 0.001], [0, 0.01], [
    #                               0.9, 0.1]], columns=['CFD_score_(highest_CFD)', 'Variant_MAF_(highest_CFD)'])
    # df_heatmap = df_heatmap.append(df_fake_target)

    # df_fake_target = pd.DataFrame([[1, 0.00001], [1, 0.0001], [1, 0.001], [1, 0.01], [
    #                               1, 0.1]], columns=['CFD_score_(highest_CFD)', 'Variant_MAF_(highest_CFD)'])
    # df_heatmap = df_heatmap.append(df_fake_target)

    # CFD score rounding to 1 decimal
    df_heatmap['CFD_score_(highest_CFD)'] = df_heatmap['CFD_score_(highest_CFD)'].astype(
        float)
    # df_heatmap['CFD_score_(highest_CFD)'] = df_heatmap['CFD_score_(highest_CFD)'].apply(
    #     lambda x: round(x, 1))
    df_heatmap['CFD_score_(highest_CFD)'] = df_heatmap['CFD_score_(highest_CFD)'].apply(
        lambda x: keep_one_decimal(x))

    print(df_heatmap)

    df_table = df_heatmap.groupby(
        ["Variant_MAF_(highest_CFD)", "CFD_score_(highest_CFD)"]).size().reset_index(name="Value")

    # df2 = pd.DataFrame([[0.1]*len(list(df_table['Variant_MAF_(highest_CFD)']))],
    #                    columns=list(df_table['Variant_MAF_(highest_CFD)']), index=[1])
    # df_table.append(df2)

    table = df_table.pivot('CFD_score_(highest_CFD)',
                           'Variant_MAF_(highest_CFD)', 'Value').fillna(0.1)

    print(table)

    cbar_ticks = [10**0, 10**1, 10**2, 10**3, 10**4, 10**5, 10**6]
    vmax = 10**6
    vmin = 10**0
    # formatter = tkr.ScalarFormatter(useMathText=True)
    log_norm = LogNorm(vmin=vmin, vmax=vmax)
    # formatter.set_scientific(True)

    figu = plt.figure()
    plt_heatmap = sns.heatmap(table, cmap="RdYlBu", fmt='.0f', annot=True, vmax=vmax, vmin=vmin, norm=log_norm,
                              cbar_kws={"ticks": cbar_ticks, 'label': 'Target sites'}, xticklabels=False, yticklabels=False)
    plt_heatmap.collections[0].colorbar.ax.yaxis.set_ticks([], minor=True)
    # plt_heatmap.collections[0].colorbar.ax.tick_params(labelsize=13)

    # set labels and position of ticks
    plt_heatmap.set_xticks(np.arange(0, 6, step=1))
    plt_heatmap.set_yticks(np.arange(0, 11, step=1))
    plt_heatmap.set_yticklabels(
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    plt_heatmap.set_xticklabels(
        labels=[0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
    plt_heatmap.invert_yaxis()

    plt.xlabel("Variant MAF")
    plt.ylabel("CFD score")
    # plt_heatmap.invert_xaxis()
    plt.tight_layout()
    plt.savefig(out_folder+'heatmap_CFD.pdf', transparent=True)
    plt.clf()
    plt.close('all')


def generate_distribution_plot_CFD(original_df, name):
    filtered_df = original_df

    plt.figure()
    for guide in filtered_df['Spacer+PAM'].unique():
        print('analyzing guide', guide)
        guide_df = filtered_df.loc[(filtered_df['Spacer+PAM'] == guide)]
        andamenti = list()
        andamento_ALT_MAF0 = list()
        altTarget_MAF0 = 0

        guide_df.sort_values(
            ['CFD_score_(highest_CFD)'], inplace=True, ascending=False)
        alt_list = guide_df['REF/ALT_origin_(highest_CFD)'].tolist()
        # print(alt_list[:10])

        # read af to select alt targets
        for target in alt_list:
            if target == 'alt':
                altTarget_MAF0 += 1
            andamento_ALT_MAF0.append(altTarget_MAF0)

        # andamenti con distribuzione andamento di ogni guida
        andamenti.append(andamento_ALT_MAF0)
        # plt line for single guide
        plt.plot(andamento_ALT_MAF0)

    andamentiArray = np.array(andamenti)
    media = np.mean(andamentiArray, axis=0)
    standarddev = np.std(andamentiArray, axis=0)
    standarderr = standarddev/np.sqrt(np.amax(andamentiArray))
    z_score = 1.96  # for confidence 95%
    lowerbound = media-(z_score*standarderr)
    upperbound = media+(z_score*standarderr)

    plt.plot(media)
    # plt.fill_between(range(len(media)), lowerbound,
    #                  upperbound, alpha=0.10)

    # plt.ylabel('ALT Targets')
    if 'log' in name:
        plt.yscale('log')
        plt.xscale('log')

    plt.xlabel('Targets')
    plt.ylabel('ALT Targets')
    plt.title('Distribution of targets with different MAFs filtered with MAF>0')
    list_labels = list(filtered_df['Spacer+PAM'].unique())
    list_labels.append('Mean distribution')

    plt.legend([])
    # plt.xticks(fontsize=13)
    # plt.yticks(fontsize=13)

    plt.tight_layout()
    plt.savefig(out_folder+name+'_distribution_plt_CFD.pdf',
                transparent=True)
    plt.clf()
    plt.close('all')
    print('done')


def generate_upset_plot_CFD(original_df):
    # CFD analysis
    # df_alt = original_df.loc[(original_df['REF/ALT_origin_(highest_CFD)']== 'alt') & (original_df["CFD_score_(highest_CFD)"] >= 0.1)]
    df_alt = original_df.loc[(
        original_df['REF/ALT_origin_(highest_CFD)'] == 'alt')]

    # create dict
    on_target_dict = dict()
    for guide in df_alt["Spacer+PAM"].unique():
        on_target_dict[str(guide)] = 'empty'

    # extract on-target chr
    try:
        on_target_chr = original_df.loc[(
            original_df['Mismatches+bulges_(highest_CFD)'] == 0)]
        on_target_guide = on_target_chr.iloc[0]['Spacer+PAM']
        on_target_chr = on_target_chr.iloc[0]['Chromosome']
        on_target_dict[str(on_target_guide)] = str(on_target_chr)
    except:
        print('on target not found for guide')

    # compute analysis to compute annotations
    df_alt['Categories'] = df_alt.apply(
        lambda row: annotation_analysis(row, on_target_dict), axis=1)

    # remove targets with empty categories
    df_alt = df_alt.loc[(df_alt['Categories'] != 'empty')]

    # collect categories per target
    categories_per_target = from_memberships(
        df_alt.Categories.str.split(','), data=df_alt)

    # create figure
    figu = plt.figure()
    upset_plot = UpSet(categories_per_target, show_counts=True,
                       sort_by='cardinality', sort_categories_by=None)
    upset_plot.plot(fig=figu)
    plt.title('ALT targets overlapping categories')
    # plt.tight_layout()
    plt.savefig(
        out_folder+'overlapping_alt_targets_categories_CFD.pdf', transparent=True)
    plt.clf()
    plt.close('all')


inTargets = sys.argv[1]  # read targets
out_folder = sys.argv[2]  # folder for output images

print('starting generating distribution and upset plots')
# create dataframe with file
start = time.time()
original_df_read = pd.read_csv(inTargets, sep="\t", index_col=False,
                               na_values=['n'])
print('read in time', time.time()-start)
# call to plot generation CFD with original data
# generate_distribution_plot_CFD(original_df_read, 'no_filter')
# generate_distribution_plot_CFD(original_df_read, 'no_filter_log')
# generate_upset_plot_CFD(original_df_read)
start = time.time()
generate_heatmap_CFD(original_df_read)
print('heatmap generated in time', time.time()-start)
# generate_upset_log_barplot_CFD()
# call to plot generation MM_BUL
# generate_distribution_plot_MMBUL(original_df)
# generate_upset_plot_MMBUL(original_df)

# call to plot generation CFD with CFD>=0.1
# print('execute with cfd>=0.1')
# cfd_01_df = original_df_read.loc[(
#     original_df_read["CFD_score_(highest_CFD)"] >= 0.1)]
# generate_distribution_plot_CFD(cfd_01_df, 'cfd>=0.1')
