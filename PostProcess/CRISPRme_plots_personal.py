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

# # Read file
# df = pd.read_csv(sys.argv[1], sep="\t",
#                  index_col=False, na_values=['n'], nrows=1200)
# out_folder = sys.argv[2]

# # Make index column that numbers the OTs starting from 1
# df = df.reset_index()
# df["index"] += 1

# Read file
df = pd.read_csv(sys.argv[1], sep="\t",
                 index_col=False, na_values=['n'], nrows=100)
out_folder = sys.argv[2]
guide = sys.argv[3]

# Make index column that numbers the OTs starting from 1
df = df.reset_index()
df["index"] += 1

# If prim_AF = 'n', then it's a ref-nominated site, so we enter a fake numerical AF
# This will cause a warning of invalid sqrt later on, but that's fine to ignore
df["prim_AF"] = df["prim_AF"].fillna(-1)

# If multiple AFs (haplotype with multiple SNPs), take min AF
# Approximation until we have haplotype frequencies
df["AF"] = df["prim_AF"].astype(str).str.split(',')
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

# # Color by annotation, with coding > accessible > neither
# # Note that coding & accessible categories are not mutually exclusive, but coding is more important than accessible, so if both, it's colored as coding
# df["color_str"] = np.where(df["gene_annotation"] ==
#                            "CDS", "red", "gray")  # coding
# df["color_str"] = np.where(df["annotation"].str.contains(
#     "HSC-1", na=False), "blue", df["color_str"])  # DHS, accessible in HSCs
# df["color"] = df["color_str"].apply(lambda x: transparent_red if x == "red" else (
#     transparent_blue if x == "blue" else transparent_gray))

# """
# Linear, annotated, top 100 (for supplement)
# """
# # Plot data
# ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)",
#                      s="ref_AF", c="color", marker='*', zorder=1)
# df.plot.scatter(x="index", y="highest_CFD_score(alt)",
#                 s="plot_AF", c="color", zorder=2, ax=ax)
# # plt.title("Top CRISPRme-identified sites for sgRNA 1617")
# plt.xlabel("Candidate off-target site")
# plt.ylabel("CFD score")

# # Boundaries
# plt.xlim(xmin=0.01, xmax=100)
# plt.ylim(ymin=0, ymax=1)

# # Arrows
# for x, y, z in zip(df["index"], df["highest_CFD_score(ref)"], df["highest_CFD_score(alt)"]-df["highest_CFD_score(ref)"]):
#     plt.arrow(x, y+0.01, 0, z-0.02, color='green', head_width=0.75,
#               head_length=0.015, length_includes_head=True, alpha=0.5, zorder=0)
#     # +/- to avoid overlap of arrow w/ points

# Size legend
s1 = mlines.Line2D([], [], marker='o', label='1', linestyle='None',
                   markersize=math.sqrt(math.sqrt((1+0.001)*1000)), color='black')
s01 = mlines.Line2D([], [], marker='o', label='0.1', linestyle='None',
                    markersize=math.sqrt(math.sqrt((0.1+0.001)*1000)), color='black')
s001 = mlines.Line2D([], [], marker='o', label='0.01', linestyle='None',
                     markersize=math.sqrt(math.sqrt((0.01+0.001)*1000)), color='black')
# plt.gca().add_artist(plt.legend(
#     handles=[s1, s01, s001], title="Allele frequency", loc=9))

# # Shape & color legend
# star = mlines.Line2D([], [], marker='*', label='Reference',
#                      linestyle='None', markersize=10, color='black')
# circle = mlines.Line2D([], [], marker='o', label='Alternative',
#                        linestyle='None', markersize=10, color='black')
# red = mpatches.Patch(color=transparent_red, label='Coding')
# blue = mpatches.Patch(color=transparent_blue, label='DHS in HSC-1')
# gray = mpatches.Patch(color=transparent_gray, label='Neither')
# plt.legend(handles=[star, circle, red, blue, gray])

# # Save
# plt.tight_layout()
# plt.savefig(
#     out_folder+f"CRISPRme_top_100_linear_annotated_for_supplement_{guide}.png")
# plt.clf()


"""
Log, ref/alt, top 1000: for main text
"""
# matplotlib plot settings
plt.rcParams["figure.dpi"] = 600
plt.rcParams["figure.figsize"] = 4.5, 2.5
plt.rcParams.update({'font.size': 7})
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# seaborn plot settings - also need to comment out line 114 and uncomment line 115
# sns.set_style('white')
# sns.set_context('talk') # alternatively, 'poster'
# k = 2 # scale factor
# fig = plt.figure(figsize=(7*k, 2*k))
# ax = fig.add_subplot(1,1,1)

# Plot data
ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)",
                     s="ref_AF", c=transparent_red, zorder=1)
# ax = df.plot.scatter(x="index", y="highest_CFD_score(ref)", s="ref_AF", c=transparent_red, zorder=1, ax=ax)
df.plot.scatter(x="index", y="highest_CFD_score(alt)",
                s="plot_AF", c=transparent_blue, zorder=2, ax=ax)
ax.set_xscale("log")
# plt.title("Top CRISPRme-identified sites for sgRNA 1617")
plt.xlabel("Candidate off-target site")
plt.ylabel("CFD score")

# Boundaries
plt.xlim(xmin=0.9, xmax=100)
plt.ylim(ymin=0, ymax=1)

# Arrows
for x, y, z in zip(df["index"], df["highest_CFD_score(ref)"], df["highest_CFD_score(alt)"]-df["highest_CFD_score(ref)"]):
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
plt.savefig(out_folder+f"CRISPRme_top_1000_log_for_main_text_{guide}.png")
plt.clf()
