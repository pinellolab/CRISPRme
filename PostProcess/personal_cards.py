#!/usr/bin/env python

# Libraries
from matplotlib.patches import FancyBboxPatch
import matplotlib.transforms as mtransforms
from operator import itemgetter
from os.path import isfile, join
from os import listdir
import os
import warnings
import glob
from itertools import islice
import sys
import numpy as np
import scipy.spatial.distance as sp
from math import pi
import pandas as pd
import math
import matplotlib
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.font_manager import FontProperties


# matplotlib.use("TkAgg")
matplotlib.use("Agg")

warnings.filterwarnings("ignore")

# argv 1 is Guide
# argv 2 is total value
# argv 3 is annotation file (new version)
# argv 4 is extended profile (can be 'no' for web server)
# argv 5 is second summary for barplot (optional, 'no' if not given in input)
# argv 6 is gecko Summary file
# argv 7 is for web server: create png instead of pdf (-ws) (optional)
# argv 8 is for sample/pop/superpop report (-sample HG001/EUR/TSI) with this option, the annotation file is .sample_annotation.GUIDE.sample.txt
# TODO if argv8 is selected, from argv6 get the directory containing the summary of the specific sample for gecko -> line 618
# NOTE al momento lo script funziona con gecko.annotation.summary; la comparison con i sample la fa su annotation.summary. Quando l'analisi per sample sarà fatta, la comparison
# con i sampla dovrà essere fatta sul file specifico
plt.style.use("seaborn-poster")
# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams["figure.figsize"] = [100, 80]
# SIZE_GECKO = 123411  # NOTE modify if new gecko annotations are done
# SIZE_GECKO = 111671
random.seed(a=None, version=2)
# final file containing all the results after post processing
personal_report = open(sys.argv[1], "r")
output_dir = sys.argv[2]

web_server = False
population = False
file_extension = "png"

if "-ws" in sys.argv[:]:
    web_server = True
if web_server:
    file_extension = "png"

personal_data = personal_report.readline().strip().split("\t")
personal_targets = personal_data[0]
PAM_creation = personal_data[1]
private_targets = personal_data[2]
# print(personal_targets, PAM_creation, private_targets)
list_targets = []
# private_targets = str(1)
if int(private_targets) > 0:
    for count, line in enumerate(personal_report):
        split = line.strip().split("\t")
        save = [split[2], split[20][0:5], split[21][0:5], split[22][0:5], split[14]]
        list_targets.append(save)
    # dataframe = pd.DataFrame(list_targets)
# print(provalist)
# dataframenew = dataframe[0:4]

# Bbox object around which the fancy box will be drawn.
# bb = mtransforms.Bbox([[0.3, 0.4], [0.7, 0.6]])
bb = mtransforms.Bbox([[0.2, 0.5], [0.8, 0.6]])
fig, ax = plt.subplots()
# x = 30*np.random.randn(10000)
# mu = x.mean()
# median = np.median(x)
# sigma = x.std()
# textstr = '\n'.join((
#     r'$\mu=%.2f$' % (mu, ),
#     r'$\mathrm{median}=%.2f$' % (median, ),
#     r'$\sigma=%.2f$' % (sigma, )))

# ax.hist(x, 50)
# these are matplotlib.patch.Patch properties
props_red = dict(boxstyle="round", facecolor="red", alpha=0.6)
props_blue = dict(boxstyle="round", alpha=0.6)
props_green = dict(boxstyle="round", facecolor="green", alpha=0.6)
props_sample = dict(boxstyle="round", facecolor="yellow")

# img = mpimg.imread('Immagine2.png')

# image = Image.open('Immagine2.png')
# image.resize((10, 10))
# imgplot = plt.imshow(image)
# # plt.imshow('')

# place a text box in upper left in axes coords
ax.text(
    0.1,
    0.7,
    "Personal targets: " + personal_targets,
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="top",
    bbox=props_green,
)
ax.text(
    0.43,
    0.7,
    "PAM creation: " + PAM_creation,
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="top",
    bbox=props_blue,
)
ax.text(
    0.75,
    0.7,
    "Private targets: " + private_targets,
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="top",
    bbox=props_red,
)
ax.text(
    0.4,
    0.78,
    "Sample name: " + sys.argv[1].split(".")[1],
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="top",
    bbox=props_sample,
)
ax.text(0.410, 0.57, "List of private targets", fontsize=12)
if len(list_targets) > 0:
    table = ax.table(
        cellText=list_targets,
        cellLoc="center",
        colLabels=[
            "Sequence",
            "CFD Individual",
            "CFD Reference",
            "Delta Risk Score",
            "Annotation",
        ],
        bbox=[0.05, 0.4, 0.9, 0.15],
    )
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(
                fontproperties=FontProperties(weight="bold"), fontsize=20
            )
else:
    ax.text(0.35, 0.5, "No private targets found for this sample", fontsize=11)

p_fancy = FancyBboxPatch(
    (bb.xmin, bb.ymin),
    abs(bb.width),
    abs(bb.height),
    boxstyle="round,pad=0.2",
    fc="none",
    ec="black",
    zorder=4,
)

ax.add_patch(p_fancy)
ax.set_frame_on(False)
# plt.show()
plt.axis("off")
# plt.tight_layout()
# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
#            hspace = 0, wspace = 0)
# plt.margins(0.01,0.01)

plt.savefig(
    output_dir + "/" + sys.argv[1].split(".")[1] + "_personal_card." + file_extension,
    format=file_extension,
    bbox_inches="tight",
    pad_inches=0,
)
