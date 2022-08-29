#!/bin/bash

#required input
#$1 is integrated file

#extract all targets that contains HGDP sample and have CFD score >=0.2 and risk score >=0.1
#$1 is bestMerge file
# awk '$14 ~ "HGDP" && $23>=0.1 && $21>=0.2' $1 >HGDP.cfd02.risk01.txt
#extract value with HGDP samples and CFD>=0.2 and RISK_SCORE>=0.1
awk '$27 ~ "HGDP" && $19>=0.2 && $2>=0.1' $1 >HGDP.cfd02.risk01.txt

#filter out all targets that are shared with 1000G
awk '$27 !~ "NA[0-9]" && $27 !~ "HG[0-9]"' HGDP.cfd02.risk01.txt >HGDPonly.cfd02.risk01.txt

#add superpopulation to each target, based on the samples sharing the target
#hg38_HGDP.samplesID.txt is the file with the data of each sample, POP, SUPERPOP, SEX
# python add_superpop.py hg38_HGDP.samplesID.txt HGDPonly.cfd02.risk01.txt >HGDPonly.cfd02.risk01.superpop.txt

#filter targets with one superpopulation
# awk '$49 !~ ","' HGDPonly.cfd02.risk01.superpop.txt >HGDPonly.cfd02.risk01.only_one_superpop.txt

#filter targets with one sample (unique)
# awk '$14 !~ ","' HGDPonly.cfd02.risk01.superpop.txt >HGDPonly.cfd02.risk01.only_one_sample.txt

# #exctract header from integrated
# #$2 is integrated file
# head -1 $2 >intergrated.header.txt
# #extract targets with HGDP01211 and risk >0
# awk '$23 ~ "HGDP01211" && $18>0' $2 >HGDP01211.cfd0.risk0.all.txt

# #filter targets with HGDP01211 and 1000G
# awk '$23 ~ "NA[0-9]" || $23 ~ "HG[0-9]"' HGDP01211.cfd0.risk0.all.txt >HGDP01211.cfd0.risk0.1000G.txt
# cat intergrated.header.txt HGDP01211.cfd0.risk0.1000G.txt >HGDP01211.cfd0.risk0.1000G.txt.tmp
# mv HGDP01211.cfd0.risk0.1000G.txt.tmp HGDP01211.cfd0.risk0.1000G.txt

# #filter targets with HGDP01211 no 1000G and other HGDP samples
# awk '$23 !~ "NA[0-9]" && $23 !~ "HG[0-9]" && $23 != "HGDP01211"' HGDP01211.cfd0.risk0.all.txt >HGDP01211.cfd0.risk0.HGDPonly.txt
# cat intergrated.header.txt HGDP01211.cfd0.risk0.HGDPonly.txt >HGDP01211.cfd0.risk0.HGDPonly.txt.tmp
# mv HGDP01211.cfd0.risk0.HGDPonly.txt.tmp HGDP01211.cfd0.risk0.HGDPonly.txt

#filter target with HGDP01211 only
# awk '$23=="HGDP01211"' HGDP01211.cfd0.risk0.all.txt >HGDP01211.cfd0.risk0.private.txt
# cat intergrated.header.txt HGDP01211.cfd0.risk0.private.txt >HGDP01211.cfd0.risk0.private.txt.tmp
# mv HGDP01211.cfd0.risk0.private.txt.tmp HGDP01211.cfd0.risk0.private.txt

# #create figures present in figure2 main paper
# #this command creates figure 2b and 2c
# #HGDPonly.cfd02.risk01.superpop.txt contains targets from HGDP with CFD >=0.2 and risk score >=0.1
# python HGDPplots.py hg38_HGDP.samplesID.txt HGDPonly.cfd02.risk01.superpop.txt 0.1 0.2
# #HGDPonly.cfd02.risk01.only_one_sample.txt contains targets from HGDP with one sample (private) with CFD >=0.2 and risk score >=0.1
# python HGDPplots.py hg38_HGDP.samplesID.txt HGDPonly.cfd02.risk01.only_one_sample.txt 0.1 0.2

# #create figures 2d
# python CRISPRme_plots.py HGDP01211.cfd0.risk0.1000G.txt ./ 1000G
# python CRISPRme_plots.py HGDP01211.cfd0.risk0.HGDPonly.txt ./ HGDPonly
# python CRISPRme_plots.py HGDP01211.cfd0.risk0.private.txt ./ HGDP01211
