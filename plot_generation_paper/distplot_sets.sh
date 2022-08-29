#!/bin/bash

#analyse 1000G set
echo "faccio 1000G"
python individual_boxplot.py /attachedvolume/scancellieri/CRISPRme/Results/sg1617-1000G_OX2V337MO4/CTAACAGTTGCTTTTATCACNNN+NGG_hg38+hg38_1000G_6+2_CFD_integrated_results.tsv /attachedvolume/scancellieri/CRISPRme/Results/sg1617-1000G+HGDP_YD92OOXH5S/CTAACAGTTGCTTTTATCACNNN+NGG_hg38+hg38_HGDP_6+2_CFD_integrated_results.tsv /attachedvolume/scancellieri/CRISPRme/Results/sg1617-1000G_OX2V337MO4/.sampleID.txt ./ 1000G

#analyse HGDP set
echo "faccio HGDP"
python individual_boxplot.py /attachedvolume/scancellieri/CRISPRme/Results/sg1617-HGDP_AMYBX6SGHQ/CTAACAGTTGCTTTTATCACNNN+NGG_hg38+hg38_HGDP_6+2_CFD_integrated_results.tsv /attachedvolume/scancellieri/CRISPRme/Results/sg1617-1000G+HGDP_YD92OOXH5S/CTAACAGTTGCTTTTATCACNNN+NGG_hg38+hg38_HGDP_6+2_CFD_integrated_results.tsv /attachedvolume/scancellieri/CRISPRme/Results/sg1617-HGDP_AMYBX6SGHQ/.sampleID.txt ./ HGDP
