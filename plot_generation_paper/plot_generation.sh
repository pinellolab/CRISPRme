#!/bin/bash

#SBATCH --job-name=revision_plot
#SBATCH --priority=TOP
#SBATCH --cpus-per-task=1
#SBATCH --mem=150000
#SBATCH --time=12:00:00
#SBATCH --partition=medium

#start creation of plots
# python gene_annotation_extraction_plots.py ../NGG_analysis/TCACTATGCTGCCGCCCAGTNNN.targets.tsv Human_TSGs.txt ../NGG_analysis/
# python ALT_analysis_plots.py ../NGG_analysis/TCACTATGCTGCCGCCCAGTNNN.targets.tsv_annotated.txt ../NGG_analysis/

#analyze NGG-complete results
# python gene_annotation_extraction_plots.py ../NGG_analysis/NGG-complete.integrated.tsv Human_TSGs.txt ../NGG_analysis/
python ALT_analysis_plots.py ../NGG_analysis/NGG-complete.integrated.tsv_annotated.txt ../NGG_analysis/
