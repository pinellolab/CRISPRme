#!/usr/bin/env python

import glob
import sys
import sqlite3
import pandas as pd
import os

GUIDE_COLUMN = 'Spacer+PAM'
CHR_COLUMN = 'Chromosome'
POS_COLUMN = 'Start_coordinate_(highest_CFD)'
MM_COLUMN = 'Mismatches_(highest_CFD)'
BLG_COLUMN = 'Bulges_(highest_CFD)'
TOTAL_COLUMN = 'Mismatches+bulges_(highest_CFD)'
BLG_T_COLUMN = 'Bulge_type_(highest_CFD)'
CFD_COLUMN = 'CFD_score_(highest_CFD)'
RISK_COLUMN = 'CFD_risk_score_(highest_CFD)'
SAMPLES_COLUMN = 'Variant_samples_(highest_CFD)'


current_working_directory = sys.argv[1]
#print('current', current_working_directory)
job_id = current_working_directory.split('/')[-1]
#print('jobid', job_id)
guide = str(sys.argv[2])
sample = str(sys.argv[3])
integrated_personal = current_working_directory + '/' + job_id + \
    '.' + sample + '.' + guide + '.personal_targets.txt'
#print('personal', integrated_personal)
integrated_private = current_working_directory + '/' + job_id + '.' + \
    sample + '.' + guide + '.private_targets.tsv'
#print('private', integrated_private)

script_path = sys.argv[4]


# print('entrato')
path_db = glob.glob(
    current_working_directory + '/' + '*.db')[0]
#print('leggo il nome db')
path_db = str(path_db)
conn = sqlite3.connect(path_db)
c = conn.cursor()
#print('connesso to db')
# extract personal targets
result_personal = pd.read_sql_query(
    "SELECT * FROM final_table WHERE \"{}\"=\'{}\' AND \"{}\" LIKE \'%{}%\'".format(GUIDE_COLUMN, guide, SAMPLES_COLUMN, sample), conn)
result_personal = result_personal.sort_values(
    [CFD_COLUMN, TOTAL_COLUMN], ascending=[False, True])
# filter private targets
result_private = result_personal[result_personal[SAMPLES_COLUMN] == sample]
conn.commit()
conn.close()
#print('fatto queries')
# save to file
result_personal.to_csv(integrated_personal, sep='\t', index=False)
result_private.to_csv(integrated_private, sep='\t', index=False)
#print('fatto files')
# generated plots
os.system(
    f"python {script_path}/CRISPRme_plots_personal.py {integrated_personal} {current_working_directory}/imgs/ {guide}.{sample}.personal > /dev/null 2>&1")
os.system(
    f"python {script_path}/CRISPRme_plots_personal.py {integrated_private} {current_working_directory}/imgs/ {guide}.{sample}.private > /dev/null 2>&1")
os.system(
    f"rm -f {integrated_personal}")
#print('fatto plots')
