#!/usr/bin/env python

import sys
import os

vcf_folder = sys.argv[1]
files = os.listdir(vcf_folder)
vcf_name = sys.argv[2]
output_folder = sys.argv[3]
ref_folder = sys.argv[4]
pam_file = sys.argv[5]
guide_file = sys.argv[6]
bMax = sys.argv[7]
mm = sys.argv[8]
bDNA = sys.argv[9]
bRNA = sys.argv[10]

for f in files:
    if 'vcf' in f:
        print("Getting indels for file", f)
        number = f.split(".")[1]
        os.system("./indels_process.sh "+vcf_folder+"/"+f+" "+number+" "+output_folder+" "+vcf_name+" "+ref_folder+" "+pam_file+" "+guide_file+" "+bMax+" "+mm+" "+bDNA+" "+bRNA)

