#!/usr/bin/env python

import os 
import sys
from datetime import datetime


vcf_dir = sys.argv[1]
files = os.listdir(vcf_dir)
vcf_name = sys.argv[2]
output_folder = sys.argv[3]
log = sys.argv[4]

for f in files:
    if "vcf.gz" in f:
        print("Creating dictionary for file", f)
        number = f.split(".")[1].replace('chr', '')
        os.system("./creazione_dizionari_zipped.py "+vcf_dir+"/"+f+" "+str(number))
        os.system("mv "+'my_dict_chr' + number + '.json'+" "+output_folder+"/dictionaries_"+vcf_name+"/") 
       
os.system('echo "Dictionaries\tEnd\t'+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+log)
