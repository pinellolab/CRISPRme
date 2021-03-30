#Convert pam and get guides from sequence
import sys
import os
import itertools
from Bio.Seq import Seq
import re
def getGuides(extracted_seq, pam, len_guide, pam_begin):
    
    len_pam = len(pam)
    #dict
    len_guide = int(len_guide)
    pam_dict = {
        'A':  "ARWMDHV",
        'C':  "CYSMBHV",
        'G':  "GRSKBDV",
        'T':  "TYWKBDH",
        'R':  "ARWMDHVSKBG",
        'Y':  "CYSMBHVWKDT",
        'S':  "CYSMBHVKDRG",
        'W':  "ARWMDHVYKBT",
        'K':  "GRSKBDVYWHT",
        'M':  "ARWMDHVYSBC",
        'B':  "CYSMBHVRKDGWT",
        'D':  "ARWMDHVSKBGYT",
        'H':  "ARWMDHVYSBCKT",
        'V':  "ARWMDHVYSBCKG",
        'N':  "ACGTRYSWKMBDHV",
    }
    list_prod = []
    for char in pam:
        list_prod.append(pam_dict[char])

    iupac_pam = []          #NNNNNNN NGG
    for element in itertools.product(*list_prod):
        iupac_pam.append(''.join(element))

    rev_pam = str(Seq(pam).reverse_complement())
    list_prod = []
    for char in rev_pam:
        list_prod.append(pam_dict[char])

    iupac_pam_reverse = []        #CCN NNNNNNN  -> results found with this pam must be reverse complemented
    for element in itertools.product(*list_prod):
        iupac_pam_reverse.append(''.join(element))

    
    extracted_seq = extracted_seq.upper()
    len_sequence = len(extracted_seq)    
    guides = []
    for pam in iupac_pam:
        pos = ([m.start() for m in re.finditer('(?=' + pam + ')', extracted_seq)])
        if pos:
            for i in pos:
                if pam_begin:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    guides.append(extracted_seq[i + len_pam : i +len_pam + len_guide])
                else:
                    if i < len_guide:
                        continue
                    #guides.append(extracted_seq[i-len_guide:i+len_pam])           # i is position where first char of pam is found, eg the N char in NNNNNN NGG
                    #print('1 for:' , extracted_seq[i-len_guide:i])
                    guides.append(extracted_seq[i-len_guide:i])
    for pam in iupac_pam_reverse:       #Negative strand
        pos = ([m.start() for m in re.finditer('(?=' + pam + ')', extracted_seq)])
        if pos:
            for i in pos:
                if pam_begin:
                    if i < len_guide:
                        continue
                    guides.append(str(Seq(extracted_seq[i-len_guide:i]).reverse_complement()))
                else:
                    if i > (len_sequence - len_guide - len_pam):
                        continue
                    #guides.append(str(Seq(extracted_seq[i:i+len_pam+len_guide]).reverse_complement()))         # i is position where first char of pam is found, eg the first C char in CCN NNNNNN
                    #print('2 for:', str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
                    guides.append(str(Seq(extracted_seq[i + len_pam : i + len_guide + len_pam]).reverse_complement()))
    return guides
    #return guides for when adding to app.py