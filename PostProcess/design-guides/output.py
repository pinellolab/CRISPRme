"""
"""

from typing import Tuple


def write_results(guide_list:list,guide_file:str,positions:Tuple,chr_name:str)->None:
    """_summary_

    Args:
        guide_list (list): _description_
        guide_file (str): _description_
        positions (Tuple): _description_
        chr_name (str): _description_
    """       
    guide_file_out=open(guide_file,'w')
    guide_summary_out=open(guide_file.replace(".txt","")+".summary.tsv",'w')
    guide_summary_out.write("sgRNA\tPAM\tsgRNA_length\tPAM_length\tchr\tstrand\n")
        
    for i,guide in enumerate(guide_list):
        guide_file_out.write(guide[0]+"\n")
        guide_summary_out.write(guide[0]+"\t"+guide[1]+"\t"+str(len(guide[0])-len(guide[1]))+"\t"+str(len(guide[1]))+"\t"+chr_name+"\t"+guide[2]+"\n")
    guide_file_out.close()
    guide_summary_out.close()
