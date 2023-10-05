from typing import Tuple
from utils import reverse_complement

def extract_guides_from_genome(positions:Tuple,genome:str,guide_len:int,pam_len:int,pam_in_start:bool)->list:
    """_summary_

    Args:
        positions (Tuple): _description_
        genome (str): _description_
        guide_len (int): _description_
        pam_len (int): _description_
        pam_in_start (bool): _description_

    Returns:
        list: _description_
    """     
    guides_list=list()
    ##pam at start == true
    if pam_in_start:
        for pos in positions[0]:
            guide=genome[pos+pam_len:pos+pam_len+guide_len] # type: ignore
            if 'N' in guide:
                continue
            guide="N"*pam_len+guide
            pam=genome[pos:pos+pam_len]
            if 'N' in pam:
                continue
            guides_list.append([guide,pam,"forward"])
            
        for pos in positions[1]:
            guide=genome[pos:pos+guide_len]
            if 'N' in guide:
                continue
            guide=reverse_complement(guide)
            guide="N"*pam_len+guide
            pam=genome[pos+guide_len:pos+pam_len+guide_len]
            if 'N' in pam:
                continue
            guides_list.append([guide,pam,"reverse"])
            
        return guides_list
    else:
        ##pam at start == false
        for pos in positions[0]:
            guide=genome[pos:pos+guide_len]
            if 'N' in guide:
                continue
            guide=guide+"N"*pam_len
            pam=genome[pos+guide_len:pos+pam_len+guide_len]
            if 'N' in pam:
                continue
            guides_list.append([guide,pam,"forward"])
            
        for pos in positions[1]:
            guide=genome[pos+pam_len:pos+pam_len+guide_len]
            if 'N' in guide:
                continue
            guide=reverse_complement(guide)
            guide=guide+"N"*pam_len
            pam=genome[pos:pos+pam_len]
            if 'N' in pam:
                continue
            guides_list.append([guide,pam,"reverse"])
            
    return guides_list
