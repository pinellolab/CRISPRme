import argparse



def read_pam(pam_file:str)->list:
    """_summary_

    Args:
        pam_file (str): _description_

    Returns:
        list: _description_
    """    
    file=open(pam_file,'r')
    
    pam=file.readline().strip().split(" ")
    pam_seq=pam[0]
    pam_len=pam[1]
    if int(pam_len)>0:
        pam_in_start=False
    else:
        pam_in_start=True
    guide_len=len(pam_seq)-abs(int(pam_len))
    file.close()
    
    return [pam_seq,abs(int(pam_len)),guide_len,pam_in_start]

def read_sequence(genome_file:str)->str:
    """_summary_

    Args:
        genome_file (str): _description_

    Returns:
        str: _description_
    """    
    file=open(genome_file,'r')
    
    genome=file.readlines()
    genome=[x.strip().upper() for x in genome]
    genome="".join(genome[1:])
    
    file.close()
    
    return genome

def read_coordinate(coordinate_file:str)->list:
    """_summary_

    Args:
        coordinate_file (str): _description_

    Returns:
        list: _description_
    """    
    file=open(coordinate_file,'r')
    coordinate_list=file.readlines()
    final_list=list()
    for coord in coordinate_list:
        coord=coord.strip().split(" ")
        final_list.append([coord[0],int(coord[1]),int(coord[2])])      
    file.close()
    
    return final_list
