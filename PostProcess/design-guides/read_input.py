import argparse



def read_pam(pam_file:str)->list:
    file=open(pam_file,'r')
    pam=file.readline().strip().split(" ")
    pam_seq=pam[0]
    pam_len=pam[1]
    if int(pam_len)>0:
        pam_in_start=False
    else:
        pam_in_start=True
    
    file.close()
    
    return [pam_seq,abs(int(pam_len)),pam_in_start]

def read_sequence(genome_file:str)->str:
    file=open(genome_file,'r')
    genome=file.readlines()
    genome="".join(genome[1:])
    
    file.close()
    
    return genome

def read_coordinate(coordinate_file:str)->list:
    file=open(coordinate_file,'r')
    coordinate_list=file.readlines()
    final_list=list()
    for coord in coordinate_list:
        coord=coord.strip().split(" ")
        final_list.append([coord[0],int(coord[1]),int(coord[2])])      
    file.close()
    
    return final_list




parser = argparse.ArgumentParser(description='Given input fasta file, chromosome coordinates file (BED format) and PAM sequence, find all possible sgRNA guides in the fasta file (SUPPORTS IUPAC NOTATIONS)')
parser.add_argument('--fasta', type=str, help='fasta file',dest='genome_file')
parser.add_argument('--pam', type=str, help='PAM file',dest='pam_file')
parser.add_argument('--coordinate', type=str, help='bed file',dest='coordinate_file')


args=parser.parse_args()
print(args.genome_file)
print(args.pam_file)
print(args.coordinate_file)

ll=read_coordinate(args.coordinate_file)
print(ll)

ll=read_pam(args.pam_file)
print(ll)

ll=read_sequence(args.genome_file)
print(ll)
