#OLD, VIEW summary_by_guide_position.py FOR NEW CHANGES

#Script che crea la tabella Summary by Guide, dato in input un file targets.txt, conta il numero di targets trovato per ogni combinazione di
#mms-bulge. Inoltre, se è presente la colonna Var_uniq, calcola il numero di Var_uniq per quella categoria di mms-bulge

#sys1 file risultati
# sys2 mms
# sys3 bulges DNA
# sys4 bulges RNA
# sys5 guida
# sys6 is jobid
# sys7 is type of post-process done ('No' -> no post process done, cannot count uniq_var | 'Uniq' -> post process done, can count uniq_var)
# inserire anche numero guide

#BUG nel conteggio  -> ricontrollando è stato risolto(?)
import sys
import numpy as np
import pandas as pd
# with open(sys.argv[1]) as targets:
#     mms = int(sys.argv[2])
#     bulges_dna = int(sys.argv[3])
#     bulges_rna = int(sys.argv[4])
#     total_count = np.zeros((mms + 1, bulges_dna +1 + bulges_rna +1 + 1))
#     print (total_count)
#     header = targets.readline()
#     for line in targets:
#         line = line.strip().split('\t')
#         total_count[int(line[6])][int(line[7])] = total_count[int(line[6])][int(line[7])] + 1
    
#     np.savetxt('summary_table.txt', total_count, delimiter='\t')
guide = sys.argv[5]
type_post = sys.argv[7]

with open(sys.argv[1]) as targets:
    mms = int(sys.argv[2])
    bulges_dna = int(sys.argv[3])
    bulges_rna = int(sys.argv[4])
    total_count_x = np.zeros((mms + 1, 1))
    total_count_dna = np.zeros((mms + 1,bulges_dna + 1))
    total_count_rna = np.zeros((mms + 1,bulges_rna + 1))
    header = targets.readline()
    if type_post == 'Uniq':
        total_count_x_uniq = np.zeros((mms + 1, 1))
        total_count_dna_uniq = np.zeros((mms + 1,bulges_dna + 1))
        total_count_rna_uniq = np.zeros((mms + 1,bulges_rna + 1))
        for line in targets:
            
            line = line.strip().split('\t')
            if line[1].replace('-','') == guide:
                if line[0] == 'X':
                    total_count_x[int(line[6])][int(line[7])] = total_count_x[int(line[6])][int(line[7])] + 1
                    if line[13] == 'y':
                        total_count_x_uniq[int(line[6])][int(line[7])] = total_count_x_uniq[int(line[6])][int(line[7])] + 1
                elif line[0] == 'DNA':
                    total_count_dna[int(line[6])][int(line[7])] = total_count_dna[int(line[6])][int(line[7])] + 1
                    if line[13] == 'y':
                        total_count_dna_uniq[int(line[6])][int(line[7])] = total_count_dna_uniq[int(line[6])][int(line[7])] + 1
                else:
                    total_count_rna[int(line[6])][int(line[7])] = total_count_rna[int(line[6])][int(line[7])] + 1
                    if line[13] == 'y':
                        total_count_rna_uniq[int(line[6])][int(line[7])] = total_count_rna_uniq[int(line[6])][int(line[7])] + 1

        
        # np.savetxt('summary_table_x.txt', total_count_x, delimiter='\t')
        # np.savetxt('summary_table_dna.txt', total_count_dna, delimiter='\t')
        # np.savetxt('summary_table_rna.txt', total_count_rna, delimiter='\t')
        tab_summary = pd.DataFrame(columns = ['Guide', 'Bulge Type', 'Bulge Size', 'Mismatches', 'Number of targets', 'Targets created by SNPs'])
        for m in range(mms + 1):
            for b_d in range(bulges_dna +1):
                tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'DNA', 'Bulge Size': b_d, 'Mismatches': m, 'Number of targets': total_count_dna[m][b_d], 'Targets created by SNPs':total_count_dna_uniq[m][b_d]  }, ignore_index = True)            

            for b_r in range(bulges_rna +1):
                tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'RNA', 'Bulge Size': b_r, 'Mismatches': m, 'Number of targets': total_count_rna[m][b_r], 'Targets created by SNPs': total_count_rna_uniq[m][b_r] }, ignore_index = True)

            tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'X', 'Bulge Size': 0, 'Mismatches': m, 'Number of targets': total_count_x[m][0], 'Targets created by SNPs':total_count_x_uniq[m][0] }, ignore_index = True)
    else:
        for line in targets:
            
            line = line.strip().split('\t')
            if line[1].replace('-','') == guide:
                if line[0] == 'X':
                    total_count_x[int(line[6])][int(line[7])] = total_count_x[int(line[6])][int(line[7])] + 1
                elif line[0] == 'DNA':
                    total_count_dna[int(line[6])][int(line[7])] = total_count_dna[int(line[6])][int(line[7])] + 1
                else:
                    total_count_rna[int(line[6])][int(line[7])] = total_count_rna[int(line[6])][int(line[7])] + 1
        
        # np.savetxt('summary_table_x.txt', total_count_x, delimiter='\t')
        # np.savetxt('summary_table_dna.txt', total_count_dna, delimiter='\t')
        # np.savetxt('summary_table_rna.txt', total_count_rna, delimiter='\t')
        tab_summary = pd.DataFrame(columns = ['Guide', 'Bulge Type', 'Bulge Size', 'Mismatches', 'Number of targets'])
        for m in range(mms + 1):
            for b_d in range(bulges_dna +1):
                tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'DNA', 'Bulge Size': b_d, 'Mismatches': m, 'Number of targets': total_count_dna[m][b_d] }, ignore_index = True)            

            for b_r in range(bulges_rna +1):
                tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'RNA', 'Bulge Size': b_r, 'Mismatches': m, 'Number of targets': total_count_rna[m][b_r] }, ignore_index = True)

            tab_summary =tab_summary.append({'Guide': guide, 'Bulge Type': 'X', 'Bulge Size': 0, 'Mismatches': m, 'Number of targets': total_count_x[m][0] }, ignore_index = True)

#print(tab_summary)
tab_summary.to_pickle(sys.argv[6] + '_summary_result_' + guide +'.txt')