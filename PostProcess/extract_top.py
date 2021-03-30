#Script extract only the top X from all clusters
# argv1 is input file ordered in cluster
# argv2 is job_id 
#NOTE 06/03  -> removed PAM Disruption calculation
import sys
header = {'Bulge Type':0, '#Bulge Type':0, 'Bulge type':0,
        'crRNA':1, 'DNA':2,
        'Chromosome':3, 'Position':4, 'Cluster Position':5 ,'Direction':6,
        'Mismatches':7, 'Bulge Size':8, 'Total':9,'Min_mismatches':10, 'Max_mismatches':11,
        #'PAM_disr':12,
        'PAM_gen':12, 'Var_uniq':13, 'Samples':14}

top_x = 1
top_x_inserted = 0
job_id = sys.argv[2]
with open(sys.argv[1]) as targets, open( job_id + '.top_' + str(top_x) + '.txt', 'w+') as result:
    line = targets.readline()
    if '#' in line:
        result.write(line)  #header
        line = targets.readline().strip().split('\t')
    else:
        line = line.strip().split('\t')
    prev_position = line[header['Cluster Position']]
    prev_chr = line[header['Chromosome']]
    result.write('\t'.join(line) + '\n')
    top_x_inserted = top_x_inserted + 1
    for line in targets:
        line = line.strip().split('\t')
        current_pos = line[header['Cluster Position']]
        current_chr = line[header['Chromosome']]
        #if i can still instert top1
        if prev_position == current_pos and prev_chr == current_chr and top_x_inserted < top_x:
            top_x_inserted = top_x_inserted + 1
            result.write('\t'.join(line) + '\n')
        elif prev_position != current_pos or prev_chr != current_chr:
            top_x_inserted = 1
            prev_position = current_pos
            prev_chr = current_chr
            result.write('\t'.join(line) + '\n')