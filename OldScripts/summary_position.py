#OLD, VIEW summary_by_guide_position.py FOR NEW CHANGES


#Script for the Summary by Position table. Given in input a targets.txt ordered by clusters, for each cluster (position) it saves: 
# -chr
# -pos
# -best target sequence
# -min mms of cluster
# -min bulge of cluster
# -Total targets with 0 mm 0 bulges
# -Total targets with 1 mm 0 bulges
# ...
# -Total targets with n mm 0 bulges
# -Total targets with 0 mm 1 bulges
# -Total targets with 1 mm 1 bulges
# ...
# -Total targets with n mm 1 bulges
# ...
# -Total targets with n mm b bulges

# argv1 is targets.txt ordered in clusters
# argv2 is job_id
# argv3 is guide
# argv4 is mms
# argv5 is bulge (max between DNA and RNA)
import sys

job_id = sys.argv[2]
guide = sys.argv[3]
job_id = '0YT6LD1ECN' #TODO cancellare, è solo temporaneo per i test
guide = 'CTAACAGTTGCTTTTATCACNNN' #TODO cancellare, è solo temporaneo per i test
mms = int(sys.argv[4])
bulge = int(sys.argv[5])
#cluster_count_mms = [0 for i in range (mms + 1)]
#cluster_count_bulges = [0 for i in range (bulge)]   #NOTE pos 0 is 1 bulge, pos 1 is 2bulges ...
#cluster_count_bulges = [[0 for i in range (mms + 1)] for i in range (bulge)] #NOTE pos 0 is 1 bulge, pos 1 is 2bulges ...
count_targets = [[0 for i in range (mms + 1)] for i in range (bulge + 1)]
with open(sys.argv[1]) as targets, open(job_id + '.summary_position.' + guide + '.txt', 'w+') as result:
    result.write('#Chromosome\tPosition\tBest Target\tMin Mismatch\tMin Bulge')
    for b in range(bulge + 1):
        for i in range(mms + 1):
            result.write('\tTargets ' + str(i) + 'MM' + str(b) + 'B' )
    # for i in range(bulge):
    #     result.write('\tTargets ' + str(i + 1) + ' Bulge')
    result.write('\n')

    for line in targets:
        line = line.strip().split('\t')
        if line[1].replace('-','') == guide:
            break    
    # line = targets.readline()
    # if '#' in line:
    #     line = targets.readline().strip().split('\t')
    # else:
    #     line = line.strip().split('\t')
     
    current_cluster = line[-1] #NOTE with new version of clusterization, must use line['Cluster Position'] and not line['Position']
    result.write(line[3] + '\t' + line[-1] + '\t' + line[2] + '\t' + line[6] + '\t' + line[7])
    mms_current_line = int(line[6])
    bulge_current_line = int(line[7])
    #cluster_count_mms[mms_current_line] = cluster_count_mms[mms_current_line] + 1
    #cluster_count_bulges[bulge_current_line - 1] = cluster_count_bulges[bulge_current_line - 1] + 1
    count_targets[bulge_current_line][mms_current_line] = count_targets[bulge_current_line][mms_current_line] + 1 
    for line in targets:
        
        line = line.strip().split('\t')
        if line[1].replace('-','') != guide:
            continue
        
        if current_cluster == line[-1]:
            mms_current_line = int(line[6])
            bulge_current_line = int(line[7])
            #cluster_count_mms[mms_current_line] = cluster_count_mms[mms_current_line] + 1
            #cluster_count_bulges[bulge_current_line - 1] = cluster_count_bulges[bulge_current_line - 1] + 1
            count_targets[bulge_current_line][mms_current_line] = count_targets[bulge_current_line][mms_current_line] + 1 
        else:
            # for m in  cluster_count_mms:
            #     result.write('\t' + str(m))
            #     cluster_count_mms = [0 for i in range (mms + 1)]
            # for b in  cluster_count_bulges:
            #     result.write('\t' + str(b))
            #     cluster_count_bulges = [0 for i in range (bulge)]
            for t in count_targets:
                for t_c in t:
                    result.write('\t' + str(t_c))
                count_targets = [[0 for i in range (mms + 1)] for i in range (bulge + 1)]
            result.write('\n')
            current_cluster = line[-1]
            result.write(line[3] + '\t' + line[-1] + '\t' + line[2] + '\t' + line[6] + '\t' + line[7])
            mms_current_line = int(line[6])
            bulge_current_line = int(line[7])
            #cluster_count_mms[mms_current_line] = cluster_count_mms[mms_current_line] + 1
            #cluster_count_bulges[bulge_current_line - 1] = cluster_count_bulges[bulge_current_line - 1] + 1
            count_targets[bulge_current_line][mms_current_line] = count_targets[bulge_current_line][mms_current_line] + 1 
    #Write result for last cluster
    # for m in  cluster_count_mms:
    #     result.write('\t' + str(m))
    # for b in  cluster_count_bulges:
    #     result.write('\t' + str(b))
    for t in count_targets:
        for t_c in t:
            result.write('\t' + str(t_c))
    result.write('\n')