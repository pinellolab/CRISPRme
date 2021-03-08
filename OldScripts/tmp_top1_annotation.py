'''
cambia l'annotazione Annotation.txt in un fac simile a Annotation.summary.txt, ma
contiene i conti solo del top1 presi da jobid.sample_annotation.guida
Il report è poi generato con questo var (uso radar_chart_new)
'''
# argv 1 is job id
# argv 2 is dir with sample annotation
import os 
import sys
from os import listdir                      #for getting directories
from os.path import isfile, isdir,join      #for getting directories
job_id = sys.argv[1]
ann_samp_file = [f for f in listdir(sys.argv[2] + '/') if isfile(join(sys.argv[2] + '/', f))]

summary_total = dict()
summary_per_guide = dict()
annotation_order = []
for f in ann_samp_file:
    if '.sample_annotation.' in f and '.population.' in f:
        guide = f.split('.sample_annotation.')[-1]
        guide = guide.split('.')[0]
        # summary_total[guide] = dict()
        summary_per_guide[guide] = dict()
        with open(f) as ann_samp:
            sum_tot = ann_samp.read().strip().split('-')[1].split('\n')[1:-1]
            # summary_per_guide[guide] = sum_tot
            for ann in sum_tot:
                tmp_ann = ann.split('\t')
                if tmp_ann[0] not in annotation_order:
                    annotation_order.append(tmp_ann[0])
                summary_per_guide[guide][tmp_ann[0]] = tmp_ann.copy()
                # print(tmp_ann)
                try:
                    summary_total[tmp_ann[0]] = [ int(summary_total[tmp_ann[0]][x]) + int(tmp_ann[x]) for x in range (1,11)]
                except:
                    summary_total[tmp_ann[0]] = tmp_ann.copy()

# print(summary_total)
# print(annotation_order)
with open(job_id + '.tmp_res.txt', 'w+') as result:
    result.write('-Summary_Total\n')
    # result.write('targets\t' + '\t'.join(str(x) for x in summary_total[tmp_ann['targets']]) + '\n')
    for pos, a in enumerate(annotation_order):
        # print(a)
        # print(summary_total[tmp_ann[a]])
        # print('tmp a',tmp_ann)
        # print('tmp a in pos',tmp_ann[pos])
        # print('sum tot a', [str(x) for x in summary_total[a][1:]])
        result.write(a + '\t' + '\t'.join([str(x) for x in summary_total[a]]) + '\n')
    for g in summary_per_guide:
        result.write('-Summary_' + g + '\n')
        for pos,a in enumerate(annotation_order):
            result.write(a + '\t' + '\t'.join(summary_per_guide[g][a][1:]) + '\n')

