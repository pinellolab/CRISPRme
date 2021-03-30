'''
Script that annotatets the samples, in order to have a fast generate-report
Input file is job_id.top_1.samples, job_id.Annotations.targets., job_id.Annotation.summary.txt, result name
Create a dict for the guides, that contains a dict for the samples, that contains a dict for the annotatio category
Eg 
{
    GUIDE1 -> {
        SAMPLE1 -> {
            EXON -> [0 0 0 0 1 0 0 8 2 0],
            INTRONS -> [0 0 0 0 1 0 0 8 2 0],
            CTCF -> [0 0 0 0 1 0 0 8 2 0]
        },
        SAMPLE2 ->{
            EXON -> [0 0 0 0 1 0 0 8 2 0],
            INTRONS -> [0 0 0 0 1 0 0 8 2 0],
            CTCF ->[0 0 0 0 1 0 0 8 2 0]
        }
    },
    GUIDE2 ->{
        SAMPLE ->{
            ANNOTATION -> [annotation count for each mm value]
        }
    }
}
'''
# argv 1 is top1.samples.txt
# argv 2 is Annotation.targets
# argv 3 is Annotation.summary.txt -> to get the name of annotations    #TODO modificare meglio, prenderle direttamente dal file
# argv 4 is result name
import sys
import os
import pandas as pd

#Dict for populations
pop_file = pd.read_excel(os.path.dirname(os.path.realpath(__file__)) + '/20130606_sample_info.xlsx')
all_samples = pop_file.Sample.to_list()
all_pop = pop_file.Population.to_list()
dict_pop = dict()
for  pos, i in enumerate(all_samples):
    try:
        dict_pop[i] = all_pop[pos]        #{'S1':'POP1', 'S2':'POP1', ...}
    except:
        dict_pop[i] = all_pop[pos]

#Dict for superpopulation
population_1000gp = {'CHB':'EAS', 'JPT':'EAS', 'CHS':'EAS', 'CDX':'EAS', 'KHV':'EAS',
                    'CEU':'EUR', 'TSI':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR',
                    'YRI':'AFR', 'LWK':'AFR', 'GWD':'AFR', 'MSL':'AFR', 'ESN':'AFR', 'ASW':'AFR', 'ACB':'AFR',
                    'MXL':'AMR', 'PUR':'AMR', 'CLM':'AMR', 'PEL':'AMR',
                    'GIH':'SAS', 'PJL':'SAS', 'BEB':'SAS', 'STU':'SAS', 'ITU':'SAS'
}
superpopulation = ['EAS', 'EUR', 'AFR', 'AMR','SAS']

result_name = sys.argv[4]
# samples_dict = {
    # GUIDE1 ->{
    #     chrXposY -> [[Sample1, sample7], []]
    #     chrXposY2 -> [[Sample5, sample7], []]
    #     chrX2posY 3-> [[Sample10, sample11, sample30], []]
    # },
    # GUIDE2 -> {
    #     CHRPOS -> [[Sample list], [List visited annotations]]                List visited annotation is empty at first, but can become -> ['exon', 'promoter',...]
    # }
# }
test_dict = {'GAGTCCGAGCAGAAGAAGAANNN':{0:0,1:0,2:0,3:0,4:0,5:0,6:0}, 'CCATCGGTGGCCGTTTGCCCNNN':{0:0,1:0,2:0,3:0,4:0,5:0,6:0}}

samples_dict = dict()
annotation_dict = dict()
with open(sys.argv[1]) as targets:
    for line in targets:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        if line[-2] == 'n':
            test_dict[line[1].replace('-','')][int(line[7])] +=1
            continue
        guide = line[1].replace('-','')
        if guide not in samples_dict:
            samples_dict[guide] = dict()
        try:
            samples_dict[guide][line[3] + line[4]][0] += line[-2].split(',')
        except:     
            samples_dict[guide][line[3] + line[4]] = [line[-2].split(','), []]
print(test_dict)
# print(samples_dict)
# print(samples_dict['CTAACAGTTGCTTTTATCACNNN']['chr2146560428'])
# print(samples_dict['TGCTTGGTCGGCACTGATAGNNN']['chr2250085897'])

ann_list = []       #TODO better way to get annotation name

with open (sys.argv[3], 'r') as ann_file:
    next(ann_file) #Skip -Summary_Total line
    next(ann_file) #Skip targets 0 0 0 ... line
    for line in ann_file:
        if '-Summary' in line:
            break
        ann_list.append(line.strip().split('\t')[0])
        


summary_targets_guide = dict()  #To save targets and annotation count for top 1
dict_pop_count = dict()         #To save targets and annotations count for populations
dict_superpop_count = dict()    #To save targets and annotations count for superpopulations

with open (sys.argv[2]) as targets:             #Count annotation for each target
    for line in targets:
        if '#' in line:
            continue
        line = line.strip().split('\t')
        guide = line[1].replace('-','')
        
        if guide not in annotation_dict.keys():
            annotation_dict[guide] = dict()
            summary_targets_guide[guide] = {'targets':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
            dict_pop_count[guide] = dict()
            dict_superpop_count[guide] = dict()
            for pop in set(all_pop):
                dict_pop_count[guide][pop] = {'targets':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
                dict_superpop_count[guide][population_1000gp[pop]] = {'targets':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
        
        # try:
        #     summary_targets_guide[guide][line[-1]][int(line[7])] += 1
        # except:
        #     summary_targets_guide[guide][line[-1]] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        #     summary_targets_guide[guide][line[-1]][int(line[7])] += 1
        # summary_targets_guide[guide]['targets'][int(line[7])] += 1
        try:
            samples_list = samples_dict[guide][line[3] + line[4]]
        except:
            samples_list = [[], ann_list]
        if line[-1] in samples_list[1]: #if target was already counted in that annotation
            continue
        #Count the annotations for the guide (only top1)
        try:
            summary_targets_guide[guide][line[-1]][int(line[7])] += 1
        except:
            summary_targets_guide[guide][line[-1]] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            summary_targets_guide[guide][line[-1]][int(line[7])] += 1
        summary_targets_guide[guide]['targets'][int(line[7])] += 1
        
        samples_dict[guide][line[3] + line[4]][1].append(line[-1])  #Get list of samples in current gudie chr pos
        visited_pop = []
        visited_superpop = []
        for sample in samples_list[0] :
            if sample not in annotation_dict[guide]:
                annotation_dict[guide][sample] = {'targets':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}     
                # print(guide, sample, line[-1], line[6])
                
            try:
                annotation_dict[guide][sample][line[-1]][int(line[7])] += 1       #increase annotation count
            except:
                annotation_dict[guide][sample][line[-1]] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                annotation_dict[guide][sample][line[-1]][int(line[7])] += 1
            
            if dict_pop[sample] in visited_pop:
                continue
            else:
                visited_pop.append(dict_pop[sample])
                dict_pop_count[guide][dict_pop[sample]]['targets'][int(line[7])] += 1
                try:
                    dict_pop_count[guide][dict_pop[sample]][line[-1]][int(line[7])] += 1
                except:
                    dict_pop_count[guide][dict_pop[sample]][line[-1]] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    dict_pop_count[guide][dict_pop[sample]][line[-1]][int(line[7])] += 1
            if population_1000gp[dict_pop[sample]] in visited_superpop:
                continue
            else:
                visited_superpop.append(population_1000gp[dict_pop[sample]])
                dict_superpop_count[guide][population_1000gp[dict_pop[sample]]]['targets'][int(line[7])] += 1
                try:
                    dict_superpop_count[guide][population_1000gp[dict_pop[sample]]][line[-1]][int(line[7])] += 1
                except:
                    dict_superpop_count[guide][population_1000gp[dict_pop[sample]]][line[-1]] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    dict_superpop_count[guide][population_1000gp[dict_pop[sample]]][line[-1]][int(line[7])] += 1
            annotation_dict[guide][sample]['targets'][int(line[7])] += 1

for guide in annotation_dict:
    with open(result_name + '.sample_annotation.' + guide +'.samples.txt', 'w+') as result:
        result.write('-Summary_Total\n')
        result.write('targets\t' + '\t'.join([str(x) for x in summary_targets_guide[guide]['targets']]) + '\n')
        for annotation in ann_list:
            try:
                result.write(annotation + '\t' + '\t'.join([str(x) for x in summary_targets_guide[guide][annotation]]) + '\n')
            except:
                result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
        
        #Write summary for each sample
        for sample in all_samples:#annotation_dict[guide]:
            result.write('-Summary_' + sample + '\n')
            try:
                result.write('targets\t' + '\t'.join([str(x) for x in annotation_dict[guide][sample]['targets']]) + '\n')
            except:
                result.write('targets' + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
            for annotation in ann_list:
                try:
                    result.write(annotation + '\t' + '\t'.join([str(x) for x in annotation_dict[guide][sample][annotation]]) + '\n')
                except:
                    result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')


for guide in annotation_dict:
    with open(result_name + '.sample_annotation.' + guide + '.population.txt', 'w+') as result:
        #Write result guide general
        result.write('-Summary_Total\n')
        result.write('targets\t' + '\t'.join([str(x) for x in summary_targets_guide[guide]['targets']]) + '\n')
        for annotation in ann_list:
            try:
                result.write(annotation + '\t' + '\t'.join([str(x) for x in summary_targets_guide[guide][annotation]]) + '\n')
            except:
                result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
        #Write result population
        for population in set(all_pop):
            result.write('-Summary_' + population + '\n')
            try:
                result.write('targets\t' + '\t'.join([str(x) for x in dict_pop_count[guide][population]['targets']]) + '\n')
            except:
                result.write('targets' + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
            for annotation in ann_list:
                try:
                    result.write(annotation + '\t' + '\t'.join([str(x) for x in dict_pop_count[guide][population][annotation]]) + '\n')
                except:
                    result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
  

#For each superpopulation, write sum of population
for guide in annotation_dict:    
    with open(result_name + '.sample_annotation.' + guide + '.superpopulation.txt', 'w+') as result:
        #Write result guide general
        result.write('-Summary_Total\n')
        result.write('targets\t' + '\t'.join([str(x) for x in summary_targets_guide[guide]['targets']]) + '\n')
        for annotation in ann_list:
            try:
                result.write(annotation + '\t' + '\t'.join([str(x) for x in summary_targets_guide[guide][annotation]]) + '\n')
            except:
                result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
        #Write result superpopulation
        for superpop in superpopulation:
            result.write('-Summary_' + superpop + '\n')
            try:
                result.write('targets\t' + '\t'.join([str(x) for x in dict_superpop_count[guide][superpop]['targets']]) + '\n')
            except:
                result.write('targets' + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
            for annotation in ann_list:
                try:
                    result.write(annotation + '\t' + '\t'.join([str(x) for x in dict_superpop_count[guide][superpop][annotation]]) + '\n')
                except:
                    result.write(annotation + '\t' + '\t'.join(['0' for i in range(10)]) + '\n')
            

            