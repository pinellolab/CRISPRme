#!/bin/sh
#NOTE ALL PATH ARE NOW ABSOLUTE
#$1 is directory of result for submitted job id (Results/jobid)
#$2 is genome_selected directory Eg Genomes/hg19_ref or Genomes/hg19_ref+1000genomeproject
#$3 is genome_ref directory     Eg Genomes/hg19_ref
#$4 is genome_idx directory     Eg genome_library/NGG_2_hg19_ref+hg19_1000genomeproject or genome_library/NGG_2_hg19_ref
#$5 is pam file             Eg Results/72C1MNXDWF/pam.txt
#$6 is guides file          Eg Results/LNHM6F3REO/guides.txt (both upload file and custom inserted)
#$7 is mms
#$8 is dna
#$9 is rna
#$10 is search_index
#$11 is search
#$12 is annotation
#$13 is generate report
#$14 is gecko comparison
#$15 is genome_ref comparison
#$16 is genme_idx_ref (for genome_ref comparison if search was done with indices)   Eg genome_library/NGG_2_hg19_ref
#$17 is send_email
#$18 is annotation file EG annotations/hg19_ref.annotations.bed
#$19 is genome type, can be 'ref', 'var', 'both'
#$20 is absolute path of directory containing the dash app
#$21 is absolute path of dictionary directory NOTE the dictionary is already specified for the selected genome
#$22 is absolute path of sample file for extracting dictionaries sample -> pop, pop-> superpop etc (using supportFunction associate sample)
#$23 is generate index ref, True or False
#$24 is generate index enr, True or False
#$25 is directory where python was called (directory where all data Genomes, Results etc is located)
#Note that if genome_selected is not enriched, the python exe will force $15 as false

cd $1
rm queue.txt
jobid=$(basename $1)
name_reference=$(basename $3)
name_enriched=$(basename $2)
dictionaries=${21}
echo 'START '$jobid ${18}
echo 'Job\tStart\t'$(date)> log.txt
used_genome_dir=$2                 
max_bulge=$([ $8 -ge $9 ] && echo "$8" || echo "$9")
total=$(($7+max_bulge))

#Start indexing reference genome
if [ ${23} = 'True' ]; then
    echo 'Indexing_Reference' >> output.txt 
    cd ${25}
    crispritz.py index-genome $name_reference $3 $1/pam_indexing.txt -bMax $max_bulge
    cd $1
fi
#Start indexing enr genome
if [ ${24} = 'True' ]; then
    echo 'Indexing_Enriched' >> output.txt
    cd ${25}
    crispritz.py index-genome $name_enriched $2 $1/pam_indexing.txt -bMax $max_bulge
    cd $1
fi

#Start search index     #NOTE new version 2.1.2 of crispritz needed
echo 'Search-index\tStart\t'$(date) >> log.txt
echo 'Search_output '${19} >  output.txt
if [ ${10} = 'True' ]; then
    #echo 'crispritz search-index'
    crispritz.py search $4 pam.txt guides.txt $jobid -mm $7 -bDNA $8 -bRNA ${9} -t #>> output.txt #TODO sistemare l'output redirection


    if [ ${15} = 'True' ]; then
        mkdir 'ref'
        echo 'Search_output_ref '${19} >>  output.txt
        crispritz.py search ${16} pam.txt guides.txt $jobid'_ref' -mm $7 -bDNA $8 -bRNA ${9} -t #>> output.txt #TODO sistemare l'output redirection
        mv ./$jobid'_ref'.*.txt 'ref'
        mv ./$jobid'_ref'.*.xls 'ref'
    fi
fi
echo 'Search-index\tDone\t'$(date) >> log.txt

#Start search
echo 'Search\tStart\t'$(date) >> log.txt
if [ ${11} = 'True' ]; then
    #echo 'crispritz search'
    crispritz.py search $used_genome_dir pam.txt guides.txt $jobid -mm $7 -var -t -th 8 #>> output.txt #-scores $3
    
    if [ ${15} = 'True' ]; then
        mkdir 'ref'
        echo 'Search_output_ref '${19} >>  output.txt 
        crispritz.py search $3 pam.txt guides.txt $jobid'_ref' -mm $7 -var -t -th 8 #>> output.txt
        mv ./$jobid'_ref'.*.txt 'ref'
        mv ./$jobid'_ref'.*.xls 'ref'
    fi
fi
echo 'Search\tDone\t'$(date) >> log.txt

#Data processing + annotation + generating report
if [ ${19} = 'ref' ]; then
    echo 'PostProcess\tStart\t'$(date) >> log.txt
    echo 'PostProcess_output' > output.txt
    echo 'Clustering... Step [1/3]' >>  output.txt
    echo 'Start cluster ref'
    python ${20}PostProcess/cluster.dict.py $jobid.targets.txt 'addGuide' 'True' 'False' guides.txt # > $jobid.targets.cluster.txt
    echo 'End cluster ref'

    echo 'Start extraction top1 ref'
    python ${20}PostProcess/extract_top.py $jobid.targets.cluster.txt $jobid  # > $jobid.top_1.txt
    echo 'End extraction top1 ref'

    echo 'Start Scoring'
    echo 'Calculating Scores... Step [2/3]' >>  output.txt
    python ${20}PostProcess/scores_guide_table.py $jobid.top_1.txt $3 pam.txt guides.txt
    echo 'End Scoring'

    echo 'Start summary by pos-guide'
    echo 'Creating Summaries... Step [3/3]' >>  output.txt
    type_post='No'      
    python ${20}PostProcess/summary_by_guide_position.py $jobid.targets.cluster.txt $7 $8 $9 guides.txt $jobid $type_post
    echo 'End summary by pos-guide'
    

    echo 'Annotation\tStart\t'$(date) >> log.txt
    # echo 'Annotate_output '${19} >  output.txt
    echo 'Annotating...' >>  output.txt
    echo 'Start annotation ref'
    if [ ${12} = 'True' ]; then
        crispritz.py annotate-results $jobid.top_1.txt ${18} $jobid >> output.txt # > $jobid.Annotation.summary.txt , $jobid.Annotation.targets.txt 
    fi
    echo 'End annotation ref'
    echo 'Annotation\tDone\t'$(date) >> log.txt
    echo 'PostProcess\tDone\t'$(date) >> log.txt
    
    #Start generate report
    echo 'Report\tStart\t'$(date) >> log.txt
    if [ ${13} = 'True' ]; then
        echo 'Generate_report' >  output.txt
        proc=$(($7 + 1))
        while IFS= read -r line || [ -n "$line" ]; do    
            if [ ${14} = 'True' ]; then         #If -gecko
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref ref/$jobid'_ref'.Annotation.summary.txt  -gecko -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -gecko -ws
                fi
            else
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref ref/$jobid'_ref'.Annotation.summary.txt -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -ws
                fi
            fi

            #Create a zip file, for each guide, containing the targets from jobid.samples.annotation.txt
            #Grep all the targets with the selected guide
            touch $jobid.targets.$line.txt
            head -1 $jobid.Annotation.targets.txt > $jobid.targets.$line.txt        #Put header
            LC_ALL=C grep $line $jobid.Annotation.targets.txt >> $jobid.targets.$line.txt

            zip $jobid.targets.$line.zip $jobid.targets.$line.txt
            rm $jobid.targets.$line.txt  #NOTE this file could be used for grepping when 'Show Targets' is selected

            echo $line >> output.txt
        done < guides.txt

        mkdir ${20}assets/Img/$jobid
        cp *.png ${20}assets/Img/$jobid/
        echo 'Report\tDone\t'$(date) >> log.txt
    fi

elif [ ${19} = 'var' ]; then
    cp ${22} ./sampleID.txt     #Copy sampleId list file to result directory
    echo 'PostProcess\tStart\t'$(date) >> log.txt
    echo 'Start cluster var'
    echo 'PostProcess_output' > output.txt
 
    echo 'Start pam analysis'
    echo 'PAM Analysis... Step [1/4]' >>  output.txt
    #Add pam creation, variant unique, samples, annotation, real guide columns
    awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide}' $jobid.targets.txt > $jobid.targets.minmaxdisr.txt
    # python ../../PostProcess/pam_analysis.py $jobid.targets.txt pam.txt 'both'     # 'both' to add PAM creation and Var Unique column -> set to n
                # > $jobid.targets.minmaxdisr.txt
    echo 'End pam analysis'

    echo 'Clustering... Step [2/4]' >>  output.txt
    #Clustering: order is by chr, cluster are not ordered among each other
    python ${20}PostProcess/cluster.dict.py $jobid.targets.minmaxdisr.txt 'no' 'True' 'False' guides.txt 'total' 'orderChr' # > $jobid.targets.minmaxdisr.cluster.txt
    echo 'End cluster var'
    
    echo 'Extracting Samples and Annotation... (This operation has a long execution time, Please Wait) Step [3/4]' >>  output.txt
    echo 'Start calc samples and annotation and scores'
    echo 'Annotation\tStart\t'$(date) >> log.txt
    python ${20}PostProcess/annotator_for_onlyvar.py ${18} $jobid.targets.minmaxdisr.cluster.txt $jobid $dictionaries pam.txt $7 $2 guides.txt $8 $9 ${22}
        # > $jobid.samples.all.annotation.txt with header AGGIORNAMENTO 11/03 QUESTO FILE NON VIENE CREATO 
        # > $jobid.samples.annotation.txt  with header AGGIORNAMENTO 11/03 Contiene top1 scomposti e top1 reference (usato per sum guide e show target guide, sample)
        # > $jobid.Annotation.summary.txt
        # > $jobid.sample_annotation.GUIDE.sample.txt
        # > $jobid.sumref.Annotation.summary.txt
        # > $jobid.cluster.tmp.txt AGGIORNAMENTO Top1 sostituito col min mms scomposto, il resto del cluster ha ancora IUPAC
        # > acfd.txt
    mv $jobid.cluster.tmp.txt $jobid.total.cluster.txt   #Now has sample and annotation (for top1, for other only blank column)

    echo 'End calc samples and annotation and scores'
    echo 'Annotation\tDone\t'$(date) >> log.txt
    
    echo 'Creating Summaries... Step [4/4]' >>  output.txt
    #Summary guide, pos #NOTE the script automatically counts only for top subclusters
    echo 'Start summary by guide and position'  #NOTE summary by guide will be overwritten
    python ${20}PostProcess/summary_by_guide_position.py $jobid.total.cluster.txt $7 $8 $9 guides.txt $jobid 'Uniq'  
    echo 'End summary by guide and position' 


    #Summary guide
    echo 'Start summary by guide'  
    python ${20}PostProcess/summary_by_guide.py $jobid.samples.annotation.txt $7 $8 $9 guides.txt $jobid 'no'
    echo 'End summary by guide'

    #Summary samples
    echo 'Start summary by samples'
    python ${20}PostProcess/summary_by_samples.py $jobid.samples.annotation.txt $jobid 'both' guides.txt ${22}
    #python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt 
    echo 'End summary by samples'

    echo 'PostProcess\tDone\t'$(date) >> log.txt

    #Start generate report
    echo 'Report\tStart\t'$(date) >> log.txt
    if [ ${13} = 'True' ]; then
        echo 'Generate_report' >  output.txt
        proc=$(($7 + 1))
        while IFS= read -r line || [ -n "$line" ]; do    
            if [ ${14} = 'True' ]; then         #If -gecko
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref ref/$jobid'_ref'.Annotation.summary.txt  -gecko -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -gecko -ws
                fi
            else
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref ref/$jobid'_ref'.Annotation.summary.txt -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -ws
                fi
            fi
            #Generate Population Distributions
            printf %s\\n $(seq 0 $total) | xargs -n 1 -P $proc -I % python ${20}PostProcess/populations_distribution.py $jobid.PopulationDistribution.txt % $line
            
            #Create a zip file, for each guide, containing the targets from jobid.samples.annotation.txt
            #Grep all the targets with the selected guide
            head -1 $jobid.samples.annotation.txt > $jobid.targets.$line.txt        #Put header
            LC_ALL=C grep $line $jobid.samples.annotation.txt >> $jobid.targets.$line.txt

            #Clusterize results
            python ${20}PostProcess/cluster.dict.py $jobid.targets.$line.txt 'no' 'True' 'True' guides.txt 'total' 'addForFinal'    #> $jobid.targets.$line.cluster.txt

            mv $jobid.targets.$line.cluster.txt $jobid.targets.$line.txt

            zip $jobid.targets.$line.zip $jobid.targets.$line.txt
            rm $jobid.targets.$line.txt  #NOTE this file could be used for grepping when 'Show Targets' is selected
            
            echo $line >> output.txt
        done < guides.txt
        mkdir ${20}assets/Img/$jobid
        cp *.png ${20}assets/Img/$jobid/
        echo 'Report\tDone\t'$(date) >> log.txt
    fi

else    #Type search = both
    cp ${22} ./sampleID.txt     #Copy sampleId list file to result directory
    echo 'Report\tStart\t'$(date) >> log.txt
    #Estract common, semicommon and unique
    echo 'Start creation semicommon, common, unique'
    echo 'PostProcess_output' > output.txt
    echo 'Processing Search Results... Step [1/6]' >>  output.txt
    ${20}PostProcess/./extraction.sh ref/$jobid'_ref.targets.txt' $jobid.targets.txt $jobid # > $jobid.common_targets.txt $jobid.semi_common_targets.txt $jobid.unique_targets.txt
    echo 'End creation semicommon, common, unique'                                          
    
    #Cluster semicommon (that also contains common -> see extraction.sh) e uniq
    # echo 'Start cluster semicommon'
    # echo 'Clustering... Step [2/6]' >>  output.txt
    # python ../../PostProcess/cluster.dict.py $jobid.semi_common_targets.txt 'no' 'True' 'False' guides.txt     #solo per aggiungere colonna total -> forse si può togliere
    # echo 'End cluster semicommon'
    # echo 'Start cluster unique'     #NOTE doing cluster separately does not create the right order of cluster (first clusters of uniq, then clusters of semi_common)
    # python ../../PostProcess/cluster.dict.py $jobid.unique_targets.txt 'no' 'True' 'False' guides.txt      #solo per aggiungere colonna total -> forse si può togliere
    # echo 'End cluster unique'
    #Colonna total già aggiunta dal search

    #Pam analysis
    echo 'Start pam analysis'
    echo 'PAM Analysis... Step [2/5]' >>  output.txt
    #Add pam creation, variant unique, samples, annotation, real guide columns
    awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\tn\tn\tn\t"real_guide}' $jobid.semi_common_targets.txt > $jobid.semi_common_targets.minmaxdisr.txt 
    # python ../../PostProcess/pam_analysis.py $jobid.semi_common_targets.txt pam.txt ${19}  # > $jobid.semi_common_targets.minmaxdisr.txt
    echo 'End pam analysis'
    echo 'Start pam creation'
    #Add pam creation, variant unique, samples, annotation, real guide columns
    awk '{real_guide=$2; gsub("-","",real_guide); print $0"\tn\ty\tn\tn\t"real_guide}' $jobid.unique_targets.txt > $jobid.unique_targets.pamcreation.txt #Add pam creation, variant unique, real guide column
    # python ../../PostProcess/pam_creation.py $jobid.unique_targets.txt pam.txt ../../$3 # > $jobid.unique_targets.pamcreation.txt
    echo 'End pam creation'
    cat $jobid.unique_targets.pamcreation.txt $jobid.semi_common_targets.minmaxdisr.txt > $jobid.total.txt

    #Cluster of jobid.total.txt and extraction of top 1
    echo 'Start cluster of total.txt'
    echo 'Clustering... Step [3/5]' >>  output.txt
    python ${20}PostProcess/cluster.dict.py $jobid.total.txt 'no' 'True' 'True' guides.txt 'total' 'orderChr'  #-> 03/03 il cluster ottenuto è ordinato per chr, uso questo per annotation
    echo 'End cluster of total.txt'
    # echo 'Start extract top1 total.txt'
    # # python ../../PostProcess/extract_top.py $jobid.total.cluster.txt $jobid # > $jobid.top_1.txt
    
    # #03/03 commentato riga sotto 
    # # awk '{guide=$2;gsub("-","",guide); print $0"\t"guide}' $jobid.total.cluster.txt > tmp_add_guide && mv tmp_add_guide $jobid.total.cluster.txt  #add real guide column for show targets on sum by pos

    # echo 'End extract top1 total.txt'

    #Scoring of top1
    # echo 'Start Scoring'
    # # python ../../PostProcess/scores_guide_table.py $jobid.top_1.txt ../../$used_genome_dir pam.txt guides.txt  #TODO da calcolare solo su target esistenti
    # echo 'End Scoring'


    #Top1 expansion
    # echo 'Start sort'
    echo 'Extracting Samples and Annotation... (This operation has a long execution time, Please Wait) Step [4/5]' >>  output.txt
    # sort -k4,4 $jobid.top_1.txt > $jobid.top_1.sort.txt && mv $jobid.top_1.sort.txt $jobid.top_1.txt 
    # echo 'End sort'
    echo 'Start calc samples and annotation and scores'
    echo 'Annotation\tStart\t'$(date) >> log.txt
    #03/03 modificato da top_1 a total.cluster
    python ${20}PostProcess/annotator_cal_sample.py ${18} $jobid.total.cluster.txt $jobid $dictionaries pam.txt $7 $3 guides.txt $8 $9 ${22}
        # > $jobid.samples.all.annotation.txt with header AGGIORNAMENTO 11/03 QUESTO FILE NON VIENE CREATO 
        # > $jobid.samples.annotation.txt  with header AGGIORNAMENTO 11/03 Contiene top1 scomposti e top1 reference (usato per sum guide e show target guide, sample)
        # > $jobid.Annotation.summary.txt
        # > $jobid.sample_annotation.GUIDE.sample.txt
        # > $jobid.sumref.Annotation.summary.txt
        # > $jobid.cluster.tmp.txt AGGIORNAMENTO Top1 sostituito col min mms scomposto, il resto del cluster ha ancora IUPAC
        # > acfd.txt
    mv $jobid.cluster.tmp.txt $jobid.total.cluster.txt   #Now has sample and annotation (for top1, for other only blank column)

    echo 'End calc samples and annotation and scores'
    echo 'Annotation\tDone\t'$(date) >> log.txt
    #Put right header into top_1.samples.all.txt
    # sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tPam_disr\tPAM_gen\tVar_uniq\tSamples\tReal Guide' $jobid.top_1.samples.all.txt
    
    
    echo 'Creating Summaries... Step [5/5]' >>  output.txt
    #Summary guide, pos #NOTE the script automatically counts only for top subclusters
    echo 'Start summary by guide and position'  #NOTE summary by guide will be overwritten
    python ${20}PostProcess/summary_by_guide_position.py $jobid.total.cluster.txt $7 $8 $9 guides.txt $jobid 'Uniq'  
    echo 'End summary by guide and position'


    #Summary guide
    echo 'Start summary by guide'  
    python ${20}PostProcess/summary_by_guide.py $jobid.samples.annotation.txt $7 $8 $9 guides.txt $jobid 'Uniq'
    echo 'End summary by guide'

    #Summary samples
    echo 'Start summary by samples'
    python ${20}PostProcess/summary_by_samples.py $jobid.samples.annotation.txt $jobid ${19} guides.txt ${22}
    #python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt 
    echo 'End summary by samples'

    echo 'PostProcess\tDone\t'$(date) >> log.txt

    #Generate report
    echo 'START Generate Report and Population Distribution'
    echo 'Report\tStart\t'$(date) >> log.txt
    if [ ${13} = 'True' ]; then
        echo 'Generate_report' >  output.txt
        proc=$(($7 + 1))
        while IFS= read -r line || [ -n "$line" ]; do    
            if [ ${14} = 'True' ]; then         #If -gecko
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref $jobid.sumref.Annotation.summary.txt  -gecko -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -gecko -ws
                fi
            else
                if [ ${15} = 'True' ]; then     #If genome_ref comparison
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -sumref $jobid.sumref.Annotation.summary.txt -ws
                else
                    printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -annotation $jobid'.Annotation.summary.txt' -extprofile *.extended_profile.xls -ws
                fi
            fi
            
            #Generate Population Distributions
            printf %s\\n $(seq 0 $total) | xargs -n 1 -P $proc -I % python ${20}PostProcess/populations_distribution.py $jobid.PopulationDistribution.txt % $line
            
            #Create a zip file, for each guide, containing the targets from jobid.samples.annotation.txt
            #Grep all the targets with the selected guide
            head -1 $jobid.samples.annotation.txt > $jobid.targets.$line.txt        #Put header
            LC_ALL=C grep $line $jobid.samples.annotation.txt >> $jobid.targets.$line.txt

            #Clusterize results
            python ${20}PostProcess/cluster.dict.py $jobid.targets.$line.txt 'no' 'True' 'True' guides.txt 'total' 'addForFinal'    #> $jobid.targets.$line.cluster.txt

            mv $jobid.targets.$line.cluster.txt $jobid.targets.$line.txt

            zip $jobid.targets.$line.zip $jobid.targets.$line.txt
            rm $jobid.targets.$line.txt  #NOTE this file could be used for grepping when 'Show Targets' is selected

            echo $line >> output.txt
        done < guides.txt
        mkdir ${20}assets/Img/$jobid
        cp *.png ${20}assets/Img/$jobid/
    fi
    echo 'Report\tDone\t'$(date) >> log.txt

    echo 'END Generate Report and Population Distribution'

fi

cd ../../
if [ ${17} = 'True' ]; then
    python ${20}send_mail.py $1
fi
echo 'Job\tDone\t'$(date)>> $1'/'log.txt
echo 'END '$jobid