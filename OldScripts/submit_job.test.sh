#!/bin/sh
#$1 is directory of result for submitted job id (Results/job_id)
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
#Note that if genome_selected is not enriched, the python exe will force $15 as false

cd $1
rm queue.txt
jobid=$(basename $1)
echo 'Job\tStart\t'$(date)> log.txt
used_genome_dir=$2                 

#Start search index     #NOTE new version 2.1.2 of crispritz needed
echo 'Search-index\tStart\t'$(date) >> log.txt
echo 'Search_output '${19} >  output.txt
if [ ${10} = 'True' ]; then
    #echo 'crispritz search-index'
    crispritz.py search ../../$4 pam.txt guides.txt $jobid -mm $7 -bDNA $8 -bRNA ${9} -t #>> output.txt #TODO sistemare l'output redirection


    if [ ${15} = 'True' ]; then
        mkdir 'ref'
        echo 'Search_output_ref '${19} >>  output.txt
        crispritz.py search ../../${16} pam.txt guides.txt $jobid'_ref' -mm $7 -bDNA $8 -bRNA ${9} -t #>> output.txt #TODO sistemare l'output redirection
        mv ./$jobid'_ref'.*.txt 'ref'
        mv ./$jobid'_ref'.*.xls 'ref'
    fi
fi
echo 'Search-index\tDone\t'$(date) >> log.txt

#Start search
echo 'Search\tStart\t'$(date) >> log.txt
if [ ${11} = 'True' ]; then
    #echo 'crispritz search'
    crispritz.py search ../../$used_genome_dir pam.txt guides.txt $jobid -mm $7 -var -t #>> output.txt #-scores $3
    
    if [ ${15} = 'True' ]; then
        mkdir 'ref'
        echo 'Search_output_ref '${19} >>  output.txt 
        crispritz.py search ../../$3 pam.txt guides.txt $jobid'_ref' -mm $7 -var -t #>> output.txt
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
    python ../../PostProcess/cluster.dict.py $jobid.targets.txt 'addGuide' 'True' 'False' guides.txt # > $jobid.targets.cluster.txt
    echo 'End cluster ref'

    echo 'Start extraction top1 ref'
    python ../../PostProcess/extract_top.py $jobid.targets.cluster.txt $jobid  # > $jobid.top_1.txt
    echo 'End extraction top1 ref'

    echo 'Start Scoring'
    echo 'Calculating Scores... Step [2/3]' >>  output.txt
    python ../../PostProcess/scores_guide_table.py $jobid.top_1.txt ../../$used_genome_dir pam.txt guides.txt
    echo 'End Scoring'

    echo 'Start summary by pos-guide'
    echo 'Creating Summaries... Step [3/3]' >>  output.txt
    type_post='No'      
    python ../../PostProcess/summary_by_guide_position.py $jobid.targets.cluster.txt $7 $8 $9 guides.txt $jobid $type_post
    echo 'End summary by pos-guide'
    echo 'PostProcess\tDone\t'$(date) >> log.txt

    echo 'Annotation\tStart\t'$(date) >> log.txt
    echo 'Annotate_output '${19} >  output.txt
    echo 'Start annotation ref'
    if [ ${12} = 'True' ]; then
        crispritz.py annotate-results $jobid.top_1.txt ../../${18} $jobid >> output.txt # > $jobid.Annotation.summary.txt , $jobid.Annotation.targets.txt 
    fi
    echo 'End annotation ref'
    echo 'Annotation\tDone\t'$(date) >> log.txt

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
            echo $line >> output.txt
        done < guides.txt

        mkdir ../../assets/Img/$jobid
        cp *.png ../../assets/Img/$jobid/
        echo 'Report\tDone\t'$(date) >> log.txt
    fi

elif [ ${19} = 'var' ]; then
    echo 'PostProcess\tStart\t'$(date) >> log.txt
    echo 'Start cluster var'
    echo 'PostProcess_output' > output.txt
    echo 'Clustering... Step [1/6]' >>  output.txt
    python ../../PostProcess/cluster.dict.py $jobid.targets.txt 'no' 'True' 'False' guides.txt # > $jobid.targets.cluster.txt
    echo 'End cluster var'

    echo 'Start extraction top1 var'
    python ../../PostProcess/extract_top.py $jobid.targets.cluster.txt $jobid  # > $jobid.top_1.txt
    echo 'End extraction top1 var'

    echo 'Start Scoring'
    echo 'Calculating Scores... Step [2/6]' >>  output.txt
    python ../../PostProcess/scores_guide_table.py $jobid.top_1.txt ../../$used_genome_dir pam.txt guides.txt
    echo 'End Scoring'

    type_post='No'
    python ../../PostProcess/summary_by_guide_position.py $jobid.targets.cluster.txt $7 $8 $9 guides.txt $jobid $type_post

    echo 'Start pam analysis'
    echo 'PAM Analysis... Step [3/6]' >>  output.txt
    python ../../PostProcess/pam_analysis.py $jobid.targets.cluster.txt pam.txt ${19}
    echo 'End pam analysis'
    
    # Extract top 1
    echo 'Start extract top1'
    echo 'Extracting Samples... (This operation has a long execution time, Please Wait) Step [4/6]' >>  output.txt
    python ../../PostProcess/extract_top.py $jobid.targets.cluster.minmaxdisr.txt $jobid # > $jobid.top_1.txt
    echo 'End extract top1'
    # Expand top 1
    echo 'Start sort'
    sort -k4,4 $jobid.top_1.txt > $jobid.top_1.sort.txt && mv $jobid.top_1.sort.txt $jobid.top_1.txt 
    echo 'End sort'
    echo 'Start calc samples'
    python ../../PostProcess/calc_samples_faster.py ../../../dictionaries $jobid.top_1.txt   #> $jobid.top_1.samples.txt $jobid.top_1.samples.all.txt
    echo 'End calc samples'
    
    #Put right header into top_1.samples.all.txt
    sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tPam_disr\tSamples\tReal Guide' $jobid.top_1.samples.all.txt
    
    # Summary by samples table
    echo 'Start summary by samples'
    echo 'Creating Summaries... Step [5/6]' >>  output.txt
    python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt
    echo 'End summary by samples'

    #Rimettere i samples nel file di cluster (solo nel top1)
    echo 'Start creating final file'
    echo 'Preparing Files... Step [6/6]' >>  output.txt
    python ../../PostProcess/reassign_sample_to_cluster.py $jobid.targets.cluster.minmaxdisr.txt $jobid.top_1.samples.txt $jobid # > $jobid.final.txt
    echo 'End creating final file'
    echo 'PostProcess\tDone\t'$(date) >> log.txt

    #Annotation of top1 with samples
    echo 'Annotation\tStart\t'$(date) >> log.txt
    echo 'Annotate_output '${19} >  output.txt
    echo 'Start Annotation'
    crispritz.py annotate-results $jobid.top_1.samples.all.txt ../../${18} $jobid >> output.txt # > $jobid.Annotation.targets.txt $jobid.Annotation.summary.txt
                                                                                # $jobid.sample_annotation.guide.sample.txt $jobid..sumref.Annotation.summary.txt
    echo 'End Annotation'
    echo 'Annotation\tDone\t'$(date) >> log.txt

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
            echo $line >> output.txt
        done < guides.txt
        mkdir ../../assets/Img/$jobid
        cp *.png ../../assets/Img/$jobid/
        echo 'Report\tDone\t'$(date) >> log.txt
    fi

else    #Type search = both
    echo 'Report\tStart\t'$(date) >> log.txt
    #Estract common, semicommon and unique
    echo 'Start creation semicommon, common, unique'
    echo 'PostProcess_output' > output.txt
    echo 'Processing Search Results... Step [1/7]' >>  output.txt
    ../../PostProcess/./extraction.sh ref/$jobid'_ref.targets.txt' $jobid.targets.txt $jobid
    echo 'End creation semicommon, common, unique'

    #Cluster common file, extract top1 and insert into semicommon
    echo 'Start cluster common'
    echo 'Clustering... Step [2/7]' >>  output.txt
    python ../../PostProcess/cluster.dict.py $jobid.common_targets.txt 'no' 'False' 'False' guides.txt # > $jobid.common_targets.cluster.txt
                                                                        #Second false to not save the colum Pos Cluster and Total -> no needed for cat into semicommon
    echo 'End cluster common'
    echo 'Start top1 extraction common'
    python ../../PostProcess/extract_top.py $jobid.common_targets.cluster.txt $jobid.common_targets # > $jobid.common_targets.top_1.txt
    cat $jobid.common_targets.top_1.txt >> $jobid.semi_common_targets.txt
    echo 'End top1 extraction common'

    #Cluster semicommon e uniq
    echo 'Start cluster semicommon'
    python ../../PostProcess/cluster.dict.py $jobid.semi_common_targets.txt 'no' 'True' 'False' guides.txt
    echo 'End cluster semicommon'
    echo 'Start cluster unique'     #NOTE doing cluster separately does not create the right order of cluster (first clusters of uniq, then clusters of semi_common)
    python ../../PostProcess/cluster.dict.py $jobid.unique_targets.txt 'no' 'True' 'False' guides.txt
    echo 'End cluster unique'

    #Pam analysis
    echo 'Start pam analysis'
    echo 'PAM Analysis... Step [3/7]' >>  output.txt
    python ../../PostProcess/pam_analysis.py $jobid.semi_common_targets.cluster.txt pam.txt ${19}  # > $jobid.semi_common_targets.cluster.minmaxdisr.txt
    echo 'End pam analysis'
    echo 'Start pam creation'
    python ../../PostProcess/pam_creation.py $jobid.unique_targets.cluster.txt pam.txt ../../$3 # > $jobid.unique_targets.cluster.pamcreation.txt
    echo 'End pam creation'
    cat $jobid.unique_targets.cluster.pamcreation.txt $jobid.semi_common_targets.cluster.minmaxdisr.txt > $jobid.total.txt

    #Cluster of jobid.total.txt and extraction of top 1
    echo 'Start cluster of total.txt'
    echo 'Calculating Scores... Step [4/7]' >>  output.txt
    python ../../PostProcess/cluster.dict.py $jobid.total.txt 'no' 'True' 'True' 'total' guides.txt
    echo 'End cluster of total.txt'
    echo 'Start extract top1 total.txt'
    python ../../PostProcess/extract_top.py $jobid.total.cluster.txt $jobid # > $jobid.top_1.txt
    echo 'End extract top1 total.txt'

    #Scoring of top1
    echo 'Start Scoring'
    python ../../PostProcess/scores_guide_table.py $jobid.top_1.txt ../../$used_genome_dir pam.txt guides.txt
    echo 'End Scoring'

    #Summary guide, pos #NOTE the script automatically counts only for top subclusters
    echo 'Start summary by guide and position'  #NOTE change to top_1 if in sum by pos want to see cluster count of top1
    python ../../PostProcess/summary_by_guide_position.py $jobid.total.cluster.txt $7 $8 $9 guides.txt $jobid 'Uniq'
    echo 'End summary by guide and position'

    #Top1 expansion
    echo 'Start sort'
    echo 'Extracting Samples... (This operation has a long execution time, Please Wait) Step [5/7]' >>  output.txt
    sort -k4,4 $jobid.top_1.txt > $jobid.top_1.sort.txt && mv $jobid.top_1.sort.txt $jobid.top_1.txt 
    echo 'End sort'
    echo 'Start calc samples'
    python ../../PostProcess/calc_samples_faster.py ../../../dictionaries $jobid.top_1.txt  #> $jobid.top_1.samples.txt $jobid.top_1.samples.all.txt
    echo 'End calc samples'

    #Put right header into top_1.samples.all.txt
    sed -i '1 i\#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\tDirection\tMismatches\tBulge_Size\tTotal\tMin_mismatches\tMax_mismatches\tPam_disr\tPAM_gen\tVar_uniq\tSamples\tReal Guide' $jobid.top_1.samples.all.txt
    
    #Summary samples
    echo 'Start summary by samples'
    echo 'Creating Summaries... Step [6/7]' >>  output.txt
    python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt 
    echo 'End summary by samples'

    #Rimettere i samples nel file di cluster (solo nel top1)
    echo 'Start creating final file'
    python ../../PostProcess/reassign_sample_to_cluster.py $jobid.total.cluster.txt $jobid.top_1.samples.txt $jobid  # > $jobid.final.txt
    echo 'End creating final file'
    echo 'Preparing Files... Step [7/7]' >>  output.txt
    echo 'PostProcess\tDone\t'$(date) >> log.txt

    #Annotation of top1 with samples
    echo 'Annotation\tStart\t'$(date) >> log.txt
    echo 'Annotate_output '${19} >  output.txt
    echo 'Start Annotation'
    crispritz.py annotate-results $jobid.top_1.samples.all.txt ../../${18} $jobid >> output.txt # > $jobid.Annotation.targets.txt $jobid.Annotation.summary.txt
                                                                                # $jobid.sample_annotation.guide.sample.txt $jobid..sumref.Annotation.summary.txt
    echo 'End Annotation'
    echo 'Annotation\tDone\t'$(date) >> log.txt

    #Generate report
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
            echo $line >> output.txt
        done < guides.txt
        mkdir ../../assets/Img/$jobid
        cp *.png ../../assets/Img/$jobid/
    fi
    echo 'Report\tDone\t'$(date) >> log.txt

fi

cd ../../
if [ ${17} = 'True' ]; then
    python send_mail.py $1
fi
echo 'Job\tDone\t'$(date)>> $1'/'log.txt
