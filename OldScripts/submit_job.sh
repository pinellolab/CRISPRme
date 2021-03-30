#!/bin/sh
#$1 is directory of result for submitted job id (Results/job_id)
#$2 is genome_selected directory Eg Genomes/hg19_ref or Genomes/hg19_ref+1000genomeproject #NOTE + or - to decide
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
jobid=$(basename $1)
echo 'Job\tStart\t'$(date)> $1'/'log.txt
used_genome_dir=$2                 

#Start search index     #NOTE new version 2.1.2 of crispritz needed
echo 'Search-index\tStart\t'$(date) >> $1'/'log.txt
echo 'Search_output '${19} >  $1/output.txt
if [ ${10} = 'True' ]; then
    #echo 'crispritz search-index'
    crispritz.py search $4 $5 $6 $jobid -index -mm $7 -bDNA $8 -bRNA ${9} -t >> $1/output.txt #TODO sistemare l'output redirection
    mv ./$jobid.*.txt $1
    mv ./$jobid.*.xls $1

    if [ ${15} = 'True' ]; then
        mkdir $1'/ref'
        echo 'Search_output_ref '${19} >>  $1/output.txt
        crispritz.py search ${16} $5 $6 $jobid'_ref' -index -mm $7 -bDNA $8 -bRNA ${9} -t >> $1/output.txt #TODO sistemare l'output redirection
        mv ./$jobid'_ref'.*.txt $1/ref
        mv ./$jobid'_ref'.*.xls $1/ref
    fi
fi
echo 'Search-index\tDone\t'$(date) >> $1'/'log.txt

#Start search
echo 'Search\tStart\t'$(date) >> $1'/'log.txt
if [ ${11} = 'True' ]; then
    #echo 'crispritz search'
    crispritz.py search $used_genome_dir $5 $6 $jobid -mm $7 -t >> $1/output.txt #-scores $3
    mv ./$jobid.*.txt $1
    mv ./$jobid.*.xls $1
    if [ ${15} = 'True' ]; then
        mkdir $1'/ref'
        echo 'Search_output_ref '${19} >>  $1/output.txt 
        crispritz.py search $3 $5 $6 $jobid'_ref' -mm $7 -t >> $1/output.txt
        mv ./$jobid'_ref'.*.txt $1/ref
        mv ./$jobid'_ref'.*.xls $1/ref
    fi
fi
echo 'Search\tDone\t'$(date) >> $1'/'log.txt


#Start annotation       #NOTE new version 2.1.2 of crispritz needed
echo 'Annotation\tStart\t'$(date) >> $1'/'log.txt
echo 'Annotate_output '${19} >  $1/output.txt
if [ ${12} = 'True' ]; then
    #echo 'crispritz annotate'
    if [ ${10} = 'True' ]; then         #Indexed search was done, read profile complete
        crispritz.py annotate-results $1'/'/$jobid.profile_complete.xls $1'/'$jobid'.targets.txt' ${18} $jobid >> $1/output.txt
    fi
    if [ ${11} = 'True' ]; then         #Normal search was done, read profile
        crispritz.py annotate-results $1'/'/$jobid.profile.xls $1'/'$jobid'.targets.txt' ${18} $jobid >> $1/output.txt
    fi
    mv ./$jobid.Annotation*.txt $1
    if [ ${15} = 'True' ]; then
        echo 'Annotate_output_ref '${19} >>  $1/output.txt
        if [ ${10} = 'True' ]; then
            crispritz.py annotate-results $1'/ref/'$jobid'_ref.profile_complete.xls' $1'/ref/'$jobid'_ref.targets.txt' ${18} $jobid'_ref' >> $1/output.txt
        fi
        if [ ${11} = 'True' ]; then
            crispritz.py annotate-results $1'/ref/'$jobid'_ref.profile.xls' $1'/ref/'$jobid'_ref.targets.txt' ${18} $jobid'_ref' >> $1/output.txt
        fi
        mv ./$jobid'_ref'.Annotation*.txt $1/ref
    fi
fi
echo 'Annotation\tDone\t'$(date) >> $1'/'log.txt

#Start generate report
echo 'Report\tStart\t'$(date) >> $1'/'log.txt
#python summary_guide.py $1 $7  #New annotated file are alredy located in right directories
profile_type='profile_complete'
if [ ${11} = 'True' ]; then
    profile_type='profile'
fi

cd $1
if [ ${13} = 'True' ]; then
    echo 'Generate_report' >  output.txt
    proc=$(($7 + 1))
    while IFS= read -r line || [ -n "$line" ]; do    
        if [ ${14} = 'True' ]; then         #If -gecko
            if [ ${15} = 'True' ]; then     #If genome_ref comparison
                
                printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -profile $jobid'.'$profile_type'.xls' -extprofile *.extended_profile.xls -annotation $jobid'.Annotation.txt' -sumref ref/$jobid'_ref'.Annotation.txt  -gecko
            else
                printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -profile $jobid'.'$profile_type'.xls' -extprofile *.extended_profile.xls -annotation $jobid'.Annotation.txt' -gecko
            fi
        else
            if [ ${15} = 'True' ]; then     #If genome_ref comparison
                printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -profile $jobid'.'$profile_type'.xls' -extprofile *.extended_profile.xls -annotation $jobid'.Annotation.txt' -sumref ref/$jobid'_ref'.Annotation.txt
            else
                printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % crispritz.py generate-report $line -mm % -profile $jobid'.'$profile_type'.xls' -extprofile *.extended_profile.xls -annotation $jobid'.Annotation.txt'
            fi
        fi
        echo $line >> output.txt
    done < guides.txt
  
fi


cd ../../

mkdir assets/Img/$jobid
cp $PWD/$1/*.png assets/Img/$jobid/
echo 'Report\tDone\t'$(date) >> $1'/'log.txt

#Temp (find better way when scores will be modified), sum of cfd per guide 
#TODO scores_test will be substituted with -scores option with crispritz 2.1.2
echo 'PostProcess\tStart\t'$(date) >> $1'/'log.txt
cd $1
echo 'Post Process start'
echo 'Start Scoring'
python ../../PostProcess/scores_guide_table.py $jobid.targets.txt ../../$used_genome_dir pam.txt guides.txt
echo 'End Scoring'
#Analysis for var/ref type ('both')
if [ ${19} = 'both' ]; then     #TODO CHECK FOR LAST COL INDICES
    #Estract common, semicommon and unique
    echo 'Start creation semicommon, common, unique'
    ../../PostProcess/./extraction.sh ref/$jobid'_ref.targets.txt' $jobid.targets.txt $jobid
    echo 'End creation semicommon, common, unique'
    #Cluster semicommon e uniq -> TODO da sistemare l'ordine dell'analisi
    echo 'Start cluster semicommon'
    python ../../PostProcess/cluster.dict.py $jobid.semi_common_targets.txt
    echo 'End cluster semicommon'
    echo 'Start cluster unique'
    python ../../PostProcess/cluster.dict.py $jobid.unique_targets.txt
    echo 'End cluster unique'
    #Pam analysis
    echo 'Start pam analysis'
    python ../../PostProcess/pam_analysis.py $jobid.semi_common_targets.cluster.txt pam.txt ${19}  # > $jobid.semi_common_targets.cluster.minmaxdisr.txt
    echo 'End pam analysis'
    echo 'Start pam creation'
    python ../../PostProcess/pam_creation.py $jobid.unique_targets.cluster.txt pam.txt ../../$3 # > $jobid.unique_targets.cluster.pamcreation.txt
    echo 'End pam creation'
    cat $jobid.unique_targets.cluster.pamcreation.txt $jobid.semi_common_targets.cluster.minmaxdisr.txt > $jobid.total.txt
    #Summary guide, pos
    echo 'Start summary by guide and position'
    python ../../PostProcess/summary_by_guide_position.py $jobid.total.txt $7 $8 $9 guides.txt $jobid 'Uniq'
    echo 'End summary by guide and position'
    # mv $jobid.total.txt $jobid.targets.cluster.txt
    #Top 1 extraction
    echo 'Start top 1 extraction'
    python ../../PostProcess/extract_top.py $jobid.total.txt $jobid # > $jobid.top_1.txt
    echo 'End top 1 extraction'
    #Top1 expansion
    echo 'Start sort'
    sort -k4,4 $jobid.top_1.txt > $jobid.top_1.sort.txt && mv $jobid.top_1.sort.txt $jobid.top_1.txt 
    echo 'End sort'
    echo 'Start calc samples'
    python ../../PostProcess/calc_samples_faster.py ../../../dictionaries $jobid.top_1.txt  #> $jobid.top_1.samples.txt
    echo 'End calc samples'
    #Summary samples
    echo 'Start summary by samples'
    python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt 
    echo 'End summary by samples'
    #Annotazioni per samples, pop, superpop
    echo 'Start annotation samples'
    python ../../PostProcess/annotation_samples.py $jobid.top_1.samples.txt $jobid.Annotation.targets.txt $jobid.Annotation.txt $jobid
    echo 'End annotation samples'
    #Rimettere i samples nel file di cluster (solo nel top1)
    echo 'Start creating final file'
    python ../../PostProcess/reassign_sample_to_cluster.py $jobid.total.txt $jobid.top_1.samples.txt $jobid  # > $jobid.final.txt
    echo 'End creating final file'

    # #TODO sistemare fare script unico
    # python ../../PostProcess/cluster.dict.py ref/$jobid'_ref'.targets.txt #jobid_ref.targets.cluster.txt
    # python ../../PostProcess/extract_top.py ref/$jobid'_ref'.targets.cluster.txt $jobid'_ref' # > $jobid_ref.top_1.txt
    # python ../../file_per_crispritz/annotator.py ../../${18} $jobid'_ref'.top_1.txt $jobid'_ref' # $jobid'_ref'.Annotation.summary.txt
    # python ../../PostProcess/tmp_top1_annotation.py $jobid ./
    # python ../../PostProcess/tmp_top1_annotation.py $jobid'_ref' ./ref/
    # if [ ${13} = 'True' ]; then
    # echo 'Generate_report' >  output.txt
    # proc=$(($7 + 1))
    # while IFS= read -r line || [ -n "$line" ]; do    
    #     if [ ${14} = 'True' ]; then         #If -gecko
    #         if [ ${15} = 'True' ]; then     #If genome_ref comparison
                
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls $jobid'_ref'.Annotation.summary.txt  /home/ubuntu/miniconda3/opt/crispritz/Python_Scripts/Plot/gecko/

    #         else
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls /home/ubuntu/miniconda3/opt/crispritz/Python_Scripts/Plot/gecko/

    #         fi
    #     else
    #         if [ ${15} = 'True' ]; then     #If genome_ref comparison
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls $jobid'_ref'.Annotation.summary.txt no
    #         else
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls no  
    #         fi
    #     fi
    #     echo $line >> output.txt
    # done < guides.txt
    # fi
fi

#Clustering for var and ref
if [ ${19} = 'ref' ]; then
    echo 'Start cluster ref'
    python ../../PostProcess/cluster.dict.py $jobid.targets.txt 'addGuide'
    echo 'End cluster ref'
elif [ ${19} = 'var' ]; then
    echo 'Start cluster var'
    python ../../PostProcess/cluster.dict.py $jobid.targets.txt
    echo 'End cluster var'
fi


if [ ${19} = 'ref' ]; then
    type_post='No'
    python ../../PostProcess/summary_by_guide_position.py $jobid.targets.cluster.txt $7 $8 $9 guides.txt $jobid $type_post
elif [ ${19} = 'var' ]; then
    type_post='No'
    python ../../PostProcess/summary_by_guide_position.py $jobid.targets.cluster.txt $7 $8 $9 guides.txt $jobid $type_post

    echo 'Start pam analysis'
    python ../../PostProcess/pam_analysis.py $jobid.targets.cluster.txt pam.txt ${19}
    echo 'End pam analysis'
    # Extract top 1
    echo 'Start extract top1'
    python ../../PostProcess/extract_top.py $jobid.targets.cluster.minmaxdisr.txt $jobid # > $jobid.top_1.txt
    echo 'End extract top1'
    # Expand top 1
    echo 'Start sort'
    sort -k4,4 $jobid.top_1.txt > $jobid.top_1.sort.txt && mv $jobid.top_1.sort.txt $jobid.top_1.txt 
    echo 'End sort'
    echo 'Start calc samples'
    python ../../PostProcess/calc_samples_faster.py ../../../dictionaries $jobid.top_1.txt   #> $jobid.top_1.samples.txt
    echo 'End calc samples'
    # Summary by samples table
    echo 'Start summary by samples'
    python ../../PostProcess/summary_by_samples.py $jobid.top_1.samples.txt $jobid ${19} guides.txt
    echo 'End summary by samples'
    echo 'Start annotation samples'
    python ../../PostProcess/annotation_samples.py $jobid.top_1.samples.txt $jobid.Annotation.targets.txt $jobid.Annotation.txt $jobid
    echo 'End annotation samples'
    #Rimettere i samples nel file di cluster (solo nel top1)
    echo 'Start creating final file'
    python ../../PostProcess/reassign_sample_to_cluster.py $jobid.targets.cluster.minmaxdisr.txt $jobid.top_1.samples.txt $jobid # > $jobid.final.txt
    echo 'End creating final file'

    # #TODO sistemare fare script unico
    
    # python ../../PostProcess/tmp_top1_annotation.py $jobid ./
    # python ../../PostProcess/tmp_top1_annotation.py $jobid'_ref' ./ref/
    # if [ ${13} = 'True' ]; then
    # echo 'Generate_report' >  output.txt
    # proc=$(($7 + 1))
    # while IFS= read -r line || [ -n "$line" ]; do    
    #     if [ ${14} = 'True' ]; then         #If -gecko
    #         if [ ${15} = 'True' ]; then     #If genome_ref comparison
                
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls ref/$jobid'_ref.tmp_res.txt'  /home/ubuntu/miniconda3/opt/crispritz/Python_Scripts/Plot/gecko/

    #         else
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls /home/ubuntu/miniconda3/opt/crispritz/Python_Scripts/Plot/gecko/

    #         fi
    #     else
    #         if [ ${15} = 'True' ]; then     #If genome_ref comparison
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls ref/$jobid'_ref.tmp_res.txt' no
    #         else
    #             printf %s\\n $(seq 0 $7) | xargs -n 1 -P $proc -I % python ../../PostProcess/radar_chart_new.py  $line % $jobid'.tmp_res.txt' *.extended_profile.xls no  
    #         fi
    #     fi
    #     echo $line >> output.txt
    # done < guides.txt
  
    # fi

    
    #TODO aggiungere terza/quarta voce nella pagina del load
fi

echo 'Post Process done'
cd ../../
echo 'PostProcess\tDone\t'$(date) >> $1'/'log.txt
if [ ${17} = 'True' ]; then
    python send_mail.py $1
fi
echo 'Job\tDone\t'$(date)>> $1'/'log.txt
