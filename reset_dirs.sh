#!/bin/bash


read -p "Do you want to remove all directories with download and user generated data? THE REMOVE IS IRREVERSIBILE (Yes/no)" answer

echo

if [[ $answer == "Yes" ]]; then

    echo "REMOVING ALL DIRECTORIES AND USER GENERATED DATA"
    rm -rf Genomes/
    rm -rf Results/
    rm -rf Dictionaries/
    rm -rf VCFs/
    rm -rf Annotations/
    rm -rf PAMs/
    rm -rf samplesIDs/

fi
