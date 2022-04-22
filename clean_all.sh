#!/bin/bash

#clean all directories to state 0

#clean genomes from variant genome and pending data
rm -rf Genomes/hg38+*
rm -rf Genomes/variants*

#clean results
rm -rf Results/*

#clean dictionaries
rm -rf Dictionaries/*

#remove genome_library
rm -rf genome_library/

