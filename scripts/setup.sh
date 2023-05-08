#!/bin/bash

#Creates needed directories for running snakemake pipeline
#Create folder struct below, cd to scripts folder an run setup.sh
#
#Screen_folder
#	|
#	|_scripts
#	|	|_setup.sh
#	|	|_other scritps
#	|
#	|_Snakefile

SCRIPT_PATH=${0%/*}

cd $SCRIPT_PATH
cd ..

mkdir dLFC_analyis
mkdir GEMINI_input
mkdir logs
mkdir metainfo
mkdir raw_sequencing_data
mkdir processed_sequencing_data
