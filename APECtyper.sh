#!/bin/bash

############################################
########### APEC Typing Pipeline ###########
############################################

# Basci steps: TBD

# Current version: 1.0 (Apr. 2022)
VERSION="APECtyper v.1.0 (Apr. 2022)"

# Required Software:
##
##
##



DIR=$( dirname "${BASH_SOURCE[0]}" )
echo $DIR


DATE=$(date +%Y-%m-%d)

# BLAST parameters
DB_APEC="${DIR}/db/apec_refs.fa"
PERC_IDENTITY=90
PERC_COVERAGE=90



function help () {
	printf "Usage: APECtyper.sh [OPTIONS] -i [FASTA or DIR] -o [DIR]\n"
	printf "\t-h\t\tprint this message\n"
	printf "\t-v\t\tprint the version\n"
	printf "\t-c\t\tcheck all dependencies in path\n"
	printf "\t-i\t\tFASTA contigs file or directory containing multiple FASTA files\n"
	printf "\t-o\t\toutput directory\n"
	printf "\t-r\t\tprint citation\n"
}

function checkDependencies () {
    if ! command -v $1
    then
        printf "ERROR: $1 could not be found." 
        exit 1
    fi
}

function mlstAnalysis () {
    echo "======== Running mlst ========"
    mlst --scheme ecoli --csv $FASTA  > ${OUTDIR}/mlst/mlst_results_${NAME}.csv
}


function blastAnalysis(){
   



}

function generateReport(){

}


# for each sample...

NAME=${x%.*}
