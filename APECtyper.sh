#!/bin/bash

############################################
########### APEC Typing Pipeline ###########
############################################

# Basci steps: TBD

# Current version: 1.0 (Apr. 2022)
VERSION="APECtyper v.1.0 (Apr. 2022)"

# Please cite:
CITATION="TBD"

# Required Software:
## NCBI BLAST+ blastn >= 2.9.0 (dependency of mlst)
## mlst (incl. dependencies)
##
##


DIR=$( dirname "${BASH_SOURCE[0]}" )
DATE=$(date +%Y-%m-%d)

# BLAST parameters
DB_APEC="${DIR}/db/apec_refs.fa"
PERC_IDENTITY=90
PERC_COVERAGE=90

#
OUTDIR=""
INPUT=""

function help () {
	printf "Usage: APECtyper.sh [OPTIONS] -i [FASTA or DIR] -o [DIR]\n"
	printf "\t-h\t\tprint this message\n"
	printf "\t-v\t\tprint the version\n"
	printf "\t-i\t\tFASTA contigs file or directory containing multiple FASTA files\n"
	printf "\t-o\t\toutput directory\n"
	printf "\t-c\t\tprint citation\n"
	exit 0
}

function checkDependencies () {
    if ! command -v $1
    then
        printf "Error: $1 could not be found." 
        exit 1
    fi
}

function mlstAnalysis () {
    echo "======== Running mlst ========"
    mlst --scheme ecoli --csv $FASTA  > ${OUTDIR}/mlst/mlst_results_${NAME}.csv
}

function makeBlastDB () {
    echo "======== Making BLAST database ========"
    rm -f $DIR/db/apec_refs.fa.*
    makeblastdb -in $DIR/db/apec_refs.fa -dbtype nucl -parse_seqids -title "APEC References"
}

function blastAnalysis () {
    echo "===== Running BLAST ====="
    blastn -query $FASTA -db $DIR/db/apec_refs.fa -perc_identity $PERC_IDENTITY -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore" -out ${OUTDIR}/blast/blast_results_${NAME}.tsv
}

function generateReport () {
   echo "This is where we will generate the report."
}


if [ $# == 0 ]
then
    help
    exit 1
fi

while getopts 'vhi:o:c' flag; do
  case "${flag}" in
    v) echo "$VERSION"
       exit 0 ;;
    h) help
       exit 0 ;;
    i) INPUT=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    c) echo -e "\n$CITATION\n" ;;
  esac
done

#### Check for empty input variables ####
[[ -z "$INPUT" ]] && { echo "Missing a contig file or directory." ; exit 1; }
[[ -z "$OUTPUT" ]] && { echo "Missing a specified output directory." ; exit 1; }

#### Check that input file/directory exists ####
if [ ! -f $INPUT ] && [ ! -d $INPUT ]
then
    printf "\nError: Input file or directory does not exist.\n"
    exit 1
fi

#### Check that output directory exists, create if does not exist ####
if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi

#### Generate list of input FASTA files ####
if [[ -f $INPUT ]]; then
    echo $INPUT > ${OUTDIR}/contigFiles.tmp
elif [[ -d $INPUT ]]; then
    ls -1 $INPUT > ${OUTDIR}/contigFiles.tmp
fi


#### MLST and BLAST of each input fasta file ####
for FASTA in $(cat ${OUTDIR}/contigFiles.tmp)
do
    NAME=${f%.*}
    echo "============== Analysis of ${NAME} =================="
    
    ##### Step 1: MLST #####
    mlstAnalysis
    if [ $? -eq 0 ]
    then
        printf "\nError when running mlst.\n"
        exit 1
    fi
    
    ##### Step 2: BLAST ##### 
    blastAnalysis
    if [ $? -eq 0 ]
    then
        printf "\nError when running BLAST.\n"
        exit 1
    fi
    
    echo "============== Analysis of ${NAME} Complete =================="
    
done
 
    
#### Remove temp files ####
rm -f ${OUTDIR}/contigFiles.tmp

echo "============== End =================="
exit 0
