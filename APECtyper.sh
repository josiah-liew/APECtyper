#!/bin/bash

########### APEC Typing Pipeline ###########

# Basic steps: TBD

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

OUTDIR=''    # name of output directory
INPUT=''     # name of input FASTA file(s)

function printHelp () {
	printf "Usage: APECtyper.sh [OPTIONS] -i [FASTA] -o [DIR]\n"
	printf "\t-h\t\tprint this message\n"
	printf "\t-v\t\tprint the version\n"
	printf "\t-i\t\tFASTA contig file(s)\n"
	printf "\t-o\t\toutput directory\n"
	printf "\t-c\t\tprint citation\n"
}

function checkDependencies () {
    if ! command -v $1
    then
        printf "Error: dependency $1 could not be found.\n" 
        exit 1
    fi
}

function mlstAnalysis () {
    echo "======== Running mlst ========"
    mkdir ${OUTDIR}/mlst
    mlst --scheme ecoli --csv $FASTA > ${OUTDIR}/mlst/mlst_results_${NAME}.csv
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

#### Parse command options and arguments ####
while getopts 'vhi:o:c' flag; do
  case "${flag}" in
    v) echo "$VERSION"
       exit 0 ;;
    h) printHelp
       exit 0 ;;
    i) INPUT=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    c) echo -e "\n$CITATION\n"
       exit 0 ;;
  esac
done

#### If no variables, print help ####
[[ $# == 0 ]] && { printHelp ; exit 1; }

echo $OUTDIR
echo $INPUT

#### Check for empty input variables ####
[[ -z "$INPUT" ]] && { echo "Error: Missing input contig file(s)." ; printHelp ; exit 1; }
[[ -z "$OUTDIR" ]] && { echo "Error: Missing a specified output directory." ; printHelp ; exit 1; }

#### Check that input file exists ####
[[ ! -f "$INPUT" ]] && { echo "Error: Input file(s) does not exist." ; exit 1; }

#### Check for dependencies ####
checkDependencies mlst
checkDependencies blastn

#### Check that output directory exists, create if does not exist ####
[[ ! -d "$OUTDIR" ]] && { mkdir "$OUTDIR" ; }

#### Generate list of input FASTA files ####
if [[ $INPUT =~ .*\*.* ]]; then
    ls -1 $INPUT > ${OUTDIR}/contigFiles.tmp
else
    echo $INPUT > ${OUTDIR}/contigFiles.tmp
fi

#### MLST and BLAST of each input fasta file ####
for FASTA in $(cat ${OUTDIR}/contigFiles.tmp); do
    
    FILE=${FASTA##*/}
    NAME=${FILE%.*}
    
    echo $FASTA
    echo $FILE
    echo $NAME
    
    echo "============== Analysis of ${NAME} =================="
    
    ##### Step 1: MLST #####
    mlstAnalysis
        # if non-zero exit status, print error and exit
        [[ $? -ne 0 ]] && { echo "Error when running mlst." ; rm -rf ${OUTDIR} ; exit 1; }
    
    ##### Step 2: BLAST ##### 
    blastAnalysis
        # if non-zero exit status, print error and exit
        [[ $? -ne 0 ]] && { echo "Error when running BLAST." ; rm -rf ${OUTDIR} ; exit 1; }

    echo "============== Analysis of ${NAME} Complete ==================" 
done
 
    
#### Remove temp files ####
rm -f ${OUTDIR}/*.tmp

echo "============== End =================="
exit 0
