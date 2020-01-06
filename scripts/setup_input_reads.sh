#!/bin/bash
#################################################################################################
##
##  Script to setup read inputs to run an entire UDP cohort through the full WDL VG pipeline
##
##  Inputs:
##
##  Assumptions:
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up input directories and downloads files needed to run a cohort through the
VG WDL pipeline on the NIH Biowulf Cluster.

Inputs:
    -l List of individuals in cohort by UDP ID (in format UDP#### UDP#### UDP#### ...)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -c PATH to where the UDP cohort raw read files are located
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 2 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "l:w:c:t:h" OPTION; do
    case $OPTION in
        l)
            COHORT_NAMES_LIST+=($OPTARG)
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        c)
            INDIVIDUALS_DATA_DIR=$OPTARG
        ;;
        t)  
            RUN_SMALL_TEST=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

if [ ! -d "${COHORT_WORKFLOW_DIR}" ]; then
    mkdir -p ${COHORT_WORKFLOW_DIR}
    chmod 2770 ${COHORT_WORKFLOW_DIR}
fi

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
if [ ! -d "${READ_DATA_DIR}" ]; then
    mkdir -p ${READ_DATA_DIR}
    chmod 2770 ${READ_DATA_DIR}
fi
cd ${READ_DATA_DIR}

if [ $RUN_SMALL_TEST == false ]; then
    for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
    do
      INDIVIDUAL_DATA_DIR="${INDIVIDUALS_DATA_DIR}/${SAMPLE_NAME}"
      if [ $(find ${INDIVIDUAL_DATA_DIR}/ -name '*_R1.fastq.gz' | wc -l) -gt 0 ]; then
        ln -s $(find ${INDIVIDUAL_DATA_DIR}/ -name '*_R1.fastq.gz') ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_1.fq.gz
        ln -s $(find ${INDIVIDUAL_DATA_DIR}/ -name '*_R2.fastq.gz') ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_2.fq.gz
      else
        PAIR_1_READS=()
        PAIR_2_READS=()
        LANE_NUMS=($(find ${INDIVIDUAL_DATA_DIR}/ \( \( \( -wholename '*WGS*Rawreads*HudsonAlpha_1*fastq.gz' -o -wholename '*fastq.gz' \) -not -wholename '*Baylor*' \) -not -wholename '*WES*' \) | awk -F'-' '{print $2}'| awk -F'_' '{print $1"_"$2}' | sort | uniq | xargs))
        for LANE_NUM in ${LANE_NUMS[@]}
        do
          PAIR_1_READS+=("$(find ${INDIVIDUAL_DATA_DIR}/ \( \( \( -wholename '*WGS*Rawreads*HudsonAlpha_1*fastq.gz' -o -wholename '*fastq.gz' \) -not -wholename '*Baylor*' \) -not -wholename '*WES*' \) | grep "${LANE_NUM}_1")")
          PAIR_2_READS+=("$(find ${INDIVIDUAL_DATA_DIR}/ \( \( \( -wholename '*WGS*Rawreads*HudsonAlpha_1*fastq.gz' -o -wholename '*fastq.gz' \) -not -wholename '*Baylor*' \) -not -wholename '*WES*' \) | grep "${LANE_NUM}_2")")
        done
        echo "INDIVIDUAL_DATA_DIR: ${INDIVIDUAL_DATA_DIR}"
        echo "LANE_NUMS: ${LANE_NUMS}"
        echo "PAIR_1_READS: ${PAIR_1_READS[@]}"
        echo "PAIR_2_READS: ${PAIR_2_READS[@]}"
        cat ${PAIR_1_READS[@]} > ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_1.fq.gz
        cat ${PAIR_2_READS[@]} > ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_2.fq.gz
      fi
    done
else
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_chr21_1.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG002_read_pair_1.fq.gz
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_chr21_2.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG002_read_pair_2.fq.gz
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003_chr21_1.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG003_read_pair_1.fq.gz
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003_chr21_2.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG003_read_pair_2.fq.gz
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004_chr21_1.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG004_read_pair_1.fq.gz
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004_chr21_2.tiny.2x250.fastq.gz -O ${READ_DATA_DIR}/HG004_read_pair_2.fq.gz
fi

exit

