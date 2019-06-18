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
    -i Cohort UDP ID (in format UDP####)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored

Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 2 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## Parse through arguments
while getopts "i:w:h" OPTION; do
    case $OPTION in
        i)
            COHORT_NAME=$OPTARG
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
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

## STEP1: COLLECT READS. Generate sample work directories, and setup input READS
READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
mkdir -p ${READ_DATA_DIR} && cd ${READ_DATA_DIR}

COHORT_NAMES_LIST=($(ls /data/Udpdata/CMarkello/${COHORT_NAME}/ | grep 'UDP' | uniq))
for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
do
  INDIVIDUAL_DATA_DIR="/data/Udpdata/CMarkello/${COHORT_NAME}/${SAMPLE_NAME}"
  PAIR_1_READS=()
  PAIR_2_READS=()
  LANE_NUMS=($(ls ${INDIVIDUAL_DATA_DIR} | awk -F'-' '{print $2}'| awk -F'_' '{print $1"_"$2}' | sort | uniq | xargs))
  for LANE_NUM in ${LANE_NUMS[@]}
  do
    PAIR_1_READS+=(${INDIVIDUAL_DATA_DIR}/"$(ls ${INDIVIDUAL_DATA_DIR} | grep "${LANE_NUM}_1")")
    PAIR_2_READS+=(${INDIVIDUAL_DATA_DIR}/"$(ls ${INDIVIDUAL_DATA_DIR} | grep "${LANE_NUM}_2")")
  done
  cat ${PAIR_1_READS[@]} > ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_1.fq.gz
  cat ${PAIR_2_READS[@]} > ${READ_DATA_DIR}/${SAMPLE_NAME}_read_pair_2.fq.gz
done

exit

