#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the 1st stage of the vg_wdl pedigree pipeline.
##
##  Inputs:
##
##  Assumptions:
##      The UDP cohort ID is the same as the UDP ID for the proband sample.
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up a bash script to run a UDP cohort through the 1st stage of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -p Proband UDP ID (in format UDP####)
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository

Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 6 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## Parse through arguments
while getopts "p:m:f:w:g:v:h" OPTION; do
    case $OPTION in
        p)
            PROBAND_SAMPLE_NAME=$OPTARG
        ;;
        m)
            MATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        f)
            PATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        g)
            WORKFLOW_INPUT_DIR=$OPTARG
        ;;
        v)
            VG_WDL_DIR=$OPTARG
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

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh
echo "module load cromwell/40 python/3.6" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh
echo "source ${VG_WDL_DIR}/miniwdl_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh
echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_trio_multi_map.wdl \
    MATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \
    MATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \
    PATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \
    PATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \
    PROBAND_INPUT_READ_FILE_1='${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_1.fq.gz' \
    PROBAND_INPUT_READ_FILE_2='${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_2.fq.gz' \
    SAMPLE_NAME_MATERNAL='${MATERNAL_SAMPLE_NAME}' \
    SAMPLE_NAME_PATERNAL='${PATERNAL_SAMPLE_NAME}' \
    SAMPLE_NAME_PROBAND='${PROBAND_SAMPLE_NAME}' \
    READS_PER_CHUNK=10000000 \
    PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt' \
    XG_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg' \
    GCSA_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa' \
    GCSA_LCP_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa.lcp' \
    REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \
    REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \
    REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \
    VG_CONTAINER='quay.io/vgteam/vg:v1.16.0' \
    SPLIT_READ_CORES=16 \
    SPLIT_READ_DISK=100 \
    MAP_CORES=32 \
    MAP_DISK=100 \
    MAP_MEM=100 \
    -c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \
    -d ${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.part_1.sh

exit

