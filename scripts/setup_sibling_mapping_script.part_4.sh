#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the 4th stage of the vg_wdl pedigree pipeline.
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

This script setups up a bash script to run a UDP cohort through the 4th stage of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 4 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "s:w:g:v:t:h" OPTION; do
    case $OPTION in
        s)
            SIBLING_SAMPLE_NAMES+=($OPTARG)
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

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"

PROBAND_SAMPLE_NAME="${SIBLING_SAMPLE_NAMES[0]}"
READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
SIB_READ_PAIR_1_LIST=""
SIB_READ_PAIR_2_LIST=""
SIB_ID_LIST=""
for SIBLING_ID in ${SIBLING_SAMPLE_NAMES[@]}
do
  SIB_READ_PAIR_1_LIST+="SIBLING_INPUT_READ_FILE_1_LIST='${READ_DATA_DIR}/${SIBLING_ID}_read_pair_1.fq.gz' "
  SIB_READ_PAIR_2_LIST+="SIBLING_INPUT_READ_FILE_2_LIST='${READ_DATA_DIR}/${SIBLING_ID}_read_pair_2.fq.gz' "
  SIB_ID_LIST+="SAMPLE_NAME_SIBLING_LIST='${SIBLING_ID}' "
done
XG_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/outputs -regex .*.xg))
GCSA_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/outputs -regex .*.gcsa))
GCSA_LCP_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/outputs -regex .*.gcsa.lcp))
GBWT_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/outputs -regex .*.gbwt))
SNARLS_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs/outputs -regex .*.snarls))

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
echo "module load cromwell/40 python/3.6" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
echo "source ${VG_WDL_DIR}/miniwdl_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
if [ $RUN_SMALL_TEST == false ]; then
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_2nd_iter_siblings_multi_map.wdl \
        ${SIB_READ_PAIR_1_LIST} \
        ${SIB_READ_PAIR_2_LIST} \
        ${SIB_ID_LIST} \
        READS_PER_CHUNK=10000000 \
        PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt' \
        XG_FILE='${XG_FILE_PATH}' \
        GCSA_FILE='${GCSA_FILE_PATH}' \
        GCSA_LCP_FILE='${GCSA_LCP_FILE_PATH}' \
        GBWT_FILE='${GBWT_FILE_PATH}' \
        SNARLS_FILE='${SNARLS_FILE_PATH}' \
        REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \
        REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \
        REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \
        VG_CONTAINER='quay.io/vgteam/vg:v1.16.0' \
        SPLIT_READ_CORES=16 \
        SPLIT_READ_DISK=50 \
        MAP_CORES=32 \
        MAP_DISK=100 \
        MAP_MEM=100 \
        -c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \
        -d ${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
else
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_2nd_iter_siblings_multi_map.wdl \
        ${SIB_READ_PAIR_1_LIST} \
        ${SIB_READ_PAIR_2_LIST} \
        ${SIB_ID_LIST} \
        READS_PER_CHUNK=100000 \
        PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_21.txt' \
        XG_FILE='${XG_FILE_PATH}' \
        GCSA_FILE='${GCSA_FILE_PATH}' \
        GCSA_LCP_FILE='${GCSA_LCP_FILE_PATH}' \
        GBWT_FILE='${GBWT_FILE_PATH}' \
        SNARLS_FILE='${SNARLS_FILE_PATH}' \
        REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \
        REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \
        REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \
        VG_CONTAINER='quay.io/vgteam/vg:v1.16.0' \
        SPLIT_READ_CORES=2 \
        SPLIT_READ_DISK=10 \
        MAP_CORES=4 \
        MAP_DISK=10 \
        MAP_MEM=10 \
        -c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \
        -d ${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.part_4.sh
fi

exit

