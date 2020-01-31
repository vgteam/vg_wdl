#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the 2nd stage of the vg_wdl pedigree pipeline.
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

This script setups up a bash script to run a UDP cohort through the 2nd stage of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -p Proband UDP ID (in format UDP####)
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 6 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "p:m:f:w:g:v:t:h" OPTION; do
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

MATERNAL_BAM_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${MATERNAL_SAMPLE_NAME}_merged.positionsorted.bam))
MATERNAL_BAM_BAI_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${MATERNAL_SAMPLE_NAME}_merged.positionsorted.bam.bai))
PATERNAL_BAM_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${PATERNAL_SAMPLE_NAME}_merged.positionsorted.bam))
PATERNAL_BAM_BAI_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${PATERNAL_SAMPLE_NAME}_merged.positionsorted.bam.bai))
PROBAND_BAM_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${PROBAND_SAMPLE_NAME}_merged.positionsorted.bam))
PROBAND_BAM_BAI_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs/output_links -name ${PROBAND_SAMPLE_NAME}_merged.positionsorted.bam.bai))

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
echo "module load cromwell/40 python/3.6" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
echo "source ${VG_WDL_DIR}/miniwdl_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
if [ $RUN_SMALL_TEST == false ]; then
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_trio_multi_call.wdl \\
MATERNAL_INPUT_BAM_FILE='${MATERNAL_BAM_PATH}' \\
MATERNAL_INPUT_BAM_FILE_INDEX='${MATERNAL_BAM_BAI_PATH}' \\
PATERNAL_INPUT_BAM_FILE='${PATERNAL_BAM_PATH}' \\
PATERNAL_INPUT_BAM_FILE_INDEX='${PATERNAL_BAM_BAI_PATH}' \\
PROBAND_INPUT_BAM_FILE='${PROBAND_BAM_PATH}' \\
PROBAND_INPUT_BAM_FILE_INDEX='${PROBAND_BAM_BAI_PATH}' \\
SAMPLE_NAME_MATERNAL='${MATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PATERNAL='${PATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PROBAND='${PROBAND_SAMPLE_NAME}' \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt' \\
XG_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
DRAGEN_MODE='true' \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
else
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_trio_multi_call.wdl \\
MATERNAL_INPUT_BAM_FILE='${MATERNAL_BAM_PATH}' \\
MATERNAL_INPUT_BAM_FILE_INDEX='${MATERNAL_BAM_BAI_PATH}' \\
PATERNAL_INPUT_BAM_FILE='${PATERNAL_BAM_PATH}' \\
PATERNAL_INPUT_BAM_FILE_INDEX='${PATERNAL_BAM_BAI_PATH}' \\
PROBAND_INPUT_BAM_FILE='${PROBAND_BAM_PATH}' \\
PROBAND_INPUT_BAM_FILE_INDEX='${PROBAND_BAM_BAI_PATH}' \\
SAMPLE_NAME_MATERNAL='${MATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PATERNAL='${PATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PROBAND='${PROBAND_SAMPLE_NAME}' \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_21.txt' \\
XG_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.xg' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
DRAGEN_MODE='true' \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.part_2.sh
fi

exit

