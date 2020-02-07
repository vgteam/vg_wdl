#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to extract outputs of the vg_wdl pedigree pipeline.
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

This script extracts file paths and copies output files of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -o PATH to the output directory
    -d (OPTIONAL, default=false) Set to 'true' if deleting leftover workflow raw output
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 5 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
DEL_RAW_OUTPUT=false

## Parse through arguments
while getopts "s:m:f:w:o:d:h" OPTION; do
    case $OPTION in
        s)
            SIBLING_SAMPLE_NAMES+=($OPTARG)
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
        o)
            OUTPUT_DIR=$OPTARG
        ;;
        d)
            DEL_RAW_OUTPUT=$OPTARG
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

cd ${COHORT_WORKFLOW_DIR}
if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p ${OUTPUT_DIR}
    chmod 2770 ${OUTPUT_DIR}
fi

PROBAND_SAMPLE_NAME="${SIBLING_SAMPLE_NAMES[0]}"
for SIBLING_ID in ${SIBLING_SAMPLE_NAMES[@]}
do
  SIB_BAM_FILE_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${SIBLING_ID}_merged.indel_realigned.bam))
  SIB_BAM_FILE_INDEX_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${SIBLING_ID}_merged.indel_realigned.bam.bai))
  SIB_GVCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.final_outputs -wholename *output_links/*${SIBLING_ID}_dragen_genotyped.hard-filtered.gvcf.gz))
  cp ${SIB_BAM_FILE_PATH} "${OUTPUT_DIR}/${SIBLING_ID}_parent_aligned.bam"
  cp ${SIB_BAM_FILE_INDEX_PATH} "${OUTPUT_DIR}/${SIBLING_ID}_parent_aligned.bam.bai"
  cp ${SIB_GVCF_PATH} "${OUTPUT_DIR}/${SIBLING_ID}_parent_aligned.gvcf.gz"
done

MATERNAL_BAM_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${MATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam))
MATERNAL_BAM_INDEX_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${MATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam.bai))
MATERNAL_GVCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs -wholename *output_links/*${MATERNAL_SAMPLE_NAME}_dragen_genotyped.hard-filtered.gvcf.gz))

cp ${MATERNAL_BAM_PATH} "${OUTPUT_DIR}/${MATERNAL_SAMPLE_NAME}_snp1kg_aligned.bam"
cp ${MATERNAL_BAM_INDEX_PATH} "${OUTPUT_DIR}/${MATERNAL_SAMPLE_NAME}_snp1kg_aligned.bam.bai"
cp ${MATERNAL_GVCF_PATH} "${OUTPUT_DIR}/${MATERNAL_SAMPLE_NAME}_snp1kg_aligned.gvcf.gz"

PATERNAL_BAM_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${PATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam))
PATERNAL_BAM_INDEX_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs -wholename *output_links/*${PATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam.bai))
PATERNAL_GVCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs -wholename *output_links/*${PATERNAL_SAMPLE_NAME}_dragen_genotyped.hard-filtered.gvcf.gz))

cp ${PATERNAL_BAM_PATH} "${OUTPUT_DIR}/${PATERNAL_SAMPLE_NAME}_snp1kg_aligned.bam"
cp ${PATERNAL_BAM_INDEX_PATH} "${OUTPUT_DIR}/${PATERNAL_SAMPLE_NAME}_snp1kg_aligned.bam.bai"
cp ${PATERNAL_GVCF_PATH} "${OUTPUT_DIR}/${PATERNAL_SAMPLE_NAME}_snp1kg_aligned.gvcf.gz"

JOINT_GENOTYPED_VCF_PATH=($(find ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.final_outputs/*_vgTrioPipeline/output_links/vgTrioPipeline.output_cohort_vcf -name ${PROBAND_SAMPLE_NAME}.snpeff.unrolled.vcf))

cp ${JOINT_GENOTYPED_VCF_PATH} "${OUTPUT_DIR}/"

tar -cvf ${PROBAND_SAMPLE_NAME}_workflow_outputs.tar ${OUTPUT_DIR}

## Delete intermediate workflow and input directories
if [ $DEL_RAW_OUTPUT == true ]; then
    if [[ $OUTPUT_DIR = *[[:space:]]* ]]; then 
        echo "ERROR: OUTPUT_DIR variable contains whitespace"
        exit 1
    fi
    if [[ $COHORT_WORKFLOW_DIR = *[[:space:]]* ]]; then 
        echo "ERROR: COHORT_WORKFLOW_DIR variable contains whitespace"
        exit 1
    fi
    if [[ $COHORT_WORKFLOW_DIR/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs = *[[:space:]]* ]]; then 
        echo "ERROR: PROBAND_SAMPLE_NAME variable contains whitespace"
        exit 1
    fi
    rm -fr ${OUTPUT_DIR}
    rm -fr ${COHORT_WORKFLOW_DIR}/input_reads
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_map.final_outputs
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_trio_call.final_outputs
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_parental_graph_construction.final_outputs
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_sibling_map.final_outputs
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_call.final_outputs
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_cohort_2nd_iter_pedigree_indel_realign.final_outputs
fi

exit

