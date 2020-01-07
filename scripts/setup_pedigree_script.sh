#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the entire vg_wdl pedigree pipeline.
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

This script setups up a bash script to run a UDP cohort through all stages of the vg_wdl
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -c PATH to .ped file containing only the mother-father-proband trio samples
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 7 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "m:f:s:c:w:g:v:t:h" OPTION; do
    case $OPTION in
        m)
            MATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        f)
            PATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        s)
            SIBLING_SAMPLE_NAMES+=($OPTARG)
        ;;
        c)
            TRIO_PEDIGREE_FILE=$OPTARG
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

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "module load cromwell/40 python/3.6" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "source ${VG_WDL_DIR}/miniwdl_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
if [ $RUN_SMALL_TEST == false ]; then
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_trio_multi_map_call.wdl \\
MATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \\
MATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \\
PATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \\
PATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \\
${SIB_READ_PAIR_1_LIST} \\
${SIB_READ_PAIR_2_LIST} \\
${SIB_ID_LIST} \\
SAMPLE_NAME_MATERNAL='${MATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PATERNAL='${PATERNAL_SAMPLE_NAME}' \\
READS_PER_CHUNK=10000000 \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt' \\
XG_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg' \\
GCSA_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa' \\
GCSA_LCP_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa.lcp' \\
GBWT_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gbwt' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
REF_FASTA_GZ='${WORKFLOW_INPUT_DIR}/hs37d5.fa.gz' \\
PED_FILE='${TRIO_PEDIGREE_FILE}' \\
GEN_MAP_FILES='${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar' \\
GRAPH_NAME='${PROBAND_SAMPLE_NAME}_parental_graph_wgs' \\
USE_HAPLOTYPES='true' \\
MAKE_SNARLS='false' \\
USE_DECOYS='true' \\
DRAGEN_MODE='false' \\
SNPEFF_ANNOTATION='true' \\
SPLIT_READ_CORES=16 \\
SPLIT_READ_DISK=10 \\
MAP_CORES=16 \\
MAP_DISK=10 \\
MAP_MEM=80 \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_map_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
else
    echo "miniwdl cromwell ${VG_WDL_DIR}/vg_wdl/workflows/vg_trio_multi_map_call.wdl \\
MATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \\
MATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \\
PATERNAL_INPUT_READ_FILE_1='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz' \\
PATERNAL_INPUT_READ_FILE_2='${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz' \\
${SIB_READ_PAIR_1_LIST} \\
${SIB_READ_PAIR_2_LIST} \\
${SIB_ID_LIST} \\
SAMPLE_NAME_MATERNAL='${MATERNAL_SAMPLE_NAME}' \\
SAMPLE_NAME_PATERNAL='${PATERNAL_SAMPLE_NAME}' \\
READS_PER_CHUNK=100000 \\
PATH_LIST_FILE='${WORKFLOW_INPUT_DIR}/path_list_21.txt' \\
XG_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.xg' \\
GCSA_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa' \\
GCSA_LCP_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa.lcp' \\
GBWT_FILE='${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gbwt' \\
REF_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa' \\
REF_INDEX_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai' \\
REF_DICT_FILE='${WORKFLOW_INPUT_DIR}/hs37d5.dict' \\
VG_CONTAINER='quay.io/vgteam/vg:v1.19.0' \\
SNPEFF_DATABASE='${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip' \\
DRAGEN_REF_INDEX_NAME='hs37d5_v7' \\
UDPBINFO_PATH='Udpbinfo' \\
HELIX_USERNAME='${USER}' \\
REF_FASTA_GZ='${WORKFLOW_INPUT_DIR}/hs37d5.fa.gz' \\
PED_FILE='${TRIO_PEDIGREE_FILE}' \\
GEN_MAP_FILES='${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar' \\
GRAPH_NAME='${PROBAND_SAMPLE_NAME}_parental_graph_wgs' \\
USE_HAPLOTYPES='true' \\
MAKE_SNARLS='false' \\
USE_DECOYS='true' \\
DRAGEN_MODE='false' \\
SNPEFF_ANNOTATION='true' \\
SPLIT_READ_CORES=2 \\
SPLIT_READ_DISK=10 \\
MAP_CORES=4 \\
MAP_DISK=10 \\
MAP_MEM=10 \\
-c ${VG_WDL_DIR}/vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf \\
-d ${PROBAND_SAMPLE_NAME}_cohort_map_call.final_outputs" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
fi

exit

