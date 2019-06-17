#!/bin/bash
#################################################################################################
##
##  Script to setup the vg wdl python environment and download workflow inputs.
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

This script setups up a vg wdl python virtual environment for easier use of vg wdl and downloads
workflow inputs on the NIH Biowulf Cluster.

Inputs:
    -g PATH to the workflow input directory
    -v PATH to the vg wdl repository

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
while getopts "g:v:h" OPTION; do
    case $OPTION in
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

module load cromwell/40 git python/3.6

## Setup vg wdl python virtualenvironment
mkdir -p ${VG_WDL_DIR} && cd ${VG_WDL_DIR}
git clone https://github.com/cmarkello/miniwdl.git
virtualenv miniwdl_venv
source miniwdl_venv/bin/activate
pip3 install ./miniwdl
deactivate

## Setup and download workflow inputs
mkdir -p ${WORKFLOW_INPUT_DIR}
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_whole_genome.txt -O ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_decoys_wgs_t289.xg -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_decoys_wgs_t289.gcsa -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_decoys_wgs_t289.gcsa.lcp -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa.lcp
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.fai -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.dict -O ${WORKFLOW_INPUT_DIR}/hs37d5.dict
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.gz -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snpEff_v4_3_GRCh37.75.zip -O ${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/genetic_map_GRCh37.tar -O ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar

exit

