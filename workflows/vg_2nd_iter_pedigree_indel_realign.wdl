version 1.0

### vg_2nd_iter_pedigree_multi_call.wdl ###
# Author: Charles Markello
# Description: Variant calling workflow for entire proband and sibling cohort.
#              Designed as the 5th step in a pedigree-backed graph alignment pipeline.

import "./vg_indel_realign.wdl" as vgIndelRealignmentWorkflow

###########################
### WORKFLOW DEFINITION ###
###########################
workflow vgTrioPipeline {
    input {
        File MATERNAL_INPUT_BAM_FILE                        # Input maternal surjected .bam file
        File MATERNAL_INPUT_BAM_FILE_INDEX                  # Input maternal .bai index of surjected .bam file
        File PATERNAL_INPUT_BAM_FILE                        # Input paternal surjected .bam file
        File PATERNAL_INPUT_BAM_FILE_INDEX                  # Input paternal .bai index of surjected .bam file
        Array[File]+ SIBLING_BAM_FILE_LIST                  # Input list of sibling surjected .bam files. Proband must be first in list.
        Array[File]+ SIBLING_BAM_FILE_INDEX_LIST            # Input list of .bai indices of surjected .bam files. Must follow same sample order as SIBLING_BAM_FILE_LIST.
        String SAMPLE_NAME_MATERNAL                         # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                         # Sample name for the father
        Array[String]+ SAMPLE_NAME_SIBLING_LIST             # Input list of sibling sample names. Must follow same order as SIBLING_BAM_FILE_LIST.
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.64.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.64.0)
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
    }
    
    
    ######################################################
    ## Run indel realignment workflows on pedigree bams ##
    ######################################################
    call vgIndelRealignmentWorkflow.vgIndelRealign as maternalIndelRealignmentWorkflow {
        input:
            INPUT_BAM_FILE=MATERNAL_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=MATERNAL_INPUT_BAM_FILE_INDEX,
            SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM
    }
    call vgIndelRealignmentWorkflow.vgIndelRealign as paternalIndelRealignmentWorkflow {
        input:
            INPUT_BAM_FILE=PATERNAL_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=PATERNAL_INPUT_BAM_FILE_INDEX,
            SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM
    }
    
    Array[Pair[File,File]] bam_pair_files_list = zip(SIBLING_BAM_FILE_LIST, SIBLING_BAM_FILE_INDEX_LIST)
    scatter (bam_pair_set in zip(bam_pair_files_list, SAMPLE_NAME_SIBLING_LIST)) {
        Pair[File,File] bam_pair_files = bam_pair_set.left
        call vgIndelRealignmentWorkflow.vgIndelRealign {
            input:
                INPUT_BAM_FILE=bam_pair_files.left,
                INPUT_BAM_FILE_INDEX=bam_pair_files.right,
                SAMPLE_NAME=bam_pair_set.right,
                VG_CONTAINER=VG_CONTAINER,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=XG_FILE,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                MAP_CORES=MAP_CORES,
                MAP_DISK=MAP_DISK,
                MAP_MEM=MAP_MEM
        }
    }
    Array[File?] output_sibling_bam_list_maybes = vgIndelRealign.output_bam
    Array[File?] output_sibling_bam_index_list_maybes = vgIndelRealign.output_bam_index
    
    output {
        File output_maternal_bam = maternalIndelRealignmentWorkflow.output_bam
        File output_maternal_bam_index = maternalIndelRealignmentWorkflow.output_bam_index
        File output_paternal_bam = paternalIndelRealignmentWorkflow.output_bam
        File output_paternal_bam_index = paternalIndelRealignmentWorkflow.output_bam_index
        Array[File] output_sibling_bam_list = select_all(output_sibling_bam_list_maybes)
        Array[File] output_sibling_bam_index_list = select_all(output_sibling_bam_index_list_maybes)
    }
}

