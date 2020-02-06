version 1.0

### vg_2nd_iter_siblings_multi_map.wdl ###
# Author: Charles Markello
# Description: Mapping workflow for proband-siblings against the parental graph reference.
#              Designed as the 4th step in a pedigree-backed graph alignment pipeline.

import "./vg_multi_map.wdl" as vgMultiMapWorkflow

workflow vgTrioPipeline {
    input {
        Array[File]+ SIBLING_INPUT_READ_FILE_1_LIST         # List of sibling 1st read pair .fastq.gz files. Proband must be first in this list.
        Array[File]+ SIBLING_INPUT_READ_FILE_2_LIST         # List of sibling 2nd read pair .fastq.gz files. Proband must be first in this list. Must follow same sample order as in SIBLING_INPUT_READ_FILE_1_LIST.
        Array[String]+ SAMPLE_NAME_SIBLING_LIST             # List of sibling sample names. Proband must be first in this list. Must follow same sample order as in SIBLING_INPUT_READ_FILE_1_LIST and SIBLING_INPUT_READ_FILE_2_LIST.
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 20000000                      # Number of reads contained in each mapping chunk (20000000 for wgs)
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File GCSA_FILE                                      # Path to .gcsa index file
        File GCSA_LCP_FILE                                  # Path to .gcsa.lcp index file
        File? GBWT_FILE                                     # (OPTIONAL) Path to .gbwt index file
        File? SNARLS_FILE                                   # (OPTIONAL) Path to .snarls index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        Int SPLIT_READ_CORES = 16
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 16
        Int MAP_DISK = 10
        Int MAP_MEM = 60
        Int MERGE_GAM_CORES = 56
        Int MERGE_GAM_DISK = 400
        Int MERGE_GAM_MEM = 100
        Int MERGE_GAM_TIME = 1200
        Boolean VGMPMAP_MODE = false                         # Set to 'false' to use "VG MAP" or set to 'true' to use "VG MPMAP" algorithm.
    }
   
    ##############################################################################
    ## Run mapping workflows on proband and siblings against the parental graph ##
    ##############################################################################
    Array[Pair[File,File]] read_pair_files_list = zip(SIBLING_INPUT_READ_FILE_1_LIST, SIBLING_INPUT_READ_FILE_2_LIST)
    scatter (read_pair_set in zip(read_pair_files_list, SAMPLE_NAME_SIBLING_LIST)) {
        Pair[File,File] read_pair_files = read_pair_set.left
        call vgMultiMapWorkflow.vgMultiMapCall {
            input:
                INPUT_READ_FILE_1=read_pair_files.left,
                INPUT_READ_FILE_2=read_pair_files.right,
                SAMPLE_NAME=read_pair_set.right,
                VG_CONTAINER=VG_CONTAINER,
                READS_PER_CHUNK=READS_PER_CHUNK,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=XG_FILE,
                GCSA_FILE=GCSA_FILE,
                GCSA_LCP_FILE=GCSA_LCP_FILE,
                GBWT_FILE=GBWT_FILE,
                SNARLS_FILE=SNARLS_FILE,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                SPLIT_READ_CORES=SPLIT_READ_CORES,
                SPLIT_READ_DISK=SPLIT_READ_DISK,
                MAP_CORES=MAP_CORES,
                MAP_DISK=MAP_DISK,
                MAP_MEM=MAP_MEM,
                MERGE_GAM_CORES=MERGE_GAM_CORES,
                MERGE_GAM_DISK=MERGE_GAM_DISK,
                MERGE_GAM_MEM=MERGE_GAM_MEM,
                MERGE_GAM_TIME=MERGE_GAM_TIME,
                VGMPMAP_MODE=VGMPMAP_MODE,
                SURJECT_MODE=true
        }
    }
    Array[File?] output_sibling_bam_list_maybes = vgMultiMapCall.output_bam
    Array[File?] output_sibling_bam_index_list_maybes = vgMultiMapCall.output_bam_index
    output {
        Array[File] output_sibling_bam_list = select_all(output_sibling_bam_list_maybes)
        Array[File] output_sibling_bam_index_list = select_all(output_sibling_bam_index_list_maybes) 
    }
}

