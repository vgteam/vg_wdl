version 1.0

### vg_2nd_iter_pedigree_multi_call.wdl ###
# Author: Charles Markello
# Description: Variant calling workflow for entire proband and sibling cohort.
#              Designed as the 5th step in a pedigree-backed graph alignment pipeline.

import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "./vg_multi_call.wdl" as vgMultiCallWorkflow

###########################
### WORKFLOW DEFINITION ###
###########################
workflow vgTrioPipeline {
    input {
        File MATERNAL_GVCF                                  # Input maternal .gvcf.gz file
        File PATERNAL_GVCF                                  # Input paternal .gvcf.gz file
        Array[File]+ SIBLING_BAM_FILE_LIST                  # Input list of sibling surjected .bam files. Proband must be first in list.
        Array[File]+ SIBLING_BAM_FILE_INDEX_LIST            # Input list of .bai indices of surjected .bam files. Must follow same sample order as SIBLING_BAM_FILE_LIST.
        Array[String]+ SAMPLE_NAME_SIBLING_LIST             # Input list of sibling sample names. Must follow same order as SIBLING_BAM_FILE_LIST.
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int CHUNK_BASES = 50000000                          # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                                  # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        File SNPEFF_DATABASE                                # Path to snpeff database .zip file
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 400
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
        String DRAGEN_REF_INDEX_NAME                        # Dragen module based reference index directory (e.g. "hs37d5_v7")
        String UDPBINFO_PATH                                # Udp data directory to use for Dragen module (e.g. "Udpbinfo", nih biowulf system only)
        String HELIX_USERNAME                               # The nih helix username which holds a user directory in UDPBINFO_PATH
        Boolean DRAGEN_MODE                                 # Set to 'true' to use the Dragen modules variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
        Boolean SNPEFF_ANNOTATION = false                   # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
    }
    
    
    ###########################################
    ## Run variant calling workflows on Trio ##
    ###########################################
    Int numSilbings = length(SAMPLE_NAME_SIBLING_LIST)
    call vgMultiCallWorkflow.vgMultiMapCall as probandCallWorkflow {
        input:
            INPUT_BAM_FILE=SIBLING_BAM_FILE_LIST[0],
            INPUT_BAM_FILE_INDEX=SIBLING_BAM_FILE_INDEX_LIST[0],
            SAMPLE_NAME=SAMPLE_NAME_SIBLING_LIST[0],
            VG_CONTAINER=VG_CONTAINER,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            SURJECT_MODE=true,
            DRAGEN_MODE=DRAGEN_MODE,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false
    }
    if (numSilbings > 1) {
        call vgMultiCallWorkflow.vgMultiMapCall as siblingCallWorkflow1 {
            input:
                INPUT_BAM_FILE=SIBLING_BAM_FILE_LIST[1],
                INPUT_BAM_FILE_INDEX=SIBLING_BAM_FILE_INDEX_LIST[1],
                SAMPLE_NAME=SAMPLE_NAME_SIBLING_LIST[1],
                VG_CONTAINER=VG_CONTAINER,
                CHUNK_BASES=CHUNK_BASES,
                OVERLAP=OVERLAP,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=XG_FILE,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                SNPEFF_DATABASE=SNPEFF_DATABASE,
                CHUNK_GAM_CORES=CHUNK_GAM_CORES,
                CHUNK_GAM_DISK=CHUNK_GAM_DISK,
                CHUNK_GAM_MEM=CHUNK_GAM_MEM,
                VGCALL_CORES=VGCALL_CORES,
                VGCALL_DISK=VGCALL_DISK,
                VGCALL_MEM=VGCALL_MEM,
                DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
                UDPBINFO_PATH=UDPBINFO_PATH,
                HELIX_USERNAME=HELIX_USERNAME,
                SURJECT_MODE=true,
                DRAGEN_MODE=DRAGEN_MODE,
                GVCF_MODE=true,
                SNPEFF_ANNOTATION=false,
                PREVIOUS_WORKFLOW_OUTPUT= if DRAGEN_MODE then probandCallWorkflow.output_vcf else "null"
        }
    }
    if (numSilbings > 2) {
        call vgMultiCallWorkflow.vgMultiMapCall as siblingCallWorkflow2 {
            input:
                INPUT_BAM_FILE=SIBLING_BAM_FILE_LIST[2],
                INPUT_BAM_FILE_INDEX=SIBLING_BAM_FILE_INDEX_LIST[2],
                SAMPLE_NAME=SAMPLE_NAME_SIBLING_LIST[2],
                VG_CONTAINER=VG_CONTAINER,
                CHUNK_BASES=CHUNK_BASES,
                OVERLAP=OVERLAP,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=XG_FILE,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                SNPEFF_DATABASE=SNPEFF_DATABASE,
                CHUNK_GAM_CORES=CHUNK_GAM_CORES,
                CHUNK_GAM_DISK=CHUNK_GAM_DISK,
                CHUNK_GAM_MEM=CHUNK_GAM_MEM,
                VGCALL_CORES=VGCALL_CORES,
                VGCALL_DISK=VGCALL_DISK,
                VGCALL_MEM=VGCALL_MEM,
                DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
                UDPBINFO_PATH=UDPBINFO_PATH,
                HELIX_USERNAME=HELIX_USERNAME,
                SURJECT_MODE=true,
                DRAGEN_MODE=DRAGEN_MODE,
                GVCF_MODE=true,
                SNPEFF_ANNOTATION=false,
                PREVIOUS_WORKFLOW_OUTPUT= if DRAGEN_MODE then siblingCallWorkflow1.output_vcf else "null"
        }
    }
    if (numSilbings > 3) {
        call vgMultiCallWorkflow.vgMultiMapCall as siblingCallWorkflow3 {
            input:
                INPUT_BAM_FILE=SIBLING_BAM_FILE_LIST[3],
                INPUT_BAM_FILE_INDEX=SIBLING_BAM_FILE_INDEX_LIST[3],
                SAMPLE_NAME=SAMPLE_NAME_SIBLING_LIST[3],
                VG_CONTAINER=VG_CONTAINER,
                CHUNK_BASES=CHUNK_BASES,
                OVERLAP=OVERLAP,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=XG_FILE,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                SNPEFF_DATABASE=SNPEFF_DATABASE,
                CHUNK_GAM_CORES=CHUNK_GAM_CORES,
                CHUNK_GAM_DISK=CHUNK_GAM_DISK,
                CHUNK_GAM_MEM=CHUNK_GAM_MEM,
                VGCALL_CORES=VGCALL_CORES,
                VGCALL_DISK=VGCALL_DISK,
                VGCALL_MEM=VGCALL_MEM,
                DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
                UDPBINFO_PATH=UDPBINFO_PATH,
                HELIX_USERNAME=HELIX_USERNAME,
                SURJECT_MODE=true,
                DRAGEN_MODE=DRAGEN_MODE,
                GVCF_MODE=true,
                SNPEFF_ANNOTATION=false,
                PREVIOUS_WORKFLOW_OUTPUT= if DRAGEN_MODE then siblingCallWorkflow2.output_vcf else "null"
        }
    }
    File proband_gvcf = probandCallWorkflow.output_vcf
    File? sibling1_gvcf = siblingCallWorkflow1.output_vcf
    File? sibling2_gvcf = siblingCallWorkflow2.output_vcf
    File? sibling3_gvcf = siblingCallWorkflow3.output_vcf
    Array[File] gvcf_files_siblings = select_all([proband_gvcf, sibling1_gvcf, sibling2_gvcf, sibling3_gvcf])
    
    #######################################################
    ## Run 2nd trio joint genotyping on new proband GVCF ##
    #######################################################
    if (!DRAGEN_MODE) {
        call runGATKCombineGenotypeGVCFs as gatkJointGenotyper2nd {
            input:
                in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
                in_gvcf_file_maternal=MATERNAL_GVCF,
                in_gvcf_file_paternal=PATERNAL_GVCF,
                in_gvcf_files_siblings=gvcf_files_siblings,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzip2ndGATKGVCF {
            input:
                in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
                in_merged_vcf_file=gatkJointGenotyper2nd.joint_genotyped_vcf,
                in_vg_container=VG_CONTAINER
        }
    }
    if (DRAGEN_MODE) {
        call runDragenJointGenotyper as dragenJointGenotyper2nd {
            input:
                in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
                in_gvcf_file_maternal=MATERNAL_GVCF,
                in_gvcf_file_paternal=PATERNAL_GVCF,
                in_gvcf_files_siblings=gvcf_files_siblings,
                in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                in_udp_data_dir=UDPBINFO_PATH,
                in_helix_username=HELIX_USERNAME
        }
    }
    File output_2nd_joint_genotyped_vcf = select_first([bgzip2ndGATKGVCF.output_merged_vcf, dragenJointGenotyper2nd.dragen_joint_genotyped_vcf])
    File output_2nd_joint_genotyped_vcf_index = select_first([bgzip2ndGATKGVCF.output_merged_vcf_index, dragenJointGenotyper2nd.dragen_joint_genotyped_vcf_index])
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION) {
        call vgMultiMapCallWorkflow.normalizeVCF as normalizeCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
                in_bgzip_vcf_file=output_2nd_joint_genotyped_vcf,
        }
        call vgMultiMapCallWorkflow.snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
                in_normalized_vcf_file=normalizeCohortVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
        }
    }
    if (!SNPEFF_ANNOTATION) {
        File final_vcf_output = output_2nd_joint_genotyped_vcf
    }
    
    output {
        File output_cohort_vcf = select_first([snpEffAnnotateCohortVCF.output_snpeff_annotated_vcf, final_vcf_output])
        Array[File] output_gvcf_files_siblings = gvcf_files_siblings
    }
}

########################
### TASK DEFINITIONS ###
########################

task runGATKCombineGenotypeGVCFs {
    input {
        String in_sample_name
        File in_gvcf_file_maternal
        File in_gvcf_file_paternal
        Array[File]+ in_gvcf_files_siblings
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
    }
    
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        gatk IndexFeatureFile \
            -F ${in_gvcf_file_maternal} \
        && gatk IndexFeatureFile \
            -F ${in_gvcf_file_paternal}
        for sibling_gvcf_file in ${sep=" " in_gvcf_files_siblings} ; do
            gatk IndexFeatureFile -F "$sibling_gvcf_file"
        done
        gatk CombineGVCFs \
          --reference ${in_reference_file} \
          -V ${in_gvcf_file_maternal} -V ${in_gvcf_file_paternal} -V ${sep=" -V " in_gvcf_files_siblings} \
          --output ${in_sample_name}_cohort.combined.gvcf \
        && gatk GenotypeGVCFs \
          --reference ${in_reference_file} \
          --variant ${in_sample_name}_cohort.combined.gvcf \
          --output ${in_sample_name}_cohort.jointgenotyped.vcf
    }
    output {
        File joint_genotyped_vcf = "${in_sample_name}_cohort.jointgenotyped.vcf"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runDragenJointGenotyper {
    input {
        String in_sample_name
        File in_gvcf_file_maternal
        File in_gvcf_file_paternal
        Array[File]+ in_gvcf_files_siblings
        String in_dragen_ref_index_name
        String in_udp_data_dir
        String in_helix_username
    }

    String maternal_gvcf_file_name = basename(in_gvcf_file_maternal)
    String paternal_gvcf_file_name = basename(in_gvcf_file_paternal)

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        
        ## Copy input GVCFs into directory that Dragen can access
        UDP_DATA_DIR_PATH="~{in_udp_data_dir}/usr/~{in_helix_username}"
        DRAGEN_SIBLING_VCF_INPUT=""
        for sibling_gvcf_file in ~{sep=" " in_gvcf_files_siblings} ; do
            SIBLING_BASENAME="$(basename $sibling_gvcf_file)"
            DRAGEN_SIBLING_VCF_INPUT+="--variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/${SIBLING_BASENAME} "
        done
        mkdir -p /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        cp ~{in_gvcf_file_maternal} ~{in_gvcf_file_paternal} ~{sep=" " in_gvcf_files_siblings} /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        JOINT_GENOTYPE_DRAGEN_WORK_DIR="/staging/~{in_helix_username}/output_cohort_joint_call_~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR} --output-file-prefix cohort_joint_genotyped_~{in_sample_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{maternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{paternal_gvcf_file_name} ${DRAGEN_SIBLING_VCF_INPUT}" && \
        mkdir /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && chmod ug+rw -R /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/. /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/" && \
        mv /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper ~{in_sample_name}_dragen_joint_genotyper && \
        rm -fr /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_joint_genotyped_vcf = "~{in_sample_name}_dragen_joint_genotyper/cohort_joint_genotyped_~{in_sample_name}.vcf.gz"
        File dragen_joint_genotyped_vcf_index = "~{in_sample_name}_dragen_joint_genotyper/cohort_joint_genotyped_~{in_sample_name}.vcf.gz.tbi"
    }
    runtime {
        memory: 50 + " GB"
    }
}

