version 1.0

### vg_trio_multi_call.wdl ###
# Author: Charles Markello
# Description: Variant calling workflow for mother-father-proband trios.
#              Designed as the 2nd step in a pedigree-backed graph alignment pipeline.

import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "./vg_multi_call.wdl" as vgMultiCallWorkflow

###########################
### WORKFLOW DEFINITION ###
###########################
workflow vgTrioPipeline {
    input {
        File MATERNAL_INPUT_BAM_FILE                        # Input maternal surjected .bam file
        File MATERNAL_INPUT_BAM_FILE_INDEX                  # Input maternal .bai index of surjected .bam file
        File PATERNAL_INPUT_BAM_FILE                        # Input paternal surjected .bam file
        File PATERNAL_INPUT_BAM_FILE_INDEX                  # Input paternal .bai index of surjected .bam file
        File PROBAND_INPUT_BAM_FILE                         # Input proband surjected .bam file
        File PROBAND_INPUT_BAM_FILE_INDEX                   # Input proband .bai index of surjected .bam file
        String SAMPLE_NAME_MATERNAL                         # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                         # Sample name for the father
        String SAMPLE_NAME_PROBAND                          # Sample name for the proband
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int CHUNK_BASES = 50000000                          # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                                  # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        File? SNPEFF_DATABASE                                # Path to snpeff database .zip file
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 400
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
        String DRAGEN_REF_INDEX_NAME                        # Dragen module based reference index directory (e.g. "hs37d5_v7")
        String UDPBINFO_PATH                                # Udp data directory to use for Dragen module (e.g. "Udpbinfo", nih biowulf system only)
        String HELIX_USERNAME                               # The nih helix username which holds a user directory in UDPBINFO_PATH
        Boolean DRAGEN_MODE = false                         # Set to 'true' to use the Dragen modules variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
    }
    
    ###########################################
    ## Run variant calling workflows on Trio ##
    ###########################################
    call vgMultiCallWorkflow.vgMultiMapCall as maternalCallWorkflow {
        input:
            INPUT_BAM_FILE=MATERNAL_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=MATERNAL_INPUT_BAM_FILE_INDEX,
            SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
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
            DRAGEN_MODE=false,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false
    }
    call vgMultiCallWorkflow.vgMultiMapCall as paternalCallWorkflow {
        input:
            INPUT_BAM_FILE=PATERNAL_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=PATERNAL_INPUT_BAM_FILE_INDEX,
            SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
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
            DRAGEN_MODE=false,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT="null"
    }
    call vgMultiCallWorkflow.vgMultiMapCall as probandCallWorkflow {
        input:
            INPUT_BAM_FILE=PROBAND_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=PROBAND_INPUT_BAM_FILE_INDEX,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
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
            DRAGEN_MODE=false,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT="null"
    }
    
    ###############################
    ## Run trio joint genotyping ##
    ###############################
    if (!DRAGEN_MODE) {
        call runGATKCombineGenotypeGVCFs as gatkJointGenotyper1st {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandCallWorkflow.output_vcf,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipGATKGVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_merged_vcf_file=gatkJointGenotyper1st.joint_genotyped_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    if (DRAGEN_MODE) {
        call runDragenJointGenotyper as dragenJointGenotyper1st {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandCallWorkflow.output_vcf,
                in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                in_udp_data_dir=UDPBINFO_PATH,
                in_helix_username=HELIX_USERNAME
        }
    }
    
    File output_joint_genotyped_vcf = select_first([bgzipGATKGVCF.output_merged_vcf, dragenJointGenotyper1st.dragen_joint_genotyped_vcf])
    File output_joint_genotyped_vcf_index = select_first([bgzipGATKGVCF.output_merged_vcf_index, dragenJointGenotyper1st.dragen_joint_genotyped_vcf_index])
    
    output {
        File output_cohort_vcf = output_joint_genotyped_vcf
        File output_cohort_vcf_index = output_joint_genotyped_vcf_index
        File output_maternal_vcf = maternalCallWorkflow.output_vcf
        File output_paternal_vcf = paternalCallWorkflow.output_vcf
        File output_proband_vcf = probandCallWorkflow.output_vcf
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
        File in_gvcf_file_proband
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
            -F ${in_gvcf_file_paternal} \
        && gatk IndexFeatureFile \
            -F ${in_gvcf_file_proband} \
        && gatk CombineGVCFs \
          --reference ${in_reference_file} \
          -V ${in_gvcf_file_maternal} -V ${in_gvcf_file_paternal} -V ${in_gvcf_file_proband} \
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
        File in_gvcf_file_proband
        String in_dragen_ref_index_name
        String in_udp_data_dir
        String in_helix_username
    }

    String maternal_gvcf_file_name = basename(in_gvcf_file_maternal)
    String paternal_gvcf_file_name = basename(in_gvcf_file_paternal)
    String proband_gvcf_file_name = basename(in_gvcf_file_proband)

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
        mkdir -p /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        cp ~{in_gvcf_file_maternal} ~{in_gvcf_file_paternal} ~{in_gvcf_file_proband} /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        JOINT_GENOTYPE_DRAGEN_WORK_DIR="/staging/~{in_helix_username}/output_cohort_joint_call_~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR} --output-file-prefix cohort_joint_genotyped_~{in_sample_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{maternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{paternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{proband_gvcf_file_name}" && \
        mkdir /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && chmod ug+rw -R /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/. /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/" && \
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

