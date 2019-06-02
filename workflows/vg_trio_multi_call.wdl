version 1.0
#version draft-2
# Apparently cromwell 36.1 doesn't parse the version line in wdl files if the wdl version is a draft version.
# Cromwell 36.1 wdl parser defaults to draft-2 if a `version` line is not given.

# TODO:
# - Use subworkflow functionality in new vg_trio_sub.wdl workflow
# - possibly define genotype phasing task
#     - Start with just the regular mpmap gcsa and xg construction

import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "./vg_multi_map.wdl" as vgMultiMapWorkflow
import "./vg_multi_call.wdl" as vgMultiCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow

workflow vgTrioPipeline {
    input {
        File PROBAND_INPUT_READ_FILE_1
        File PROBAND_INPUT_READ_FILE_2
        File MATERNAL_INPUT_BAM_FILE
        File MATERNAL_INPUT_BAM_FILE_INDEX
        File PATERNAL_INPUT_BAM_FILE
        File PATERNAL_INPUT_BAM_FILE_INDEX
        File PROBAND_INPUT_BAM_FILE
        File PROBAND_INPUT_BAM_FILE_INDEX
        String SAMPLE_NAME_MATERNAL
        String SAMPLE_NAME_PATERNAL
        String SAMPLE_NAME_PROBAND
        String VG_CONTAINER
        Int READS_PER_CHUNK
        Int CHUNK_BASES
        Int OVERLAP
        File PATH_LIST_FILE
        File XG_FILE
        File GCSA_FILE
        File GCSA_LCP_FILE
        File GBWT_FILE
        File SNARLS_FILE
        File REF_FILE
        File REF_INDEX_FILE
        File REF_DICT_FILE
        File SNPEFF_DATABASE
        Int SPLIT_READ_CORES
        Int SPLIT_READ_DISK
        Int MAP_CORES
        Int MAP_DISK
        Int MAP_MEM
        Int MERGE_GAM_CORES
        Int MERGE_GAM_DISK
        Int MERGE_GAM_MEM
        Int MERGE_GAM_TIME
        Int CHUNK_GAM_CORES
        Int CHUNK_GAM_DISK
        Int CHUNK_GAM_MEM
        Int VGCALL_CORES
        Int VGCALL_DISK
        Int VGCALL_MEM
        String DRAGEN_REF_INDEX_NAME
        String UDPBINFO_PATH
        String HELIX_USERNAME
        Boolean RUN_VGMPMAP_ALGORITHM
        Boolean RUN_DRAGEN_CALLER
        Array[String]+ CONTIGS = ["20", "21"]
        File REF_FASTA_GZ
        File PED_FILE
        File GEN_MAP_FILES
        String GRAPH_NAME
        Boolean USE_HAPLOTYPES = true
        Boolean MAKE_SNARLS = true
        Boolean USE_DECOYS = true
        String DECOY_REGEX = ">GL"
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
            RUN_LINEAR_CALLER=true,
            RUN_DRAGEN_CALLER=RUN_DRAGEN_CALLER,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false
    }
    call vgMultiCallWorkflow.vgMultiMapCall as paternalCallWorkflow {
        input:
            INPUT_BAM_FILE=PATERNAL_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=PATERNAL_INPUT_BAM_FILE,
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
            RUN_LINEAR_CALLER=true,
            RUN_DRAGEN_CALLER=RUN_DRAGEN_CALLER,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT= if RUN_DRAGEN_CALLER then maternalCallWorkflow.output_vcf else "null"
    }
    call vgMultiCallWorkflow.vgMultiMapCall as probandCallWorkflow {
        input:
            INPUT_BAM_FILE=PROBAND_INPUT_BAM_FILE,
            INPUT_BAM_FILE_INDEX=PROBAND_INPUT_BAM_FILE,
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
            RUN_LINEAR_CALLER=true,
            RUN_DRAGEN_CALLER=RUN_DRAGEN_CALLER,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT= if RUN_DRAGEN_CALLER then paternalCallWorkflow.output_vcf else "null"
    }
    
    ###############################
    ## Run trio joint genotyping ##
    ###############################
    if (!RUN_DRAGEN_CALLER) {
        call runGATKCombineGenotypeGVCFs {
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
                in_merged_vcf_file=runGATKCombineGenotypeGVCFs.joint_genotyped_vcf,
                in_vg_container=VG_CONTAINER
        }
    }
    if (RUN_DRAGEN_CALLER) {
        call runDragenJointGenotyper {
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
    
    #####################################
    ## Run parental graph construction ##
    #####################################
    File output_joint_genotyped_vcf = select_first([bgzipGATKGVCF.output_merged_vcf, runDragenJointGenotyper.dragen_joint_genotyped_vcf])
    File output_joint_genotyped_vcf_index = select_first([bgzipGATKGVCF.output_merged_vcf_index, runDragenJointGenotyper.dragen_joint_genotyped_vcf_index])
    call runSplitJointGenotypedVCF as splitJointGenotypedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=output_joint_genotyped_vcf,
            joint_genotyped_vcf_index=output_joint_genotyped_vcf_index,
            contigs=CONTIGS,
            filter_parents=false
    }
    scatter (contig_pair in zip(CONTIGS, splitJointGenotypedVCF.contig_vcfs)) {
        call runWhatsHapPhasing {
            input:
                in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                joint_genotyped_vcf=contig_pair.right,
                in_maternal_bam=MATERNAL_INPUT_BAM_FILE,
                in_maternal_bam_index=MATERNAL_INPUT_BAM_FILE_INDEX,
                in_paternal_bam=PATERNAL_INPUT_BAM_FILE,
                in_paternal_bam_index=PATERNAL_INPUT_BAM_FILE_INDEX,
                in_proband_bam=PROBAND_INPUT_BAM_FILE,
                in_proband_bam_index=PROBAND_INPUT_BAM_FILE_INDEX,
                in_ped_file=PED_FILE,
                in_genetic_map=GEN_MAP_FILES,
                in_contig=contig_pair.left,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
    }
    call vgMultiMapCallWorkflow.concatClippedVCFChunks as concatCohortPhasedVCFs {
        input:
            in_sample_name="${SAMPLE_NAME_PROBAND}_cohort",
            in_clipped_vcf_chunk_files=runWhatsHapPhasing.phased_cohort_vcf
    }
    call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipConcatPhasedVCF {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_merged_vcf_file=concatCohortPhasedVCFs.output_merged_vcf,
            in_vg_container=VG_CONTAINER
    }
    call runSplitJointGenotypedVCF as splitPhasedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=bgzipConcatPhasedVCF.output_merged_vcf,
            joint_genotyped_vcf_index=bgzipConcatPhasedVCF.output_merged_vcf_index,
            contigs=CONTIGS,
            filter_parents=true
    }
    call vgConstructWorkflow.vg_construct_and_index as constructGraphIndexWorkflow {
        input:
            graph_name=GRAPH_NAME,
            ref_fasta_gz=REF_FASTA_GZ,
            contigs=CONTIGS,
            contigs_vcf_gz=splitPhasedVCF.contig_vcfs,
            use_haplotypes=USE_HAPLOTYPES,
            make_snarls=MAKE_SNARLS,
            use_decoys=USE_DECOYS,
            include_mt=USE_DECOYS,
            decoy_regex=DECOY_REGEX,
            vg_docker=VG_CONTAINER
    }
    
    ##################################################################################
    ## Run mapping and calling workflow of proband reads against the parental graph ##
    ##################################################################################
    call vgMultiMapCallWorkflow.vgMultiMapCall as probandMapCallWorkflow2ndIteration {
        input:
            INPUT_READ_FILE_1=PROBAND_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PROBAND_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=constructGraphIndexWorkflow.xg,
            GCSA_FILE=constructGraphIndexWorkflow.gcsa,
            GCSA_LCP_FILE=constructGraphIndexWorkflow.gcsa_lcp,
            GBWT_FILE=constructGraphIndexWorkflow.gbwt,
            SNARLS_FILE=constructGraphIndexWorkflow.snarls,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            RUN_VGMPMAP_ALGORITHM=RUN_VGMPMAP_ALGORITHM,
            RUN_DRAGEN_CALLER=RUN_DRAGEN_CALLER,
            RUN_LINEAR_CALLER=true,
            RUN_GVCF_OUTPUT=false,
            RUN_SNPEFF_ANNOTATION=true
    }
    
    output {
        File output_cohort_vcf = output_joint_genotyped_vcf
        File output_final_proband_vcf = probandMapCallWorkflow2ndIteration.output_vcf
        File? output_final_proband_bam = probandMapCallWorkflow2ndIteration.output_bam
        File? output_final_proband_bam_index = probandMapCallWorkflow2ndIteration.output_bam_index
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
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR} --output-file-prefix cohort_joint_genotyped_~{in_sample_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{maternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{paternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{proband_gvcf_file_name}" && \
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

# Split a vcf into a list of vcfs each representing a contig
task runSplitJointGenotypedVCF {
    input {
        String in_proband_sample_name
        String in_maternal_sample_name
        String in_paternal_sample_name
        File joint_genotyped_vcf
        File joint_genotyped_vcf_index
        Array[String]+ contigs
        Boolean filter_parents
    }

    command <<<
        set -exu -o pipefail
        
        if [ ~{filter_parents} == true ]; then
          SAMPLE_FILTER_STRING="-s ~{in_maternal_sample_name},~{in_paternal_sample_name}"
        else
          SAMPLE_FILTER_STRING=""
        fi

        while read -r contig; do
            bcftools view -O z -r "${contig}" ${SAMPLE_FILTER_STRING} ~{joint_genotyped_vcf} > "${contig}.vcf.gz"
            echo "${contig}.vcf.gz" >> contig_vcf_list.txt
        done < "~{write_lines(contigs)}"
    >>>
    output {
        Array[File]+ contig_vcfs = read_lines("contig_vcf_list.txt")
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task runWhatsHapPhasing {
    input {
        String in_cohort_sample_name
        File joint_genotyped_vcf
        File? in_maternal_bam
        File? in_maternal_bam_index
        File? in_paternal_bam
        File? in_paternal_bam_index
        File? in_proband_bam
        File? in_proband_bam_index
        File in_ped_file
        File in_genetic_map
        String in_contig
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
    }

    command <<<
        set -exu -o pipefail

        tar -xvf ~{in_genetic_map}
        if [[ ~{in_contig} == "Y" || ~{in_contig} == "MT" ]]; then
            GENMAP_OPTION_STRING=""
        elif [ ~{in_contig} == "X" ]; then
            GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chrX_nonPAR_combined_b37.txt --chromosome X"
        else
            GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chr~{in_contig}_combined_b37.txt --chromosome ~{in_contig}"
        fi
        whatshap phase \
            --reference ~{in_reference_file} \
            --indels \
            --ped ~{in_ped_file} \
            ${GENMAP_OPTION_STRING} \
            -o ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf \
            ~{joint_genotyped_vcf} ~{in_proband_bam} ~{in_maternal_bam} ~{in_paternal_bam} \
        && bgzip ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf
    >>>

    output {
        File phased_cohort_vcf = "~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf.gz"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/whatshap:0.18--py37h6bb024c_0"
    }
}


