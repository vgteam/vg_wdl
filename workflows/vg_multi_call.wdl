version 1.0

### vg_multi_call.wdl ###
## Author: Charles Markello
## Description: Core VG variant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMapCall {
    input {
        Boolean SURJECT_MODE = true                 # Set to 'true' to run pipeline using alignmed BAM files surjected from GAM. Set to 'false' to output graph aligned GAM files.
        Boolean DRAGEN_MODE = false                 # Set to 'true' to use the Dragen modules variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
        Boolean GVCF_MODE = false                   # Set to 'true' to process and output gVCFs instead of VCFs.
        Boolean SNPEFF_ANNOTATION = true            # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        Boolean SV_CALLER_MODE = false              # Set to 'true' to run structural variant calling from graph aligned GAMs (SURJECT_MODE must be 'false' for this feature to be used)
        Boolean GOOGLE_CLEANUP_MODE = false         # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        File? INPUT_BAM_FILE                        # Input sample surjected .bam file
        File? INPUT_BAM_FILE_INDEX                  # Input sample .bai index of surjected .bam file.
        File? INPUT_GAM_FILE                        # Input sample .gam file
        File? INPUT_GAM_FILE_INDEX                  # Input sample .gai index of .gam file
        String PREVIOUS_WORKFLOW_OUTPUT = "null"    # Dummy input for iterative workflow dependency functionality.
        String SAMPLE_NAME                          # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        String PCR_INDEL_MODEL = "CONSERVATIVE"     # PCR indel model used in GATK Haplotypecaller (NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE)
        Int CHUNK_BASES = 50000000                  # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                          # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                        # (OPTIONAL) Text file where each line is a path name in the XG index
        File? XG_FILE                               # Path to .xg index file
        File REF_FILE                               # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                         # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                          # Path to .dict file of the REF_FILE fasta reference
        File SNPEFF_DATABASE                        # Path to snpeff database .zip file for snpEff annotation functionality.
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 400
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
        String DRAGEN_REF_INDEX_NAME                # Dragen module based reference index directory (e.g. "hs37d5_v7")
        String UDPBINFO_PATH                        # Udp data directory to use for Dragen module (e.g. "Udpbinfo", nih biowulf system only)
        String HELIX_USERNAME                       # The nih helix username which holds a user directory in UDPBINFO_PATH
    }
    
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    if (!defined(PATH_LIST_FILE) && defined(XG_FILE)) {
        call extractPathNames {
            input:
                in_xg_file=XG_FILE,
                in_vg_container=VG_CONTAINER
        }
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractPathNames.output_path_list])
    
    File input_alignment_file = select_first([INPUT_BAM_FILE, INPUT_GAM_FILE])
    File input_alignment_file_index = select_first([INPUT_BAM_FILE_INDEX, INPUT_GAM_FILE_INDEX])
    ##############################################
    # Run the linear alignment calling procedure #
    ##############################################
    if (SURJECT_MODE) {
        if ((!GVCF_MODE) && (!DRAGEN_MODE)) {
            # Split merged alignment by contigs list
            call splitBAMbyPath {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_merged_bam_file=input_alignment_file,
                    in_merged_bam_file_index=input_alignment_file_index,
                    in_path_list_file=pipeline_path_list_file
            }
            # Run distributed GATK linear variant calling
            scatter (gatk_caller_input_files in splitBAMbyPath.bams_and_indexes_by_contig) {
                # Run regular VCF genotypers
                call runGATKHaplotypeCaller {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=gatk_caller_input_files.left,
                        in_bam_file_index=gatk_caller_input_files.right,
                        in_pcr_indel_model=PCR_INDEL_MODEL,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
                # Cleanup intermediate variant calling files after use
                if (GOOGLE_CLEANUP_MODE) {
                    call cleanUpGoogleFilestore as cleanUpLinearCallerInputsGoogle {
                        input:
                            previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                            current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                    }
                }
                if (!GOOGLE_CLEANUP_MODE) {
                    call cleanUpUnixFilesystem as cleanUpLinearCallerInputsUnix {
                        input:
                            previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                            current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                    }
                }
            }

            # Merge linear-based called VCFs
            Array[File?] final_contig_vcf_output = runGATKHaplotypeCaller.genotyped_vcf
            Array[File] final_contig_vcf_output_list = select_all(final_contig_vcf_output)
            call concatClippedVCFChunks as concatLinearVCFChunks {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_clipped_vcf_chunk_files=final_contig_vcf_output_list
            }
            call bgzipMergedVCF as bgzipGATKCalledVCF {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_merged_vcf_file=concatLinearVCFChunks.output_merged_vcf,
                    in_vg_container=VG_CONTAINER
            }
        }
        # Run VCF variant calling using the Dragen module
        if ((!GVCF_MODE) && (DRAGEN_MODE)) {
            call runDragenCaller {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=input_alignment_file,
                    in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                    in_udp_data_dir=UDPBINFO_PATH,
                    in_helix_username=HELIX_USERNAME
            }
        }
        # Run GVCF variant calling procedure
        if (GVCF_MODE) {
            # Run GVCF genotypers
            if (!DRAGEN_MODE) {
                call runGATKHaplotypeCallerGVCF {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=input_alignment_file,
                        in_bam_file_index=input_alignment_file_index,
                        in_pcr_indel_model=PCR_INDEL_MODEL,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
            }
            # Run GVCF varaint calling using the Dragen module
            if (DRAGEN_MODE) {
                call runDragenCallerGVCF {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=input_alignment_file,
                        in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                        in_udp_data_dir=UDPBINFO_PATH,
                        in_helix_username=HELIX_USERNAME
                }
            }

            File final_gvcf_output = select_first([runGATKHaplotypeCallerGVCF.rawLikelihoods_gvcf, runDragenCallerGVCF.dragen_genotyped_gvcf])
        }
    }

    ################################
    # Run the VG calling procedure #
    ################################
    if ((!SURJECT_MODE) && defined(XG_FILE)) {
        # Set chunk context amount depending on structural variant or default variant calling mode
        Int default_or_sv_chunk_context = if SV_CALLER_MODE then 200 else 50
        # Chunk GAM alignments by contigs list
        call chunkAlignmentsByPathNames {
            input:
            in_sample_name=SAMPLE_NAME,
            in_merged_sorted_gam=input_alignment_file,
            in_merged_sorted_gam_gai=input_alignment_file_index,
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_chunk_context=default_or_sv_chunk_context,
            in_chunk_bases=CHUNK_BASES,
            in_overlap=OVERLAP,
            in_vg_container=VG_CONTAINER,
            in_chunk_gam_cores=CHUNK_GAM_CORES,
            in_chunk_gam_disk=CHUNK_GAM_DISK,
            in_chunk_gam_mem=CHUNK_GAM_MEM
        }
        Array[Pair[File,File]] vg_caller_input_files_list = zip(chunkAlignmentsByPathNames.output_vg_chunks, chunkAlignmentsByPathNames.output_gam_chunks)
        # Run distributed graph-based variant calling
        scatter (vg_caller_input_files in vg_caller_input_files_list) {
            call runVGCaller {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_chunk_bed_file=chunkAlignmentsByPathNames.output_bed_chunk_file,
                    in_vg_file=vg_caller_input_files.left,
                    in_gam_file=vg_caller_input_files.right,
                    in_chunk_bases=CHUNK_BASES,
                    in_overlap=OVERLAP,
                    in_vg_container=VG_CONTAINER,
                    in_vgcall_cores=VGCALL_CORES,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM,
                    in_sv_mode=SV_CALLER_MODE
            }
            call runVCFClipper {
                input:
                    in_chunk_vcf=runVGCaller.output_vcf,
                    in_chunk_vcf_index=runVGCaller.output_vcf_index,
                    in_chunk_clip_string=runVGCaller.clip_string
            }
            # Cleanup vg call input files after use
            if (GOOGLE_CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpVGCallInputsGoogle {
                    input:
                        previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                        current_task_output = runVCFClipper.output_clipped_vcf
                }
            }
            if (!GOOGLE_CLEANUP_MODE) {
                call cleanUpUnixFilesystem as cleanUpVGCallInputsUnix {
                    input:
                        previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                        current_task_output = runVCFClipper.output_clipped_vcf
                }
            }
        }
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksVGCall {
            input:
                in_sample_name=SAMPLE_NAME,
                in_clipped_vcf_chunk_files=runVCFClipper.output_clipped_vcf
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_vcf_file=concatVCFChunksVGCall.output_merged_vcf,
                in_vg_container=VG_CONTAINER
        }
    }

    # Extract either the linear-based or graph-based VCF
    File variantcaller_vcf_output = select_first([runDragenCaller.dragen_genotyped_vcf, bgzipVGCalledVCF.output_merged_vcf, bgzipGATKCalledVCF.output_merged_vcf, final_gvcf_output])
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION) {
        call normalizeVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bgzip_vcf_file=variantcaller_vcf_output,
        }
        call snpEffAnnotateVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_normalized_vcf_file=normalizeVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
        }
    }
    if (!SNPEFF_ANNOTATION) {
        File final_vcf_output = variantcaller_vcf_output
    }

    output {
        File output_vcf = select_first([snpEffAnnotateVCF.output_snpeff_annotated_vcf, final_vcf_output])
    }
}

########################
### TASK DEFINITIONS ###
########################

# Tasks for intermediate file cleanup
task cleanUpUnixFilesystem {
    input {
        Array[String] previous_task_outputs
        String current_task_output
    }
    command <<<
        set -eux -o pipefail
        cat ~{write_lines(previous_task_outputs)} | sed 's/.*\(\/cromwell-executions\)/\1/g' | xargs -I {} ls -li {} | cut -f 1 -d ' ' | xargs -I {} find ../../../ -xdev -inum {} | xargs -I {} rm -v {}
    >>>
    runtime {
        time: 10
        docker: "null"
        continueOnReturnCode: true
    }
}

task cleanUpGoogleFilestore {
    input {
        Array[String] previous_task_outputs
        String current_task_output
    }
    command {
        set -eux -o pipefail
        gsutil rm -I < ${write_lines(previous_task_outputs)}
    }
    runtime {
        docker: "google/cloud-sdk"
        continueOnReturnCode: true
    }
}

task extractPathNames {
    input {
        File? in_xg_file
        String in_vg_container
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt
    }
    output {
        File output_path_list = "path_list.txt"
    }
    runtime {
        memory: "50 GB"
        disks: "local-disk 50 SSD"
        docker: in_vg_container
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai
        
        while IFS=$'\t' read -ra path_list_line; do
            path_name="${path_list_line[0]}"
            samtools view \
              -@ 32 \
              -h -O BAM \
              input_bam_file.bam ${path_name} > ~{in_sample_name}.${path_name}.bam \
            && samtools index \
              ~{in_sample_name}.${path_name}.bam
        done < ~{in_path_list_file}
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
        Array[Pair[File, File]] bams_and_indexes_by_contig = zip(bam_contig_files, bam_contig_files_index)
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task runGATKHaplotypeCaller {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        String in_pcr_indel_model
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
        
        ln -s ${in_bam_file} input_bam_file.bam
        ln -s ${in_bam_file_index} input_bam_file.bam.bai
        
        gatk HaplotypeCaller \
          --native-pair-hmm-threads 32 \
          --pcr-indel-model ${in_pcr_indel_model} \
          --reference ${in_reference_file} \
          --input input_bam_file.bam \
          --output ${in_sample_name}.vcf \
        && bgzip ${in_sample_name}.vcf
    }
    output {
        File genotyped_vcf = "${in_sample_name}.vcf.gz"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runGATKHaplotypeCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        String in_pcr_indel_model
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
        
        ln -s ${in_bam_file} input_bam_file.bam
        ln -s ${in_bam_file_index} input_bam_file.bam.bai

        gatk HaplotypeCaller \
          --native-pair-hmm-threads 32 \
          -ERC GVCF \
          --pcr-indel-model ${in_pcr_indel_model} \
          --reference ${in_reference_file} \
          --input input_bam_file.bam \
          --output ${in_sample_name}.rawLikelihoods.gvcf \
        && bgzip ${in_sample_name}.rawLikelihoods.gvcf
    }
    output {
        File rawLikelihoods_gvcf = "${in_sample_name}.rawLikelihoods.gvcf.gz"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runDragenCaller {
    input {
        String in_sample_name
        File in_bam_file
        String in_dragen_ref_index_name
        String in_udp_data_dir
        String in_helix_username
    }

    String bam_file_name = basename(in_bam_file)

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

        UDP_DATA_DIR_PATH="~{in_udp_data_dir}/usr/~{in_helix_username}"
        mkdir -p /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/ && \
        cp ~{in_bam_file} /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/ && \
        DRAGEN_WORK_DIR_PATH="/staging/~{in_helix_username}/~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${DRAGEN_WORK_DIR_PATH}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} -b /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/~{bam_file_name} --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-sample-name ~{in_sample_name} --intermediate-results-dir ${TMP_DIR} --output-directory ${DRAGEN_WORK_DIR_PATH} --output-file-prefix ~{in_sample_name}_dragen_genotyped" && \
        mkdir /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper && chmod ug+rw -R /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${DRAGEN_WORK_DIR_PATH}/. /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${DRAGEN_WORK_DIR_PATH}/" && \
        mv /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper ~{in_sample_name}_dragen_genotyper && \
        rm -fr /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_genotyped_vcf = "~{in_sample_name}_dragen_genotyper/~{in_sample_name}_dragen_genotyped.vcf.gz"
    }
    runtime {
        memory: 50 + " GB"
    }
}

task runDragenCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        String in_dragen_ref_index_name
        String in_udp_data_dir
        String in_helix_username
    }

    String bam_file_name = basename(in_bam_file)

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

        UDP_DATA_DIR_PATH="~{in_udp_data_dir}/usr/~{in_helix_username}"
        mkdir -p /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/ && \
        cp ~{in_bam_file} /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/ && \
        DRAGEN_WORK_DIR_PATH="/staging/~{in_helix_username}/~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${DRAGEN_WORK_DIR_PATH}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} -b /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/~{bam_file_name} --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-emit-ref-confidence GVCF --vc-sample-name ~{in_sample_name} --intermediate-results-dir ${TMP_DIR} --output-directory ${DRAGEN_WORK_DIR_PATH} --output-file-prefix ~{in_sample_name}_dragen_genotyped" && \
        mkdir /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper && chmod ug+rw -R /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${DRAGEN_WORK_DIR_PATH}/. /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper" && \
        ssh ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${DRAGEN_WORK_DIR_PATH}/" && \
        mv /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_genotyper ~{in_sample_name}_dragen_genotyper && \
        rm -fr /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_genotyped_gvcf = "~{in_sample_name}_dragen_genotyper/~{in_sample_name}_dragen_genotyped.gvcf.gz"
    }
    runtime {
        memory: 50 + " GB"
    }
}

task chunkAlignmentsByPathNames {
    input {
        String in_sample_name
        File in_merged_sorted_gam
        File in_merged_sorted_gam_gai
        File? in_xg_file
        Int in_chunk_context
        File in_path_list_file
        Int in_chunk_bases
        Int in_overlap
        String in_vg_container
        Int in_chunk_gam_cores
        Int in_chunk_gam_disk
        Int in_chunk_gam_mem
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

        vg chunk \
            -x ${in_xg_file} \
            -a ${in_merged_sorted_gam} \
            -c ${in_chunk_context} \
            -P ${in_path_list_file} \
            -g \
            -s ${in_chunk_bases} -o ${in_overlap} \
            -b call_chunk \
            -t ${in_chunk_gam_cores} \
            -E output_bed_chunks.bed -f
    }
    output {
        Array[File] output_vg_chunks = glob("call_chunk*.vg")
        Array[File] output_gam_chunks = glob("call_chunk*.gam")
        File output_bed_chunk_file = "output_bed_chunks.bed"
    }
    runtime {
        time: 900
        memory: in_chunk_gam_mem + " GB"
        cpu: in_chunk_gam_cores
        disks: "local-disk " + in_chunk_gam_disk + " SSD"
        docker: in_vg_container
    }
}

task runVGCaller {
    input {
        String in_sample_name
        File in_chunk_bed_file
        File in_vg_file
        File in_gam_file
        Int in_chunk_bases
        Int in_overlap
        String in_vg_container
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
        Boolean in_sv_mode
    }

    String chunk_tag = basename(in_vg_file, ".vg")
    command <<<
        set -eux -o pipefail

        VG_INDEX_XG_COMMAND=""
        VG_FILTER_COMMAND=""
        VG_AUGMENT_SV_OPTIONS=""
        VG_CALL_SV_OPTIONS=""
        if [ ~{in_sv_mode} == false ]; then
            VG_INDEX_XG_COMMAND="vg index ~{in_vg_file} -x ~{chunk_tag}.xg -t ~{in_vgcall_cores} && \\"
            VG_FILTER_COMMAND="vg filter ~{in_gam_file} -t 1 -r 0.9 -fu -s 1000 -m 1 -q 15 -D 999 -x ~{chunk_tag}.xg > ~{chunk_tag}.filtered.gam && \\"
            VG_AUGMENT_SV_OPTIONS="~{chunk_tag}.filtered_gam"
        else
            VG_AUGMENT_SV_OPTIONS="~{in_gam_file} --recall"
            VG_CALL_SV_OPTIONS="-u -n 0 -e 1000 -G 3"
        fi

        PATH_NAME="$(echo ~{chunk_tag} | cut -f 4 -d '_')"
        OFFSET=""
        BED_CHUNK_LINE_RECORD=""
        FIRST_Q=""
        LAST_Q=""
        CHUNK_INDEX_COUNTER="0"
        CHUNK_TAG_INDEX="0"
        CHR_LENGTH=""
        while IFS=$'\t' read -ra bed_chunk_line; do
            bed_line_chunk_tag="$(echo ${bed_chunk_line[3]} | cut -f 1 -d '.')"
            if [ "~{chunk_tag}" = "${bed_line_chunk_tag}" ]; then
                OFFSET="${bed_chunk_line[1]}"
                BED_CHUNK_LINE_RECORD="${bed_chunk_line[@]//$'\t'/ }"
                CHUNK_TAG_INDEX="${CHUNK_INDEX_COUNTER}"
            fi
            bed_line_path_name="$(echo ${bed_line_chunk_tag} | cut -f 4 -d '_')"
            if [ "${PATH_NAME}" = "${bed_line_path_name}" ]; then
                if [ "${FIRST_Q}" = "" ]; then
                    FIRST_Q="${bed_chunk_line[@]//$'\t'/ }"
                fi
                LAST_Q="${bed_chunk_line[@]//$'\t'/ }"
                CHUNK_INDEX_COUNTER="$(( ${CHUNK_INDEX_COUNTER} + 1 ))"
                CHR_LENGTH="${bed_chunk_line[2]}"
            fi
        done < ~{in_chunk_bed_file}

        ${VG_INDEX_XG_COMMAND}
        ${VG_FILTER_COMMAND}
        vg augment \
            ~{in_vg_file} \
            ${VG_AUGMENT_SV_OPTIONS} \
            -t ~{in_vgcall_cores} \
            -a pileup \
            -Z ~{chunk_tag}.trans \
            -S ~{chunk_tag}.support > ~{chunk_tag}.aug.vg && \
        vg call \
            ~{chunk_tag}.aug.vg \
            -t ~{in_vgcall_cores} \
            -S ~{in_sample_name} \
            -z ~{chunk_tag}.trans \
            -s ~{chunk_tag}.support \
            -b ~{in_vg_file} \
            -r ${PATH_NAME} \
            -c ${PATH_NAME} \
            -l ${CHR_LENGTH} \
            ${VG_CALL_SV_OPTIONS} \
            -o ${OFFSET} > ~{chunk_tag}.vcf && \
        head -10000 ~{chunk_tag}.vcf | grep "^#" >> ~{chunk_tag}.sorted.vcf && \
        if [ "$(cat ~{chunk_tag}.vcf | grep -v '^#')" ]; then
            cat ~{chunk_tag}.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ~{chunk_tag}.sorted.vcf
        fi
        bgzip ~{chunk_tag}.sorted.vcf && \
        tabix -f -p vcf ~{chunk_tag}.sorted.vcf.gz

        # Compile clipping string
        clipped_chunk_offset="$(( (${CHUNK_TAG_INDEX} * ~{in_chunk_bases}) - (${CHUNK_TAG_INDEX} * ~{in_overlap}) ))"
        if [ "${BED_CHUNK_LINE_RECORD}" = "${FIRST_Q}" ]; then
          START="$(( ${clipped_chunk_offset} + 1 ))"
        else
          START="$(( ${clipped_chunk_offset} + 1 + ~{in_overlap}/2 ))"
        fi
        if [ "${BED_CHUNK_LINE_RECORD}" = "${LAST_Q}" ]; then
          END="$(( ${clipped_chunk_offset} + ~{in_chunk_bases}))"
        else
          END="$(( ${clipped_chunk_offset} + ~{in_chunk_bases} - (~{in_overlap}/2) ))"
        fi
        echo "${PATH_NAME}:${START}-${END}" > ~{chunk_tag}_clip_string.txt
    >>>
    output {
        File output_vcf = "~{chunk_tag}.sorted.vcf.gz"
        File output_vcf_index = "~{chunk_tag}.sorted.vcf.gz.tbi"
        String clip_string = read_lines("~{chunk_tag}_clip_string.txt")[0]
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: in_vg_container
    }
}

task runVCFClipper {
    input {
        File in_chunk_vcf
        File in_chunk_vcf_index
        String in_chunk_clip_string
    }

    String chunk_tag = basename(in_chunk_vcf, ".sorted.vcf.gz")
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

        bcftools view -O b ${in_chunk_vcf} -t ${in_chunk_clip_string} > ${chunk_tag}.clipped.vcf.gz
    }
    output {
        File output_clipped_vcf = "${chunk_tag}.clipped.vcf.gz"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
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

        for vcf_file in ${sep=" " in_clipped_vcf_chunk_files} ; do
            bcftools index "$vcf_file"
        done
        bcftools concat -a ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort - > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
    }

    # TODO:
    #   If GVCF in in_merged_vcf_file then output_vcf_extension="gvcf" else output_vcf_extension="vcf"
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

        bgzip -c ${in_merged_vcf_file} > ${in_sample_name}_merged.vcf.gz && \
        tabix -f -p vcf ${in_sample_name}_merged.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}_merged.vcf.gz.tbi"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
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

        bcftools norm -m-both --threads 16 -o ${in_sample_name}.unrolled.vcf ${in_bgzip_vcf_file}
    }
    output {
        File output_normalized_vcf = "${in_sample_name}.unrolled.vcf"
    }
    runtime {
        cpu: 16
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File in_snpeff_database
    }

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

        unzip ~{in_snpeff_database}
        snpEff -Xmx40g -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data GRCh37.75 ~{in_normalized_vcf_file} > ~{in_sample_name}.snpeff.unrolled.vcf
    >>>
    output {
        File output_snpeff_annotated_vcf = "~{in_sample_name}.snpeff.unrolled.vcf"
    }
    runtime {
        cpu: 16
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/snpeff:4.3.1t--2"
    }
}


