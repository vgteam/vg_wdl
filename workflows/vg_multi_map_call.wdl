version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "./vg_multi_map.wdl" as vgMultiMapWorkflow

### vg_multi_map_call.wdl ###
workflow vgMultiMapCall {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        description: "Core VG mapping and variant calling workflow for single sample datasets."
    }

    input {
        String MAPPER = "GIRAFFE"                   # Set to 'MAP' to use the "VG MAP" algorithm, set to 'MPMAP' to use "VG MPMAP" algorithm, set to 'GIRAFFE' to use "VG GIRAFFE".
        Boolean SURJECT_MODE = true             # Set to 'true' to run pipeline using alignmed BAM files surjected from GAM. Set to 'false' to output graph aligned GAM files.
        Boolean DEEPVARIANT_MODE = false        # Set to 'true' to use the DeepVariant variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
        Boolean GVCF_MODE = false               # Set to 'true' to process and output gVCFs instead of VCFs.
        Boolean SNPEFF_ANNOTATION = true        # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        Boolean SV_CALLER_MODE = false          # Set to 'true' to run structural variant calling from graph aligned GAMs (SURJECT_MODE must be 'false' for this feature to be used)
        Boolean CLEANUP_FILES = true            # Set to 'false' to turn off intermediate file cleanup.
        Boolean GOOGLE_CLEANUP_MODE = true     # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        File INPUT_READ_FILE_1                  # Input sample 1st read pair fastq.gz
        File INPUT_READ_FILE_2                  # Input sample 2nd read pair fastq.gz
        String SAMPLE_NAME                      # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.64.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.64.0)
        String PCR_INDEL_MODEL = "CONSERVATIVE" # PCR indel model used in GATK Haplotypecaller (NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE)
        Int READS_PER_CHUNK = 20000000          # Number of reads contained in each mapping chunk (20000000 for wgs).
        Int CHUNK_BASES = 50000000              # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                      # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                    # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                            # Path to .xg index file
        File? GCSA_FILE                         # Path to .gcsa index file
        File? GCSA_LCP_FILE                     # Path to .gcsa.lcp index file
        File? GBWT_FILE                         # (OPTIONAL) Path to .gbwt index file
        File? GGBWT_FILE                        # (OPTIONAL) Path to .gg index file
        File? DIST_FILE                         # (OPTIONAL) Path to .dist index file
        File? MIN_FILE                          # (OPTIONAL) Path to .min index file
        File? SNARLS_FILE                       # (OPTIONAL) Path to .snarls index file
        File REF_FILE                           # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                     # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                      # Path to .dict file of the REF_FILE fasta reference
        String REFERENCE_PREFIX = ""            # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File? SNPEFF_DATABASE                   # Path to snpeff database .zip file for snpEff annotation functionality.
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 32
        Int MAP_DISK = 10
        Int MAP_MEM = 60
        Int MERGE_GAM_CORES = 56
        Int MERGE_GAM_DISK = 100
        Int MERGE_GAM_MEM = 40
        Int MERGE_GAM_TIME = 2400
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 100
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 80
    }

    # Call vgMultiMap workflow for read mapping and alignment
    call vgMultiMapWorkflow.vgMultiMap {
        input:
            MAPPER=MAPPER,
            SURJECT_MODE=SURJECT_MODE,
            INPUT_READ_FILE_1=INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CLEANUP_FILES=CLEANUP_FILES,
            GOOGLE_CLEANUP_MODE=GOOGLE_CLEANUP_MODE
    }

    ##############################################
    # Run the linear alignment calling procedure #
    ##############################################
    if (SURJECT_MODE) {
        # Run VCF variant calling procedure
        # Split merged alignment by contigs list
        call splitBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_bam_file=select_first([vgMultiMap.output_merged_bam]),
                in_merged_bam_file_index=select_first([vgMultiMap.output_merged_bam_index]),
                in_path_list_file=select_first([PATH_LIST_FILE, vgMultiMap.output_path_list]),
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        # Run distributed GATK linear variant calling
        scatter (gatk_caller_input_files in splitBAMbyPath.bams_and_indexes_by_contig) {
            if (!GVCF_MODE) {
                # Run regular VCF genotypers
                if (!DEEPVARIANT_MODE) {
                    call runGATKHaplotypeCaller {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_pcr_indel_model=PCR_INDEL_MODEL,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpGoogleFilestore as cleanUpGATKCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpUnixFilesystem as cleanUpGATKCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                            }
                        }
                    }
                }
                if (DEEPVARIANT_MODE) {
                    call runDeepVariantCaller {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpGoogleFilestore as cleanUpDeepVariantCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCaller.genotyped_vcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpUnixFilesystem as cleanUpDeepVariantCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCaller.genotyped_vcf
                            }
                        }
                    }
                }
            }
            # Run GVCF variant calling procedure
            if (GVCF_MODE) {
                # Run GVCF genotypers
                if (!DEEPVARIANT_MODE) {
                    call runGATKHaplotypeCallerGVCF {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_pcr_indel_model=PCR_INDEL_MODEL,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpGoogleFilestore as cleanUpGATKGVCFCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCallerGVCF.genotyped_gvcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpUnixFilesystem as cleanUpGATKGVCFCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCallerGVCF.genotyped_gvcf
                            }
                        }
                    }
                }
                if (DEEPVARIANT_MODE) {
                    call runDeepVariantCallerGVCF {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpGoogleFilestore as cleanUpDeepVariantGVCFCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCallerGVCF.genotyped_gvcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call vgMultiMapWorkflow.cleanUpUnixFilesystem as cleanUpDeepVariantGVCFCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCallerGVCF.genotyped_gvcf
                            }
                        }
                    }
                }
            }
            File final_contig_vcf = select_first([runGATKHaplotypeCaller.genotyped_vcf, runDeepVariantCaller.genotyped_vcf, runGATKHaplotypeCallerGVCF.genotyped_gvcf, runDeepVariantCallerGVCF.genotyped_gvcf])
        }
        # Merge linear-based called VCFs
        Array[File?] final_contig_vcf_output = final_contig_vcf
        Array[File] final_contig_vcf_output_list = select_all(final_contig_vcf_output)
        call concatClippedVCFChunks as concatLinearVCFChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_clipped_vcf_chunk_files=final_contig_vcf_output_list,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call bgzipMergedVCF as bgzipLinearCalledVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_vcf_file=concatLinearVCFChunks.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }

        File final_vcf_output = bgzipLinearCalledVCF.output_merged_vcf
    }

    ################################
    # Run the VG calling procedure #
    ################################
    if (!SURJECT_MODE) {
        # Set chunk context amount depending on structural variant or default variant calling mode
        Int default_or_sv_chunk_context = if SV_CALLER_MODE then 2500 else 50
        # Chunk GAM alignments by contigs list
        call chunkAlignmentsByPathNames {
            input:
            in_sample_name=SAMPLE_NAME,
            in_merged_sorted_gam=select_first([vgMultiMap.output_merged_gam]),
            in_merged_sorted_gam_gai=select_first([vgMultiMap.output_merged_gam_index]),
            in_xg_file=XG_FILE,
            in_path_list_file=select_first([PATH_LIST_FILE, vgMultiMap.output_path_list]),
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
                    in_chunk_clip_string=runVGCaller.clip_string,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM
            }
            # Cleanup vg call input files after use
            if (CLEANUP_FILES) {
                if (GOOGLE_CLEANUP_MODE) {
                    call vgMultiMapWorkflow.cleanUpGoogleFilestore as cleanUpVGCallInputsGoogle {
                        input:
                            previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                            current_task_output = runVCFClipper.output_clipped_vcf
                    }
                }
                if (!GOOGLE_CLEANUP_MODE) {
                    call vgMultiMapWorkflow.cleanUpUnixFilesystem as cleanUpVGCallInputsUnix {
                        input:
                            previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                            current_task_output = runVCFClipper.output_clipped_vcf
                    }
                }
            }
        }
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksVGCall {
            input:
                in_sample_name=SAMPLE_NAME,
                in_clipped_vcf_chunk_files=runVCFClipper.output_clipped_vcf,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_vcf_file=concatVCFChunksVGCall.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }

    # Extract either the linear-based or graph-based VCF
    File variantcaller_vcf_output = select_first([bgzipVGCalledVCF.output_merged_vcf, final_vcf_output])
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION && defined(SNPEFF_DATABASE)) {
        call normalizeVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bgzip_vcf_file=variantcaller_vcf_output,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call snpEffAnnotateVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_normalized_vcf_file=normalizeVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    output {
        File output_vcf = select_first([snpEffAnnotateVCF.output_snpeff_annotated_vcf, variantcaller_vcf_output])
        File? output_bam = vgMultiMap.output_merged_bam
        File? output_bam_index = vgMultiMap.output_merged_bam_index
        File? output_gam = vgMultiMap.output_merged_gam
        File? output_gam_index = vgMultiMap.output_merged_gam_index
    }
}

########################
### TASK DEFINITIONS ###
########################

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        while IFS=$'\t' read -ra path_list_line; do
            path_name="${path_list_line[0]}"
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${path_name} \
              -o ~{in_sample_name}.${path_name}.bam \
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
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
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
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        gatk HaplotypeCaller \
          --native-pair-hmm-threads ~{in_vgcall_cores} \
          --pcr-indel-model ~{in_pcr_indel_model} \
          -L ${CONTIG_ID} \
          --reference ~{in_reference_file} \
          --input input_bam_file.bam \
          --output ~{in_sample_name}.vcf \
        && bgzip ~{in_sample_name}.vcf
    >>>
    output {
        File genotyped_vcf = "~{in_sample_name}.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "broadinstitute/gatk@sha256:cc8981d0527e716775645b04a7f59e96a52ad59a7ae9788ddc47902384bf35aa"
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
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        gatk HaplotypeCaller \
          --native-pair-hmm-threads ~{in_vgcall_cores} \
          -ERC GVCF \
          -L ${CONTIG_ID} \
          --pcr-indel-model ~{in_pcr_indel_model} \
          --reference ~{in_reference_file} \
          --input input_bam_file.bam \
          --output ~{in_sample_name}.rawLikelihoods.gvcf \
        && bgzip ~{in_sample_name}.rawLikelihoods.gvcf
    >>>
    output {
        File genotyped_gvcf = "~{in_sample_name}.rawLikelihoods.gvcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "broadinstitute/gatk@sha256:cc8981d0527e716775645b04a7f59e96a52ad59a7ae9788ddc47902384bf35aa"
    }
}

task runDeepVariantCaller {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))

        /opt/deepvariant/bin/run_deepvariant \
        --num_shards ~{in_vgcall_cores} \
        --model_type WGS \
        --regions ${CONTIG_ID} \
        --ref ~{in_reference_file} \
        --sample_name ~{in_sample_name} \
        --reads input_bam_file.bam \
        --output_vcf ~{in_sample_name}.vcf \
        && bgzip ~{in_sample_name}.vcf
    >>>
    output {
        File genotyped_vcf = "~{in_sample_name}.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:1.0.0"
    }
}

task runDeepVariantCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))

        /opt/deepvariant/bin/run_deepvariant \
        --num_shards ~{in_vgcall_cores} \
        --model_type WGS \
        --regions ${CONTIG_ID} \
        --ref ~{in_reference_file} \
        --sample_name ~{in_sample_name} \
        --reads input_bam_file.bam \
        --output_gvcf ~{in_sample_name}.rawLikelihoods.gvcf \
        --output_vcf ~{in_sample_name}.rawLikelihoods.vcf \
        && bgzip ~{in_sample_name}.rawLikelihoods.gvcf
    >>>
    output {
        File genotyped_vcf = "~{in_sample_name}.rawLikelihoods.vcf.gz"
        File genotyped_gvcf = "~{in_sample_name}.rawLikelihoods.gvcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:1.0.0"
    }
}

task chunkAlignmentsByPathNames {
    input {
        String in_sample_name
        File in_merged_sorted_gam
        File in_merged_sorted_gam_gai
        File in_xg_file
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
        VG_CALL_SV_OPTIONS=""
        

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

        vg index -x ~{chunk_tag}.xg -t ~{in_vgcall_cores} ~{in_vg_file} && \
        vg filter ~{in_gam_file} -t 1 -q 15 -r 0.9 -fu -s 1000 -D 999 -x ~{chunk_tag}.xg > ~{chunk_tag}.filtered.gam

        vg pack \
            -x ~{chunk_tag}.xg \
            -t ~{in_vgcall_cores} \
            -g ~{chunk_tag}.filtered.gam -Q 5 -o ~{chunk_tag}.pack && \
            
        vg call \
            ~{chunk_tag}.xg \
            -t ~{in_vgcall_cores} \
            -k ~{chunk_tag}.pack \
            -p ${PATH_NAME} \
            -o ${OFFSET} \
            -l ${CHR_LENGTH} \
            ${VG_CALL_SV_OPTIONS} > ~{chunk_tag}.vcf && \
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
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        time: 60
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        time: 30
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: in_vg_container
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        bcftools norm -m-both --threads ~{in_vgcall_cores} -o ~{in_sample_name}.unrolled.vcf ~{in_bgzip_vcf_file}
    >>>
    output {
        File output_normalized_vcf = "~{in_sample_name}.unrolled.vcf"
    }
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File? in_snpeff_database
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:4.3.1t--2"
    }
}


