version 1.0

### vg_deeptrio_calling_workflow.wdl ###
## Author: Charles Markello
## Description: Core VG DEEPTRIO workflow for calling variants on trio datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgDeepTrioCall {
    input {
        Boolean GOOGLE_CLEANUP_MODE = false                         # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup
        Boolean ABRA_REALIGN = false                                # Set to 'true' to use ABRA2 IndelRealigner instead of GATK for indel realignment
        File? MATERNAL_BAM_FILE                                     # Input maternal surjected .bam file
        File? MATERNAL_BAM_FILE_INDEX                               # Input maternal .bai index of surjected .bam file.
        Array[File]? MATERNAL_BAM_CONTIG_LIST                       # Input maternal bam per contig in a list of files.
        Array[File]? MATERNAL_BAM_INDEX_CONTIG_LIST                 # Input maternal bam index per contig in a list of files following the same contig order as MATERNAL_BAM_CONTIG_LIST.
        File? PATERNAL_BAM_FILE                                     # Input paternal surjected .bam file
        File? PATERNAL_BAM_FILE_INDEX                               # Input paternal .bai index of surjected .bam file.
        Array[File]? PATERNAL_BAM_CONTIG_LIST                       # Input paternal bam per contig in a list of files.
        Array[File]? PATERNAL_BAM_INDEX_CONTIG_LIST                 # Input paternal bam index per contig in a list of files following the same contig order as PATERNAL_BAM_CONTIG_LIST.
        File? CHILD_BAM_FILE                                         # Input child surjected .bam file
        File? CHILD_BAM_FILE_INDEX                                   # Input child .bai index of surjected .bam file.
        String SAMPLE_NAME_CHILD                                    # The child sample name
        String SAMPLE_NAME_MATERNAL                                 # The maternal sample name
        String SAMPLE_NAME_PATERNAL                                 # The paternal sample name
        File PATH_LIST_FILE                                         # Text file where each line is a path name in the XG index used in 
        File REF_FILE                                               # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                         # Path to .fai index of the REF_FILE fasta reference
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.19.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int VGCALL_CORES = 32
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
    }
    # Split merged alignment by contigs list
    call splitBAMbyPath as splitChildBAMbyPath {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_merged_bam_file=CHILD_BAM_FILE,
            in_merged_bam_file_index=CHILD_BAM_FILE_INDEX,
            in_path_list_file=PATH_LIST_FILE,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }
    if (!defined(MATERNAL_BAM_CONTIG_LIST)) {
        call splitBAMbyPath as splitMaternalBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_merged_bam_file=MATERNAL_BAM_FILE,
                in_merged_bam_file_index=MATERNAL_BAM_FILE_INDEX,
                in_path_list_file=PATH_LIST_FILE,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
    }
    if (!defined(PATERNAL_BAM_CONTIG_LIST)) {
        call splitBAMbyPath as splitPaternalBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_merged_bam_file=PATERNAL_BAM_FILE,
                in_merged_bam_file_index=PATERNAL_BAM_FILE_INDEX,
                in_path_list_file=PATH_LIST_FILE,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
    }
    # Run distributed DeepTRIO linear variant calling for each chromosomal contig
    Array[File] maternal_bams_by_contig = select_first([MATERNAL_BAM_CONTIG_LIST, splitMaternalBAMbyPath.bam_contig_files])
    Array[File] maternal_bam_indexes_by_contig = select_first([MATERNAL_BAM_INDEX_CONTIG_LIST, splitMaternalBAMbyPath.bam_contig_files_index])
    Array[Pair[File, File]] maternal_bams_and_indexes_by_contig = zip(maternal_bams_by_contig, maternal_bam_indexes_by_contig)
    Array[File] paternal_bams_by_contig = select_first([PATERNAL_BAM_CONTIG_LIST, splitPaternalBAMbyPath.bam_contig_files])
    Array[File] paternal_bam_indexes_by_contig = select_first([PATERNAL_BAM_INDEX_CONTIG_LIST, splitPaternalBAMbyPath.bam_contig_files_index])
    Array[Pair[File, File]] paternal_bams_and_indexes_by_contig = zip(paternal_bams_by_contig, paternal_bam_indexes_by_contig)
    Array[Pair[File, File]] child_bams_and_indexes_by_contig = zip(splitChildBAMbyPath.bam_contig_files, splitChildBAMbyPath.bam_contig_files_index)
    scatter (deeptrio_caller_input_files in zip(child_bams_and_indexes_by_contig, zip(maternal_bams_and_indexes_by_contig, paternal_bams_and_indexes_by_contig))) {
        Pair[File,File] child_deeptrio_caller_input_files = deeptrio_caller_input_files.left
        Pair[File,File] maternal_deeptrio_caller_input_files = deeptrio_caller_input_files.right.left
        Pair[File,File] paternal_deeptrio_caller_input_files = deeptrio_caller_input_files.right.right
        if (ABRA_REALIGN) {
            call runGATKRealignerTargetCreator as realignTargetCreateChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_deeptrio_caller_input_files,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE
            }
            call runAbraRealigner as abraRealignChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_deeptrio_caller_input_files,
                    in_target_bed_file=realignTargetCreateChild.realigner_target_bed,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE
            }
            # Only do indel realignment of parents if there isn't already a supplied indel realigned bam set
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKRealignerTargetCreator as realignTargetCreateMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_deeptrio_caller_input_files,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
                call runAbraRealigner as abraRealignMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_deeptrio_caller_input_files,
                        in_target_bed_file=realignTargetCreateMaternal.realigner_target_bed,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKRealignerTargetCreator as realignTargetCreatePaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_deeptrio_caller_input_files,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
                call runAbraRealigner as abraRealignPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_deeptrio_caller_input_files,
                        in_target_bed_file=realignTargetCreatePaternal.realigner_target_bed,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
        }
        if (!ABRA_REALIGN) {
            call runGATKIndelRealigner as gatkRealignChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_deeptrio_caller_input_files,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE
            }
            # Only do indel realignment of parents if there isn't already a supplied indel realigned bam set
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKIndelRealigner as gatkRealignMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_deeptrio_caller_input_files,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKIndelRealigner as gatkRealignPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_deeptrio_caller_input_files,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
        }
        File child_indel_realigned_bam = select_first([abraRealignChild.indel_realigned_bam, gatkRealignChild.indel_realigned_bam])
        File child_indel_realigned_bam_index = select_first([abraRealignChild.indel_realigned_bam_index, gatkRealignChild.indel_realigned_bam_index])
        File maternal_indel_realigned_bam = select_first([abraRealignMaternal.indel_realigned_bam, gatkRealignMaternal.indel_realigned_bam, maternal_deeptrio_caller_input_files.left])
        File maternal_indel_realigned_bam_index = select_first([abraRealignMaternal.indel_realigned_bam_index, gatkRealignMaternal.indel_realigned_bam_index, maternal_deeptrio_caller_input_files.right])
        File paternal_indel_realigned_bam = select_first([abraRealignPaternal.indel_realigned_bam, gatkRealignPaternal.indel_realigned_bam, paternal_deeptrio_caller_input_files.left])
        File paternal_indel_realigned_bam_index = select_first([abraRealignPaternal.indel_realigned_bam_index, gatkRealignPaternal.indel_realigned_bam_index, paternal_deeptrio_caller_input_files.right])

        call runDeepTrioMakeExamples {
            input:
                in_child_name=SAMPLE_NAME_CHILD,
                in_maternal_name=SAMPLE_NAME_MATERNAL,
                in_paternal_name=SAMPLE_NAME_PATERNAL,
                in_child_bam_file=child_indel_realigned_bam,
                in_child_bam_file_index=child_indel_realigned_bam_index,
                in_maternal_bam_file=maternal_indel_realigned_bam,
                in_maternal_bam_file_index=maternal_indel_realigned_bam_index,
                in_paternal_bam_file=paternal_indel_realigned_bam,
                in_paternal_bam_file_index=paternal_indel_realigned_bam_index,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call runDeepTrioCallVariants as callVariantsChild {
            input:
                in_sample_name=SAMPLE_NAME_CHILD,
                in_sample_type="child",
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_examples_file=runDeepTrioMakeExamples.child_examples_file,
                in_nonvariant_site_tf_file=runDeepTrioMakeExamples.child_nonvariant_site_tf_file,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
            call runDeepTrioCallVariants as callVariantsMaternal {
                input:
                    in_sample_name=SAMPLE_NAME_MATERNAL,
                    in_sample_type="parent1",
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_examples_file=runDeepTrioMakeExamples.maternal_examples_file,
                    in_nonvariant_site_tf_file=runDeepTrioMakeExamples.maternal_nonvariant_site_tf_file,
                    in_vgcall_cores=VGCALL_CORES,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM
            }
        }
        if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
            call runDeepTrioCallVariants as callVariantsPaternal {
                input:
                    in_sample_name=SAMPLE_NAME_PATERNAL,
                    in_sample_type="parent2",
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_examples_file=runDeepTrioMakeExamples.paternal_examples_file,
                    in_nonvariant_site_tf_file=runDeepTrioMakeExamples.paternal_nonvariant_site_tf_file,
                    in_vgcall_cores=VGCALL_CORES,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM
            }
        }
    }
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks as concatVCFChunksChild {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_clipped_vcf_chunk_files=callVariantsChild.output_gvcf_file,
            in_vgcall_disk=VGCALL_DISK,
            in_vgcall_mem=VGCALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF as bgzipVGCalledChildVCF {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_merged_vcf_file=concatVCFChunksChild.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_vgcall_disk=VGCALL_DISK,
            in_vgcall_mem=VGCALL_MEM
    }
    if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] maternal_contig_gvcf_output_list = select_all(callVariantsMaternal.output_gvcf_file)
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksMaternal {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_clipped_vcf_chunk_files=maternal_contig_gvcf_output_list,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledMaternalVCF {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_merged_vcf_file=concatVCFChunksMaternal.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] paternal_contig_gvcf_output_list = select_all(callVariantsPaternal.output_gvcf_file)
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksPaternal {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_clipped_vcf_chunk_files=paternal_contig_gvcf_output_list,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledPaternalVCF {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_merged_vcf_file=concatVCFChunksPaternal.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    
    # Collect parental indel-realigned BAM contig lists
    if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] maternal_indelrealigned_bams_by_contig = maternal_indel_realigned_bam
        Array[File] maternal_indelrealigned_bam_indexes_by_contig = maternal_indel_realigned_bam_index
    }
    if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] paternal_indelrealigned_bams_by_contig = paternal_indel_realigned_bam
        Array[File] paternal_indelrealigned_bam_indexes_by_contig = paternal_indel_realigned_bam_index
    }
    
    output {
        File output_child_gvcf = bgzipVGCalledChildVCF.output_merged_vcf
        File output_child_gvcf_index = bgzipVGCalledChildVCF.output_merged_vcf_index
        File? output_maternal_gvcf = bgzipVGCalledMaternalVCF.output_merged_vcf
        File? output_maternal_gvcf_index = bgzipVGCalledMaternalVCF.output_merged_vcf_index
        File? output_paternal_gvcf = bgzipVGCalledPaternalVCF.output_merged_vcf
        File? output_paternal_gvcf_index = bgzipVGCalledPaternalVCF.output_merged_vcf_index
        Array[File] output_child_indelrealigned_bams = child_indel_realigned_bam
        Array[File] output_child_indelrealigned_bam_indexes = child_indel_realigned_bam_index
        Array[File]? output_maternal_indelrealigned_bams = maternal_indelrealigned_bams_by_contig
        Array[File]? output_maternal_indelrealigned_bam_indexes = maternal_indelrealigned_bam_indexes_by_contig
        Array[File]? output_paternal_indelrealigned_bams = paternal_indelrealigned_bams_by_contig
        Array[File]? output_paternal_indelrealigned_bam_indexes = paternal_indelrealigned_bam_indexes_by_contig
    }
}
########################
### TASK DEFINITIONS ###
########################
task splitBAMbyPath {
    input {
        String in_sample_name
        File? in_merged_bam_file
        File? in_merged_bam_file_index
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
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKRealignerTargetCreator { 
    input { 
        String in_sample_name 
        Pair[File, File] in_bam_file 
        File in_reference_file 
        File in_reference_index_file 
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
 
        ln -f -s ~{in_bam_file.left} input_bam_file.bam 
        ln -f -s ~{in_bam_file.right} input_bam_file.bam.bai 
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}')) 
         
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \ 
          --remove_program_records \ 
          -drf DuplicateRead \ 
          --disable_bam_indexing \ 
          -nt "32" \ 
          -R ~{in_reference_file} \ 
          -L ${CONTIG_ID} \ 
          -I input_bam_file.bam \ 
          --out forIndelRealigner.intervals 
         
        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG}.intervals.bed 
    >>> 
    output { 
        File realigner_target_bed = glob("*.bed")[0] 
    } 
    runtime { 
        time: 180 
        memory: 20 + " GB" 
        cpu: 32 
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b" 
    } 
}

task runAbraRealigner {
    input {
        String in_sample_name
        Pair[File, File] in_bam_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
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

        ln -f -s ~{in_bam_file.left} input_bam_file.bam
        ln -f -s ~{in_bam_file.right} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        java -Xmx20G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref ~{in_reference_file} \
          --index \
          --threads 32
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bam.bai")[0]
    }
    runtime {
        time: 180
        memory: 20 + " GB"
        cpu: 32
        docker: "dceoy/abra2:latest"
    }
}

task runGATKIndelRealigner {
    input {
        String in_sample_name
        Pair[File, File] in_bam_file
        File in_reference_file
        File in_reference_index_file
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

        ln -f -s ~{in_bam_file.left} input_bam_file.bam
        ln -f -s ~{in_bam_file.right} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "$(nproc --all)" \
          -R ~{in_reference_file} \
          -I input_bam_file.bam \
          -L ${CONTIG_ID} \
          --out forIndelRealigner.intervals \
        && java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
          --remove_program_records \
          -R ~{in_reference_file} \
          --targetIntervals forIndelRealigner.intervals \
          -I input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bam.bai")[0]
    }
    runtime {
        time: 180
        memory: 20 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}

task runDeepTrioMakeExamples {
    input {
        String in_child_name
        String in_maternal_name
        String in_paternal_name
        File in_child_bam_file
        File in_child_bam_file_index
        File in_maternal_bam_file
        File in_maternal_bam_file_index
        File in_paternal_bam_file
        File in_paternal_bam_file_index
        File in_reference_file
        File in_reference_index_file
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

        ln -s ~{in_child_bam_file} input_bam_file.child.bam
        ln -s ~{in_child_bam_file_index} input_bam_file.child.bam.bai
        ln -s ~{in_maternal_bam_file} input_bam_file.maternal.bam
        ln -s ~{in_maternal_bam_file_index} input_bam_file.maternal.bam.bai
        ln -s ~{in_paternal_bam_file} input_bam_file.paternal.bam
        ln -s ~{in_paternal_bam_file_index} input_bam_file.paternal.bam.bai
        CONTIG_ID=($(ls ~{in_child_bam_file} | awk -F'.' '{print $(NF-1)}'))
        
        seq 0 ~{in_vgcall_cores} | \
        parallel -q halt 2 --line-buffer /opt/deepvariant/bin/deeptrio/make_examples \
        --mode calling \
        --ref ~{in_reference_file} \
        --reads_parent1 input_bam_file.paternal.bam \
        --reads_parent2 input_bam_file.maternal.bam \
        --reads input_bam_file.child.bam \
        --examples ./make_examples.tfrecord@~{in_vgcall_cores}.gz \
        --sample_name in_child_name \
        --sample_name_parent1 in_paternal_name \
        --sample_name_parent2 in_maternal_name \
        --gvcf ./gvcf.tfrecord@~{in_vgcall_cores}.gz \
        --min_mapping_quality 1 \
        --pileup_image_height_child 60 \
        --pileup_image_height_parent 40 \
        --regions ${CONTIG_ID} \
        --task {}

        ls | grep 'make_examples_child.tfrecord' | tar -czf 'make_examples_child.tfrecord.tar.gz' -T -
        ls | grep 'make_examples_parent1.tfrecord' | tar -czf 'make_examples_parent1.tfrecord.tar.gz' -T -
        ls | grep 'make_examples_parent2.tfrecord' | tar -czf 'make_examples_parent2.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_child.tfrecord' | tar -czf 'gvcf_child.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_parent1.tfrecord' | tar -czf 'gvcf_parent1.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_parent2.tfrecord' | tar -czf 'gvcf_parent2.tfrecord.tar.gz' -T -
    >>>
    output {
        File child_examples_file = "make_examples_child.tfrecord.tar.gz"
        File paternal_examples_file = "make_examples_parent1.tfrecord.tar.gz"
        File maternal_examples_file = "make_examples_parent2.tfrecord.tar.gz"
        File child_nonvariant_site_tf_file = "gvcf_child.tfrecord.tar.gz"
        File paternal_nonvariant_site_tf_file = "gvcf_parent1.tfrecord.tar.gz"
        File maternal_nonvariant_site_tf_file = "gvcf_parent2.tfrecord.tar.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:deeptrio-1.1.0"
    }
}

task runDeepTrioCallVariants {
    input {
        String in_sample_name
        String in_sample_type
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
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
        tar -xzf ~{in_examples_file}
        tar -xzf ~{in_nonvariant_site_tf_file}
        if [ ~{in_sample_type} == "child" ]; then
            EXAMPLES_FILE="make_examples_child.tfrecord~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_child.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_child.tfrecord~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/child/model.ckpt"
        elif [ ~{in_sample_type} == "parent1" ]; then
            EXAMPLES_FILE="make_examples_parent1.tfrecord~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_parent1.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_parent1.tfrecord~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/parent/model.ckpt"
        elif [ ~{in_sample_type} == "parent2" ]; then
            EXAMPLES_FILE="make_examples_parent2.tfrecord~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_parent2.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_parent2.tfrecord~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/parent/model.ckpt"
        fi
        /opt/deepvariant/bin/call_variants \
        --outfile ${OUTPUT_FILE} \
        --examples ${EXAMPLES_FILE} \
        --checkpoint ${DEEPTRIO_MODEL} && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref ~{in_reference_file} \
        --infile ${OUTPUT_FILE} \
        --nonvariant_site_tfrecord_path ${NONVARIANT_SITE_FILE} \
        --outfile "~{in_sample_name}_deeptrio.vcf.gz" \
        --gvcf_outfile "~{in_sample_name}_deeptrio.g.vcf.gz"
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deeptrio.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deeptrio.g.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:deeptrio-1.1.0"
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
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
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


