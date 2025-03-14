version 1.0

### vg_indel_realign.wdl ###
## Author: Charles Markello
## Description: Core VG indel realignment workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMapCall {
    input {
        Boolean GOOGLE_CLEANUP_MODE = false         # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        File? INPUT_BAM_FILE                         # Input sample surjected .bam file
        File? INPUT_BAM_FILE_INDEX                   # Input sample .bai index of surjected .bam file.
        String SAMPLE_NAME                          # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.64.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.64.0)
        File? PATH_LIST_FILE                        # (OPTIONAL) Text file where each line is a path name in the XG index
        File? XG_FILE                                # Path to .xg index file
        File REF_FILE                               # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                         # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                          # Path to .dict file of the REF_FILE fasta reference
        Int MAP_CORES = 32 
        Int MAP_DISK = 100 
        Int MAP_MEM = 100
    }
    
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    if (!defined(PATH_LIST_FILE)) {
        call extractPathNames {
            input:
                in_xg_file=XG_FILE,
                in_vg_container=VG_CONTAINER
        }
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractPathNames.output_path_list])
    
    ##############################################
    # Run the indel realignment procedure #
    ##############################################
    # Split merged alignment by contigs list
    call splitBAMbyPath {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=INPUT_BAM_FILE,
            in_merged_bam_file_index=INPUT_BAM_FILE_INDEX,
            in_path_list_file=pipeline_path_list_file
    }
    # Run distributed Indel Realignment on contig BAMs
    scatter (gatk_caller_input_files in splitBAMbyPath.bams_and_indexes_by_contig) {
        call runGATKIndelRealigner {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=gatk_caller_input_files,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        # Call md tags
        call sortMDTagBAMFile {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_chunk_file=runGATKIndelRealigner.indel_realigned_bam,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
        }
        # Cleanup intermediate variant calling files after use
        if (GOOGLE_CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpLinearCallerInputsGoogle {
                input:
                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right, runGATKIndelRealigner.indel_realigned_bam],
                    current_task_output = sortMDTagBAMFile.sorted_bam_file
            }
        }
        if (!GOOGLE_CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpLinearCallerInputsUnix {
                input:
                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right, runGATKIndelRealigner.indel_realigned_bam],
                    current_task_output = sortMDTagBAMFile.sorted_bam_file
            }
        }
    }
    # Merge Indel Realigned BAM
    Array[File?] indel_realigned_bam_files_maybes = sortMDTagBAMFile.sorted_bam_file
    Array[File] indel_realigned_bam_files = select_all(indel_realigned_bam_files_maybes)
    call mergeIndelRealignedBAMs {
        input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=indel_realigned_bam_files
    }
    # Cleanup indel realigned bam files after use
    if (GOOGLE_CLEANUP_MODE) {
        call cleanUpGoogleFilestore as cleanUpIndelRealignedBamsGoogle {
            input:
                previous_task_outputs = indel_realigned_bam_files,
                current_task_output = mergeIndelRealignedBAMs.merged_indel_realigned_bam_file
        }
    }
    if (!GOOGLE_CLEANUP_MODE) {
        call cleanUpUnixFilesystem as cleanUpIndelRealignedBamsUnix {
            input:
                previous_task_outputs = indel_realigned_bam_files,
                current_task_output = mergeIndelRealignedBAMs.merged_indel_realigned_bam_file
        }
    }

    output {
        File output_bam = mergeIndelRealignedBAMs.merged_indel_realigned_bam_file
        File output_bam_index = mergeIndelRealignedBAMs.merged_indel_realigned_bam_file_index
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
        time: 20
        memory: 2 + " GB"
        cpu: 2
        disks: "local-disk 10 SSD"
        docker: "ubuntu@sha256:2695d3e10e69cc500a16eae6d6629c803c43ab075fa5ce60813a0fc49c47e859"
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
        docker: "google/cloud-sdk&sha256:4ef6b0e969fa96f10acfd893644d100469e979f4384e5e70f58be5cb80593a8a"
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
        File? in_merged_bam_file
        File? in_merged_bam_file_index
        File in_path_list_file
    }

    command <<<
        set -eux -o pipefail
        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai
        while IFS=$'\t' read -ra path_list_line; do
            path_name="${path_list_line[0]}"
            samtools view \
              -@ "$(nproc --all)" \
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
        time: 400
        memory: 10 + " GB"
        cpu: 16
        disks: "local-disk 10 SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKIndelRealigner {
    input {
        String in_sample_name
        Pair[File, File] in_bam_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
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
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "$(nproc --all)" \
          -R ~{in_reference_file} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals \
        && java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
          --remove_program_records \
          --disable_bam_indexing \
          -R ~{in_reference_file} \
          --targetIntervals forIndelRealigner.intervals \
          -I input_bam_file.bam \
          --out ~{in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.indel_realigned.bam
    >>>
    output {
        File indel_realigned_bam = "~{in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.indel_realigned.bam"
    }
    runtime {
        time: 180
        memory: 20 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}

task mergeIndelRealignedBAMs {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
        samtools merge \
          -f -p -c --threads "$(nproc --all)" \
          ~{in_sample_name}_merged.indel_realigned.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.indel_realigned.bam
    >>>
    output {
        File merged_indel_realigned_bam_file = "~{in_sample_name}_merged.indel_realigned.bam"
        File merged_indel_realigned_bam_file_index = "~{in_sample_name}_merged.indel_realigned.bam.bai"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task sortMDTagBAMFile {
    input {
        String in_sample_name
        File in_bam_chunk_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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
        samtools sort \
          --threads ${in_map_cores} \
          ${in_bam_chunk_file} \
          -O BAM \
        | samtools calmd \
          -b \
          - \
          ${in_reference_file} \
          > ${in_sample_name}_positionsorted.mdtag.bam
    }
    output {
        File sorted_bam_file = "${in_sample_name}_positionsorted.mdtag.bam"
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

