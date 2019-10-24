version 1.0

### vg_multi_map.wdl ###
## Author: Charles Markello
## Description: Core VG mapping workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMapCall {
    input {
        Boolean VGMPMAP_MODE = true                     # Set to 'false' to use "VG MAP" or set to 'true' to use "VG MPMAP" algorithm
        Boolean SURJECT_MODE = true                     # Set to 'true' to run pipeline using alignmed BAM files surjected from GAM. Set to 'false' to output graph aligned GAM files.
        Boolean GOOGLE_CLEANUP_MODE = false             # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        File INPUT_READ_FILE_1                          # Input sample 1st read pair fastq.gz
        File INPUT_READ_FILE_2                          # Input sample 2nd read pair fastq.gz
        String SAMPLE_NAME                              # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                    # Path to .xg index file
        File GCSA_FILE                                  # Path to .gcsa index file
        File GCSA_LCP_FILE                              # Path to .gcsa.lcp index file
        File? GBWT_FILE                                 # (OPTIONAL) Path to .gbwt index file
        File? SNARLS_FILE                               # (OPTIONAL) Path to .snarls index file
        File REF_FILE                                   # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                             # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                              # Path to .dict file of the REF_FILE fasta reference
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAM_CORES = 56
        Int MERGE_GAM_DISK = 400
        Int MERGE_GAM_MEM = 100
        Int MERGE_GAM_TIME = 1200
    }

    # Split input reads into chunks for parallelized mapping
    call splitReads as firstReadPair {
        input:
            in_read_file=INPUT_READ_FILE_1,
            in_pair_id="1",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }
    call splitReads as secondReadPair {
        input:
            in_read_file=INPUT_READ_FILE_2,
            in_pair_id="2",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
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

    ################################################################
    # Distribute vg mapping opperation over each chunked read pair #
    ################################################################
    Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
    scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
        if (VGMPMAP_MODE) {
            call runVGMPMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_snarls_file=SNARLS_FILE,
                    in_sample_name=SAMPLE_NAME,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        if (!VGMPMAP_MODE) {
            call runVGMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        
        # Surject GAM alignment files to BAM if SURJECT_MODE set to true
        File vg_map_algorithm_chunk_gam_output = select_first([runVGMAP.chunk_gam_file, runVGMPMAP.chunk_gam_file])
        # Cleanup input reads after use
        if (GOOGLE_CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpVGMapperInputsGoogle {
                input:
                    previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                    current_task_output = vg_map_algorithm_chunk_gam_output
            }
        }
        if (!GOOGLE_CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpVGMapperInputsUnix {
                input:
                    previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                    current_task_output = vg_map_algorithm_chunk_gam_output
            }
        }
        if (SURJECT_MODE) {
            call runSurject {
                input:
                    in_gam_chunk_file=vg_map_algorithm_chunk_gam_output,
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_sample_name=SAMPLE_NAME,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
            call sortMDTagBAMFile {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_chunk_file=runSurject.chunk_bam_file,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM,
                    in_vgmpmap_mode=VGMPMAP_MODE
            }
            call runPICARD {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=sortMDTagBAMFile.sorted_bam_file,
                    in_bam_file_index=sortMDTagBAMFile.sorted_bam_file_index,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE,
                    in_map_cores=MAP_CORES,
                    in_map_mem=MAP_MEM
            }
            # Cleanup intermediate surject files after use
            if (GOOGLE_CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpVGSurjectInputsGoogle {
                    input:
                        previous_task_outputs = [vg_map_algorithm_chunk_gam_output, runSurject.chunk_bam_file, sortMDTagBAMFile.sorted_bam_file, sortMDTagBAMFile.sorted_bam_file_index],
                        current_task_output = runPICARD.mark_dupped_reordered_bam
                }
            }
            if (!GOOGLE_CLEANUP_MODE) {
                call cleanUpUnixFilesystem as cleanUpVGSurjectInputsUnix {
                    input:
                        previous_task_outputs = [vg_map_algorithm_chunk_gam_output, runSurject.chunk_bam_file, sortMDTagBAMFile.sorted_bam_file, sortMDTagBAMFile.sorted_bam_file_index],
                        current_task_output = runPICARD.mark_dupped_reordered_bam
                }
            }
        }
    }   

    if (SURJECT_MODE) {
        # Merge chunked alignments from surjected GAM files
        Array[File?] alignment_chunk_bam_files_maybes = runPICARD.mark_dupped_reordered_bam
        Array[File] alignment_chunk_bam_files_valid = select_all(alignment_chunk_bam_files_maybes)
        call mergeAlignmentBAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=alignment_chunk_bam_files_valid
        }
        File merged_bam_file_output = mergeAlignmentBAMChunks.merged_bam_file
        File merged_bam_file_index_output = mergeAlignmentBAMChunks.merged_bam_file_index
        if (GOOGLE_CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpAlignmentBAMChunksGoogle {
                input:
                    previous_task_outputs = alignment_chunk_bam_files_valid,
                    current_task_output = merged_bam_file_output
            }
        }
        if (!GOOGLE_CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpAlignmentBAMChunksUnix {
                input:
                    previous_task_outputs = alignment_chunk_bam_files_valid,
                    current_task_output = merged_bam_file_output
            }
        }
    } 
    if (!SURJECT_MODE) {
        call mergeAlignmentGAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_vg_container=VG_CONTAINER,
                in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
                in_merge_gam_cores=MERGE_GAM_CORES,
                in_merge_gam_disk=MERGE_GAM_DISK,
                in_merge_gam_mem=MERGE_GAM_MEM,
                in_merge_gam_time=MERGE_GAM_TIME
        }
        # Cleanup gam chunk files after use
        if (GOOGLE_CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpGAMChunksGoogle {
                input:
                    previous_task_outputs = vg_map_algorithm_chunk_gam_output,
                    current_task_output = mergeAlignmentGAMChunks.merged_sorted_gam_file
            }
        }
        if (!GOOGLE_CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpGAMChunksUnix {
                input:
                    previous_task_outputs = vg_map_algorithm_chunk_gam_output,
                    current_task_output = mergeAlignmentGAMChunks.merged_sorted_gam_file
            }
        }
    }
    output {
        File? output_bam = merged_bam_file_output
        File? output_bam_index = merged_bam_file_index_output
        File? output_gam = mergeAlignmentGAMChunks.merged_sorted_gam_file
        File? output_gam_index = mergeAlignmentGAMChunks.merged_sorted_gam_gai_file
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

task splitReads {
    input {
        File in_read_file
        String in_pair_id
        String in_vg_container
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk
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

        CHUNK_LINES=$(( ~{in_reads_per_chunk} * 4 ))
        gzip -cd ~{in_read_file} | split -l $CHUNK_LINES --filter='pigz -p ~{in_split_read_cores} > ${FILE}.fq.gz' - "fq_chunk_~{in_pair_id}.part."
    >>>
    output {
        Array[File] output_read_chunks = glob("fq_chunk_~{in_pair_id}.part.*")
    }
    runtime {
        cpu: in_split_read_cores
        memory: "40 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: in_vg_container
    }
}

task extractPathNames {
    input {
        File in_xg_file
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

task runVGMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gcsa_file
        File in_gcsa_lcp_file
        File? in_gbwt_file
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    Boolean gbwt_options = defined(in_gbwt_file)

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

        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        GBWT_OPTION_STRING=""
        if [ ~{gbwt_options} == true ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file}"
        fi
        ln -s ~{in_gcsa_file} input_gcsa_file.gcsa
        ln -s ~{in_gcsa_lcp_file} input_gcsa_file.gcsa.lcp
        vg map \
          -x ~{in_xg_file} \
          -g input_gcsa_file.gcsa \
          ${GBWT_OPTION_STRING} \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task runVGMPMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gcsa_file
        File in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_snarls_file
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    Boolean gbwt_options = defined(in_gbwt_file)
    Boolean snarl_options = defined(in_snarls_file)

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

        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        GBWT_OPTION_STRING=""
        if [ ~{gbwt_options} == true ] && [ ~{snarl_options} == false ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file} --recombination-penalty 5.0"
        elif [ ~{gbwt_options} == true ] && [ ~{snarl_options} == true ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file} -s ~{in_snarls_file} --recombination-penalty 5.0"
        fi
        ln -s ~{in_gcsa_file} input_gcsa_file.gcsa
        ln -s ~{in_gcsa_lcp_file} input_gcsa_file.gcsa.lcp
        vg mpmap \
          -S \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -g input_gcsa_file.gcsa \
          ${GBWT_OPTION_STRING} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --sample "~{in_sample_name}" \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task runSurject {
    input {
        File in_gam_chunk_file
        File in_xg_file
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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

        READ_CHUNK_ID=($(ls ~{in_gam_chunk_file} | awk -F'.' '{print $(NF-1)}'))
        vg surject \
          -i \
          -x ~{in_xg_file} \
          -t ~{in_map_cores} \
          -b ~{in_gam_chunk_file} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
    >>>
    output {
        File chunk_bam_file = glob("*.bam")[0]
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
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
        Boolean in_vgmpmap_mode
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
        if [ ~{in_vgmpmap_mode} == false ]; then
            samtools sort \
              --threads ~{in_map_cores} \
              -n \
              ~{in_bam_chunk_file} \
              -O BAM \
            | samtools fixmate \
              -O BAM \
              - \
              - \
            | samtools sort \
              --threads ~{in_map_cores} \
              - \
              -O BAM \
            | samtools addreplacerg \
              -O BAM \
              -r ID:1 -r LB:lib1 -r SM:~{in_sample_name} -r PL:illumina -r PU:unit1 \
              - \
            | samtools view \
              -@ 32 \
              -h -O SAM \
              - \
            | samtools view \
              -@ 32 \
              -h -O BAM \
              - \
            | samtools calmd \
              -b \
              - \
              ~{in_reference_file} \
              > ~{in_sample_name}_positionsorted.mdtag.bam \
            && samtools index \
              ~{in_sample_name}_positionsorted.mdtag.bam
        else
            samtools sort \
              --threads ~{in_map_cores} \
              ~{in_bam_chunk_file} \
              -O BAM \
            | samtools calmd \
              -b \
              - \
              ~{in_reference_file} \
              > ~{in_sample_name}_positionsorted.mdtag.bam \
            && samtools index \
              ~{in_sample_name}_positionsorted.mdtag.bam
        fi
    >>>
    output {
        File sorted_bam_file = "${in_sample_name}_positionsorted.mdtag.bam"
        File sorted_bam_file_index = "${in_sample_name}_positionsorted.mdtag.bam.bai"
    }
    runtime {
        time: 200
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task runPICARD {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_map_cores
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

        java -Xmx${in_map_mem}g -XX:ParallelGCThreads=${in_map_cores} -jar /usr/picard/picard.jar MarkDuplicates \
          PROGRAM_RECORD_ID=null \
          VALIDATION_STRINGENCY=LENIENT \
          I=${in_bam_file} \
          O=${in_sample_name}.mdtag.dupmarked.bam \
          M=marked_dup_metrics.txt 2> mark_dup_stderr.txt \
        && java -Xmx${in_map_mem}g -XX:ParallelGCThreads=${in_map_cores} -jar /usr/picard/picard.jar ReorderSam \
            VALIDATION_STRINGENCY=LENIENT \
            REFERENCE=${in_reference_file} \
            INPUT=${in_sample_name}.mdtag.dupmarked.bam \
            OUTPUT=${in_sample_name}.mdtag.dupmarked.reordered.bam \
        && rm -f ${in_sample_name}.mdtag.dupmarked.bam
    }
    output {
        File mark_dupped_reordered_bam = "${in_sample_name}.mdtag.dupmarked.reordered.bam"
    }
    runtime {
        time: 180
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        docker: "broadinstitute/picard:2.20.4"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
        samtools merge \
          -f -p -c --threads 32 \
          ${in_sample_name}_merged.positionsorted.bam \
          ${sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ${in_sample_name}_merged.positionsorted.bam
    }
    output {
        File merged_bam_file = "${in_sample_name}_merged.positionsorted.bam"
        File merged_bam_file_index = "${in_sample_name}_merged.positionsorted.bam.bai"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

# TODO: Duplicate task of mergeAlignmentBAMChunksVGMPMAP
task mergeIndelRealignedBAMs {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
        samtools merge \
          -f -p -c --threads 32 \
          ${in_sample_name}_merged.indel_realigned.bam \
          ${sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ${in_sample_name}_merged.indel_realigned.bam
    }
    output {
        File merged_indel_realigned_bam_file = "${in_sample_name}_merged.indel_realigned.bam"
        File merged_indel_realigned_bam_file_index = "${in_sample_name}_merged.indel_realigned.bam.bai"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task mergeAlignmentGAMChunks {
    input {
        String in_sample_name
        String in_vg_container
        Array[File] in_alignment_gam_chunk_files
        Int in_merge_gam_cores
        Int in_merge_gam_disk
        Int in_merge_gam_mem
        Int in_merge_gam_time
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

        cat ${sep=" " in_alignment_gam_chunk_files} > ${in_sample_name}_merged.gam \
        && vg gamsort \
            ${in_sample_name}_merged.gam \
            -i ${in_sample_name}_merged.sorted.gam.gai \
            -t ${in_merge_gam_cores} > ${in_sample_name}_merged.sorted.gam
    }
    output {
        File merged_sorted_gam_file = "${in_sample_name}_merged.sorted.gam"
        File merged_sorted_gam_gai_file = "${in_sample_name}_merged.sorted.gam.gai"
    }
    runtime {
        memory: in_merge_gam_mem + " GB"
        cpu: in_merge_gam_cores
        disks: "local-disk " + in_merge_gam_disk  + " SSD"
        time: in_merge_gam_time
        docker: in_vg_container
    }
}




