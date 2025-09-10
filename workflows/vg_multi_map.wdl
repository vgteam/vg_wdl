version 1.0

### vg_multi_map.wdl ###
## Author: Charles Markello
## Description: Core VG mapping workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMap {
    input {
        String MAPPER = "GIRAFFE"                       # Set to 'MAP' to use the "VG MAP" algorithm, set to 'MPMAP' to use "VG MPMAP" algorithm, set to 'GIRAFFE' to use "VG GIRAFFE".
        Boolean SURJECT_MODE = true                     # Set to 'true' to run pipeline using alignmed BAM files surjected from GAM. Set to 'false' to output graph aligned GAM files.
        Boolean CLEANUP_FILES = false                   # Set to 'false' to turn off intermediate file cleanup.
        Boolean GOOGLE_CLEANUP_MODE = true              # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        Boolean SMALL_RESOURCES = false                 # Set to 'true' to use small resources for tiny test dataset
        File INPUT_READ_FILE_1                          # Input sample 1st read pair fastq.gz
        File INPUT_READ_FILE_2                          # Input sample 2nd read pair fastq.gz
        String SAMPLE_NAME                              # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.64.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.64.0)
        Int READS_PER_CHUNK = 200000000                 # Number of reads contained in each mapping chunk (20000000 for wgs)
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                    # Path to .xg index file
        File? GCSA_FILE                                 # (OPTIONAL) Path to .gcsa index file
        File? GCSA_LCP_FILE                             # (OPTIONAL) Path to .gcsa.lcp index file
        File? GBWT_FILE                                 # (OPTIONAL) Path to .gbwt index file
        File? GGBWT_FILE                                # (OPTIONAL) Path to .gg index file
        File? DIST_FILE                                 # (OPTIONAL) Path to .dist index file
        File? MIN_FILE                                  # (OPTIONAL) Path to .min index file
        File? SNARLS_FILE                               # (OPTIONAL) Path to .snarls index file
        File REF_FILE                                   # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                             # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                              # Path to .dict file of the REF_FILE fasta reference
    }

    # Split input reads into chunks for parallelized mapping
    call splitReads as firstReadPair {
        input:
            in_read_file=INPUT_READ_FILE_1,
            in_pair_id="1",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_small_resources=SMALL_RESOURCES
    }
    call splitReads as secondReadPair {
        input:
            in_read_file=INPUT_READ_FILE_2,
            in_pair_id="2",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_small_resources=SMALL_RESOURCES
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
        if (MAPPER == "MPMAP") {
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
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_small_resources=SMALL_RESOURCES
            }
        }
        if (MAPPER == "MAP") {
            call runVGMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_small_resources=SMALL_RESOURCES
            }
        }
        if (MAPPER == "GIRAFFE") {
            call runVGGIRAFFE {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_ggbwt_file=GGBWT_FILE,
                    in_dist_file=DIST_FILE,
                    in_min_file=MIN_FILE,
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_small_resources=SMALL_RESOURCES
            }
        }
        
        # Surject GAM alignment files to BAM if SURJECT_MODE set to true
        File vg_map_algorithm_chunk_gam_output = select_first([runVGMAP.chunk_gam_file, runVGMPMAP.chunk_gam_file, runVGGIRAFFE.chunk_gam_file])
        # Cleanup input reads after use
        if (CLEANUP_FILES) {
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
        }
        if (SURJECT_MODE) {
            call sortMDTagBAMFile {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_chunk_file=vg_map_algorithm_chunk_gam_output,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_mapper_used=MAPPER,
                    in_small_resources=SMALL_RESOURCES
            }
            # Cleanup intermediate surject files after use
            if (CLEANUP_FILES) {
                if (GOOGLE_CLEANUP_MODE) {
                    call cleanUpGoogleFilestore as cleanUpVGSurjectInputsGoogle {
                        input:
                            previous_task_outputs = [vg_map_algorithm_chunk_gam_output],
                            current_task_output = sortMDTagBAMFile.mark_dupped_reordered_bam
                    }
                }
                if (!GOOGLE_CLEANUP_MODE) {
                    call cleanUpUnixFilesystem as cleanUpVGSurjectInputsUnix {
                        input:
                            previous_task_outputs = [vg_map_algorithm_chunk_gam_output],
                            current_task_output = sortMDTagBAMFile.mark_dupped_reordered_bam
                    }
                }
            }
        }
    }   

    if (SURJECT_MODE) {
        # Merge chunked alignments from surjected GAM files
        Array[File?] alignment_chunk_bam_files_maybes = sortMDTagBAMFile.mark_dupped_reordered_bam
        Array[File] alignment_chunk_bam_files_valid = select_all(alignment_chunk_bam_files_maybes)
        call mergeAlignmentBAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=alignment_chunk_bam_files_valid,
                in_small_resources=SMALL_RESOURCES
        }
        File merged_bam_file_output = mergeAlignmentBAMChunks.merged_bam_file
        File merged_bam_file_index_output = mergeAlignmentBAMChunks.merged_bam_file_index
        if (CLEANUP_FILES) {
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
    } 
    if (!SURJECT_MODE) {
        call mergeAlignmentGAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_vg_container=VG_CONTAINER,
                in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
                in_small_resources=SMALL_RESOURCES
        }
        # Cleanup gam chunk files after use
        if (CLEANUP_FILES) {
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
        docker: "ubuntu:latest"
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
        Boolean in_small_resources
    }
    
    Int in_split_read_cores = if in_small_resources then 2 else 32
    Int in_split_read_disk = if in_small_resources then 1 else 200
    
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
        preemptible: 2
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: "quay.io/glennhickey/pigz:2.3.1"
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
        preemptible: 2
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
        File? in_gcsa_file
        File? in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Boolean in_small_resources
    }

    Boolean gbwt_options = defined(in_gbwt_file)
    
    Int in_map_cores = if in_small_resources then 8 else 32
    Int in_map_disk = if in_small_resources then 80 else 100
    String in_map_mem = if in_small_resources then "80" else "100"

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
        if [ ~{surject_output} == false ]; then
            vg map \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else
           vg map \
              --ref-paths ~{in_ref_dict} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --surject-to bam \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        preemptible: 2
        time: 300
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
        File? in_gcsa_file
        File? in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_snarls_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Boolean in_small_resources
    }

    Boolean gbwt_options = defined(in_gbwt_file)
    Boolean snarl_options = defined(in_snarls_file)
    Int in_map_cores = if in_small_resources then 8 else 32
    Int in_map_disk = if in_small_resources then 80 else 100
    String in_map_mem = if in_small_resources then "80" else "100"

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
        if [ ~{surject_output} == false ]; then
            vg mpmap \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -S \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else
            vg mpmap \
              --ref-paths ~{in_ref_dict} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --output-fmt BAM \
              -S \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task runVGGIRAFFE {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File? in_gbwt_file
        File? in_ggbwt_file
        File? in_dist_file
        File? in_min_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Boolean in_small_resources
    }
    
    Int in_map_cores = if in_small_resources then 4 else 16
    Int in_map_disk = if in_small_resources then 80 else 100
    String in_map_mem = if in_small_resources then "80" else "100"

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
        if [ ~{surject_output} == false ]; then
            vg giraffe \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -H ~{in_gbwt_file} \
              -g ~{in_ggbwt_file} \
              -d ~{in_dist_file} \
              -m ~{in_min_file} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else
            vg giraffe \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --output-format BAM \
              --ref-paths ~{in_ref_dict} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -H ~{in_gbwt_file} \
              -g ~{in_ggbwt_file} \
              -d ~{in_dist_file} \
              -m ~{in_min_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        preemptible: 2
        time: 300
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
        String in_mapper_used
        Boolean in_small_resources
    }
    
    Int in_map_cores = if in_small_resources then 8 else 8
    Int in_map_disk = if in_small_resources then 10 else 50
    String in_map_mem = if in_small_resources then "10" else "32"

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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        if [[ ~{in_mapper_used} == *"GIRAFFE"* ]]; then
            samtools sort \
              --threads ~{in_map_cores} \
              ~{in_bam_chunk_file} \
              -O BAM \
            | samtools calmd \
              -b \
              - \
              ref.fna \
            | samtools addreplacerg \
                - \
                -O BAM \
                -o ~{in_sample_name}.mdtag.dupmarked.bam \
                -r ID:1 \
                -r LB:lib1 \
                -r SM:~{in_sample_name} \
                -r PL:illumina \
                -r PU:unit1
        else
            samtools sort \
              --threads ~{in_map_cores} \
              ~{in_bam_chunk_file} \
              -O BAM \
            | samtools calmd \
              -b \
              - \
              ref.fna \
              > ~{in_sample_name}_positionsorted.mdtag.bam \
            && samtools index \
              ~{in_sample_name}_positionsorted.mdtag.bam
            java -Xmx~{in_map_mem}g -XX:ParallelGCThreads=~{in_map_cores} -jar /usr/picard/picard.jar MarkDuplicates \
              PROGRAM_RECORD_ID=null \
              VALIDATION_STRINGENCY=LENIENT \
              I=~{in_sample_name}_positionsorted.mdtag.bam \
              O=~{in_sample_name}.mdtag.dupmarked.bam \
              M=marked_dup_metrics.txt 2> mark_dup_stderr.txt \
            && rm -f ~{in_sample_name}_positionsorted.mdtag.bam ~{in_sample_name}_positionsorted.mdtag.bam.bai
        fi
    >>>
    output {
        File mark_dupped_reordered_bam = "~{in_sample_name}.mdtag.dupmarked.bam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: if in_mapper_used == "GIRAFFE" then "quay.io/ucsc_cgl/samtools:latest" else "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Boolean in_small_resources
    }
    
    Int in_merge_bam_cores = if in_small_resources then 4 else 12
    Int in_merge_bam_disk = if in_small_resources then 10 else 100
    String in_merge_bam_mem = if in_small_resources then "10" else "20"
    Int in_merge_bam_time = if in_small_resources then 30 else 240
    
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
          -f -p -c --threads ~{in_merge_bam_cores} \
          ~{in_sample_name}_merged.positionsorted.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.positionsorted.bam
    >>>
    output {
        File merged_bam_file = "~{in_sample_name}_merged.positionsorted.bam"
        File merged_bam_file_index = "~{in_sample_name}_merged.positionsorted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: in_merge_bam_time
        memory: in_merge_bam_mem + " GB"
        cpu: in_merge_bam_cores
        disks: "local-disk " + in_merge_bam_disk  + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task mergeAlignmentGAMChunks {
    input {
        String in_sample_name
        String in_vg_container
        Array[File] in_alignment_gam_chunk_files
        Boolean in_small_resources
    }
    
    Int in_merge_gam_cores = if in_small_resources then 1 else 56
    Int in_merge_gam_disk = if in_small_resources then 20 else 400
    String in_merge_gam_mem = if in_small_resources then "20" else "100"
    Int in_merge_gam_time = if in_small_resources then 30 else 1200

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
        preemptible: 2
        memory: in_merge_gam_mem + " GB"
        cpu: in_merge_gam_cores
        disks: "local-disk " + in_merge_gam_disk  + " SSD"
        time: in_merge_gam_time
        docker: in_vg_container
    }
}




