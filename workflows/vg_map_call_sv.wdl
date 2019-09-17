version 1.0

### vg_map_call_sv.wdl ###
# Authors: Charles Markello, Jean Monlong
# Description: Core VG mapping and variant calling workflow for single sample datasets.
# Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMapCall {
    input {
        String SAMPLE_NAME                      # The sample name
        File INPUT_READ_FILE_1                  # Input sample 1st read pair fastq.gz
        File INPUT_READ_FILE_2                  # Input sample 2nd read pair fastq.gz
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        File XG_FILE                            # Path to .xg index file
        File GCSA_FILE                          # Path to .gcsa index file
        File GCSA_LCP_FILE                      # Path to .gcsa.lcp index file
        File? GBWT_FILE                         # (OPTIONAL) Path to .gbwt index file
        File? PATH_LIST_FILE                    # (OPTIONAL) Text file where each line is a path name in the XG index
        Int READS_PER_CHUNK = 20000000          # Number of reads contained in each mapping chunk (20000000 for wgs).
        Int CHUNK_BASES = 50000000              # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                      # Number of overlapping bases between each .gam chunk
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAM_CORES = 64
        Int MERGE_GAM_DISK = 400
        Int MERGE_GAM_MEM = 100
        Int MERGE_GAM_TIME = 1200
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 400
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 32
        Int VGCALL_DISK = 100
        Int VGCALL_MEM = 64
        Boolean PACK_MODE = true               # Set to 'true' to run vg call from packed files and avoids GAM sorting (SURJECT_MODE must be 'false' for this feature to be used)
        Boolean GOOGLE_CLEANUP_MODE = false     # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        Boolean CLEANUP_MODE = true             # Set to 'false' to disable cleanup (useful to be able to use use call caching to restart workflows)
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
    scatter (read_pair_chunk_files in read_pair_chunk_files_list){
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
        
        File vg_map_algorithm_chunk_gam_output = runVGMAP.chunk_gam_file
        # Cleanup input reads after use
        if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpVGMapperInputsGoogle {
                input:
                previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                current_task_output = vg_map_algorithm_chunk_gam_output
            }
        }
        if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpVGMapperInputsUnix {
                input:
                previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                current_task_output = vg_map_algorithm_chunk_gam_output
            }
        }
    }
    
    ################################
    # Run the VG calling procedure #
    ################################
    # Merge chunked graph alignments
    call mergeAlignmentGAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_vg_container=VG_CONTAINER,
        in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
        in_merge_gam_cores=MERGE_GAM_CORES,
        in_merge_gam_disk=MERGE_GAM_DISK,
        in_merge_gam_mem=MERGE_GAM_MEM,
        in_merge_gam_time=MERGE_GAM_TIME,
        sort_gam=!PACK_MODE
    }
    # Cleanup gam chunk files after use
    if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpGoogleFilestore as cleanUpGAMChunksGoogle {
            input:
            previous_task_outputs = vg_map_algorithm_chunk_gam_output,
            current_task_output = mergeAlignmentGAMChunks.merged_gam_file
        }
    }
    if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpUnixFilesystem as cleanUpGAMChunksUnix {
            input:
            previous_task_outputs = vg_map_algorithm_chunk_gam_output,
            current_task_output = mergeAlignmentGAMChunks.merged_gam_file
        }
    }
    if (PACK_MODE) {
        call runVGPackCaller {
            input: 
            in_sample_name=SAMPLE_NAME, 
            in_xg_file=XG_FILE, 
            in_gam_file=mergeAlignmentGAMChunks.merged_gam_file,
            in_vg_container=VG_CONTAINER,
            in_vgcall_cores=VGCALL_CORES,
            in_vgcall_disk=VGCALL_DISK,
            in_vgcall_mem=VGCALL_MEM
        }
        # Cleanup vg call input files after use
        if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpVGPackCallInputsGoogle {
                input:
                previous_task_outputs = [mergeAlignmentGAMChunks.merged_gam_file],
                current_task_output = runVGPackCaller.output_vcf
            }
        }
        if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpVGPackCallInputsUnix {
                input:
                previous_task_outputs = [mergeAlignmentGAMChunks.merged_gam_file],
                current_task_output = runVGPackCaller.output_vcf
            }
        }
    }
    if (!PACK_MODE) {
        # Set chunk context amount depending on structural variant or default variant calling mode
        Int default_or_sv_chunk_context = if SV_CALLER_MODE then 2500 else 50
        # Chunk GAM alignments by contigs list
        call chunkAlignmentsByPathNames {
            input:
            in_sample_name=SAMPLE_NAME,
            in_merged_sorted_gam=mergeAlignmentGAMChunks.merged_sorted_gam_file,
            in_merged_sorted_gam_gai=mergeAlignmentGAMChunks.merged_sorted_gam_gai_file,
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
            if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpVGCallInputsGoogle {
                    input:
                    previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                    current_task_output = runVCFClipper.output_clipped_vcf
                }
            }
            if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
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
    File final_vcf_output = select_first([bgzipVGCalledVCF.output_merged_vcf, runVGPackCaller.output_vcf])
    output {
        File output_vcf = final_vcf_output
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
        docker: "ubuntu"
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
        memory: "20 GB"
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

task mergeAlignmentGAMChunks {
    input {
        String in_sample_name
        String in_vg_container
        Array[File] in_alignment_gam_chunk_files
        Int in_merge_gam_cores
        Int in_merge_gam_disk
        Int in_merge_gam_mem
        Int in_merge_gam_time
        Boolean sort_gam = true
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

    VG_GAMSORT_COMMAND="touch ~{in_sample_name}_merged.sorted.gam ~{in_sample_name}_merged.sorted.gam.gai"
    if [ ~{sort_gam} == true ]; then
        VG_GAMSORT_COMMAND="vg gamsort ~{in_sample_name}_merged.gam -i ~{in_sample_name}_merged.sorted.gam.gai -t ~{in_merge_gam_cores} > ~{in_sample_name}_merged.sorted.gam"
    fi
    
    cat ~{sep=" " in_alignment_gam_chunk_files} > ~{in_sample_name}_merged.gam
    ${VG_GAMSORT_COMMAND}
    >>>
    output {
        File merged_sorted_gam_file = "~{in_sample_name}_merged.sorted.gam"
        File merged_sorted_gam_gai_file = "~{in_sample_name}_merged.sorted.gam.gai"
        File merged_gam_file = "~{in_sample_name}_merged.gam"
    }
    runtime {
        memory: in_merge_gam_mem + " GB"
        cpu: in_merge_gam_cores
        disks: "local-disk " + in_merge_gam_disk  + " SSD"
        time: in_merge_gam_time
        docker: in_vg_container
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
    >>>
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

task runVGPackCaller {
    input {
        String in_sample_name
        File in_xg_file
        File in_gam_file
        String in_vg_container
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    String graph_tag = basename(in_xg_file, ".xg")
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -eux -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        vg pack \
           -x ~{in_xg_file} \
           -g ~{in_gam_file} \
           -q \
           -t ~{in_vgcall_cores} \
           -o ~{graph_tag}.pack
        
        vg call \
           -k ~{graph_tag}.pack \
           -t ~{in_vgcall_cores} \
           -s ~{in_sample_name} \
           ~{in_xg_file} > ~{graph_tag}.vcf

        head -10000 ~{graph_tag}.vcf | grep "^#" >> ~{graph_tag}.sorted.vcf
        if [ "$(cat ~{graph_tag}.vcf | grep -v '^#')" ]; then
            cat ~{graph_tag}.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ~{graph_tag}.sorted.vcf
        fi
        bgzip ~{graph_tag}.sorted.vcf && \
            tabix -f -p vcf ~{graph_tag}.sorted.vcf.gz
    >>>
    output {
        File output_vcf = "~{graph_tag}.sorted.vcf.gz"
        File output_vcf_index = "~{graph_tag}.sorted.vcf.gz.tbi"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        maxRetries: 3
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
