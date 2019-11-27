version 1.0

workflow vgMapCallSV {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Read mapping and SV genotyping using vg. It takes as inputs reads in FASTQ and graphs containing the structural variants to genotype. The graphs files required include the XG and GCSA indexes. Including the snarls index is optional but speeds up the computation. It outputs a VCF file." 
    }
    input {
        String SAMPLE_NAME                      # The sample name
        File? INPUT_READ_FILE_1                 # Input sample 1st read pair fastq.gz
        File? INPUT_READ_FILE_2                 # Input sample 2nd read pair fastq.gz
        File? INPUT_CRAM_FILE                   # Input CRAM file
        File? CRAM_REF                          # Genome fasta file associated with the CRAM file
        File? CRAM_REF_INDEX                    # Index of the fasta file associated with the CRAM file
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.19.0"   # VG Container used in the pipeline
        File XG_FILE                            # Path to .xg index file
        File GCSA_FILE                          # Path to .gcsa index file
        File GCSA_LCP_FILE                      # Path to .gcsa.lcp index file
        File? SNARL_FILE                        # (OPTIONAL) Path to snarl index file
        File? GBWT_FILE                         # (OPTIONAL) Path to .gbwt index file
        File? PATH_LIST_FILE                    # (OPTIONAL) Text file where each line is a path name in the XG index
        Int READS_PER_CHUNK = 20000000          # Number of reads contained in each mapping chunk (20000000 for wgs).
        Int CRAM_CONVERT_CORES = 16
        Int CRAM_CONVERT_DISK = 200
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAM_CORES = 64
        Int MERGE_GAM_DISK = 400
        Int MERGE_GAM_MEM = 100
        Int MERGE_GAM_TIME = 1200
        Int VGCALL_CORES = 10
        Int VGCALL_DISK = 200
        Int VGCALL_MEM = 100
        Int PREEMPTIBLE = 0
        Boolean GOOGLE_CLEANUP_MODE = false     # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        Boolean CLEANUP_MODE = true             # Set to 'false' to disable cleanup (useful to be able to use use call caching to restart workflows)
    }

    if(!defined(INPUT_READ_FILE_1) && defined(INPUT_CRAM_FILE) && defined(CRAM_REF_INDEX) && defined(CRAM_REF)){
        call convertCram {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_cram_convert_cores=CRAM_CONVERT_CORES,
            in_cram_convert_disk=CRAM_CONVERT_DISK,
            in_preemptible=0
        }
    }
    File in_read_file1 = select_first([INPUT_READ_FILE_1, convertCram.fastq1])
    File in_read_file2 = select_first([INPUT_READ_FILE_2, convertCram.fastq2])
    # Split input reads into chunks for parallelized mapping
    call splitReads as firstReadPair {
        input:
        in_read_file=in_read_file1,
        in_pair_id="1",
        in_vg_container=VG_CONTAINER,
        in_reads_per_chunk=READS_PER_CHUNK,
        in_split_read_cores=SPLIT_READ_CORES,
        in_split_read_disk=SPLIT_READ_DISK,
        in_preemptible=PREEMPTIBLE
    }
    call splitReads as secondReadPair {
        input:
        in_read_file=in_read_file2,
        in_pair_id="2",
        in_vg_container=VG_CONTAINER,
        in_reads_per_chunk=READS_PER_CHUNK,
        in_split_read_cores=SPLIT_READ_CORES,
        in_split_read_disk=SPLIT_READ_DISK,
        in_preemptible=PREEMPTIBLE
    }
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    if (!defined(PATH_LIST_FILE)) {
        call extractPathNames {
            input:
            in_xg_file=XG_FILE,
            in_vg_container=VG_CONTAINER,
            in_preemptible=PREEMPTIBLE
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
            in_map_mem=MAP_MEM,
            in_preemptible=PREEMPTIBLE
        }
        
        File vg_map_algorithm_chunk_gam_output = runVGMAP.chunk_gam_file
        # Cleanup input reads after use
        if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpGoogleFilestore as cleanUpVGMapperInputsGoogle {
                input:
                previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                current_task_output = vg_map_algorithm_chunk_gam_output,
                in_preemptible=PREEMPTIBLE
            }
        }
        if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpVGMapperInputsUnix {
                input:
                previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                current_task_output = vg_map_algorithm_chunk_gam_output,
                in_preemptible=PREEMPTIBLE
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
        in_preemptible=PREEMPTIBLE
    }
    # Cleanup gam chunk files after use
    if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpGoogleFilestore as cleanUpGAMChunksGoogle {
            input:
            previous_task_outputs = vg_map_algorithm_chunk_gam_output,
            current_task_output = mergeAlignmentGAMChunks.merged_gam_file,
            in_preemptible=PREEMPTIBLE
        }
    }
    if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpUnixFilesystem as cleanUpGAMChunksUnix {
            input:
            previous_task_outputs = vg_map_algorithm_chunk_gam_output,
            current_task_output = mergeAlignmentGAMChunks.merged_gam_file,
            in_preemptible=PREEMPTIBLE
        }
    }
    call runVGPackCaller {
        input: 
        in_sample_name=SAMPLE_NAME, 
        in_xg_file=XG_FILE, 
        in_gam_file=mergeAlignmentGAMChunks.merged_gam_file,
        in_snarl_file=SNARL_FILE,
        in_vg_container=VG_CONTAINER,
        in_vgcall_cores=VGCALL_CORES,
        in_vgcall_disk=VGCALL_DISK,
        in_vgcall_mem=VGCALL_MEM,
        in_preemptible=PREEMPTIBLE
    }
    # Cleanup vg call input files after use
    if (GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpGoogleFilestore as cleanUpVGPackCallInputsGoogle {
            input:
            previous_task_outputs = [mergeAlignmentGAMChunks.merged_gam_file],
            current_task_output = runVGPackCaller.output_vcf,
            in_preemptible=PREEMPTIBLE
        }
    }
    if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
        call cleanUpUnixFilesystem as cleanUpVGPackCallInputsUnix {
            input:
            previous_task_outputs = [mergeAlignmentGAMChunks.merged_gam_file],
            current_task_output = runVGPackCaller.output_vcf,
            in_preemptible=PREEMPTIBLE
        }
    }
    
    # Extract either the linear-based or graph-based VCF
    File final_vcf_output = runVGPackCaller.output_vcf
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
        Int in_preemptible
    }
    command <<<
        set -eux -o pipefail
        cat ~{write_lines(previous_task_outputs)} | sed 's/.*\(\/cromwell-executions\)/\1/g' | xargs -I {} ls -li {} | cut -f 1 -d ' ' | xargs -I {} find ../../../ -xdev -inum {} | xargs -I {} rm -v {}
    >>>
    runtime {
        docker: "ubuntu"
        continueOnReturnCode: true
        preemptible: in_preemptible
    }
}

task cleanUpGoogleFilestore {
    input {
        Array[String] previous_task_outputs
        String current_task_output
        Int in_preemptible
    }
    command {
        set -eux -o pipefail
        gsutil rm -I < ${write_lines(previous_task_outputs)}
    }
    runtime {
        docker: "google/cloud-sdk"
        continueOnReturnCode: true
        preemptible: in_preemptible
    }
}

task convertCram {
    input {
        File? in_cram_file
        File? in_ref_file
        File? in_ref_index_file
        Int in_cram_convert_cores
        Int in_cram_convert_disk
        Int in_preemptible
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

        VIEW_CORES=""
        if [ ~{in_cram_convert_cores} -gt 3 ]; then
            VIEW_CORES="-@ $(( ~{in_cram_convert_cores} - 3 ))"
        fi

        samtools view $VIEW_CORES -h -T ~{in_ref_file} ~{in_cram_file} | samtools collate -Ouf - . | samtools fastq -1 R1.fastq.gz -2 R2.fastq.gz -s singletons.fastq.gz -0 weird.fastq.gz -N -c 1 -
    >>>
    output {
        File fastq1='R1.fastq.gz'
        File fastq2='R2.fastq.gz'
    }
    runtime {
        cpu: in_cram_convert_cores
        memory: "50 GB"
        disks: "local-disk " + in_cram_convert_disk + " SSD"
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        preemptible: in_preemptible
    }
}

task splitReads {
    input {
        File? in_read_file
        String in_pair_id
        String in_vg_container
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk
        Int in_preemptible
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
        preemptible: in_preemptible
    }
}

task extractPathNames {
    input {
        File in_xg_file
        String in_vg_container
        Int in_preemptible
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
        preemptible: in_preemptible
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
        Int in_preemptible
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
        preemptible: in_preemptible
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
        Boolean sort_gam = false
        Int in_preemptible
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
        preemptible: in_preemptible
    }
}

task runVGPackCaller {
    input {
        String in_sample_name
        File in_xg_file
        File in_gam_file
        File? in_snarl_file
        String in_vg_container
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
        Int in_preemptible
    }

    Boolean snarl_options = defined(in_snarl_file)
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

        SNARL_OPTION_STRING=""
        if [ ~{snarl_options} == true ]; then
          SNARL_OPTION_STRING="--snarls ~{in_snarl_file}"
        fi
        
        vg call \
           -k ~{graph_tag}.pack \
           -t ~{in_vgcall_cores} \
           -s ~{in_sample_name} ${SNARL_OPTION_STRING} \
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
        preemptible: in_preemptible
    }
}
