version 1.0

workflow vgMapCallSV {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Read mapping and SV genotyping using vg. This workflow uses the new and faster mapper (giraffe/gaffe). It takes a CRAM file and graphs containing the structural variants to genotype. The XG, minimizer and distance graph indexes as required, as well as the original VCF used to create the graph. Including the snarls index is optional but speeds up the computation. It outputs a VCF file. The CRAM files are converted to FASTQ in chunks. This is used to parallelize the CRAM conversion jobs and the mapping jobs. It also makes them short enough that they can be run on cheaper instances (e.g. pre-emptible). The NB_CRAM_CHUNKS parameter controls the number of chunks (~15 recommended for 20x Illumina WGS). The MAX_CRAM_CHUNKS parameter can be used to down-sample the reads by using only chunks up to MAX_CRAM_CHUNKS. For example, with NB_CRAM_CHUNKS=15 and MAX_CRAM_CHUNKS=5 only the first 5 chunks will be used, leading to downsampling 1/3 of the reads. To allow for preemptible instance, increase the PREEMPTIBLE attempts number parameter, e.g. PREEMPTIBLE=3." 
    }
    input {
        String SAMPLE_NAME                                  # The sample name
        File INPUT_CRAM_FILE                                # Input CRAM file
        File CRAM_REF                                       # Genome fasta file associated with the CRAM file
        File CRAM_REF_INDEX                                 # Index of the fasta file associated with the CRAM file
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.22.0"   # VG Container used in the pipeline
        File XG_FILE                                        # Path to .xg index file
        File? VCF_FILE                                      # Path to the VCF used to create the graph
        File GBWT_FILE                                      # Path to .gbwt index file
        File DIST_IDX_FILE                                  # Path to .dist index file
        File MIN_IDX_FILE                                   # Path to .min index file
        File? SNARL_FILE                                    # (OPTIONAL) Path to snarl index file
        Int NB_CRAM_CHUNKS = 15                             # Number of chunks to split the reads in
        Int MAX_CRAM_CHUNKS = 15                            # Number of chunks to actually analyze
        Int CRAM_CONVERT_CORES = 3                          # Resources for the different tasks
        Int CRAM_CONVERT_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAM_DISK = 400
        Int VGCALL_CORES = 10
        Int VGCALL_DISK = 200
        Int VGCALL_MEM = 100
        Int PREEMPTIBLE = 0                                # Number of attempts to use pre-emptible instances
        Boolean GOOGLE_CLEANUP_MODE = false                # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        Boolean CLEANUP_MODE = true                        # Set to 'false' to disable cleanup (useful to be able to use use call caching to restart workflows)
    }

    # Split the reads (CRAM file) in chunks (FASTQ) and map each chunk separately.
    scatter (chunk_id in range(MAX_CRAM_CHUNKS)){
        call splitCramFastq {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_cram_convert_cores=CRAM_CONVERT_CORES,
            in_cram_convert_disk=CRAM_CONVERT_DISK,
            in_chunk_id=chunk_id,
            in_nb_chunks=NB_CRAM_CHUNKS,
            in_preemptible=PREEMPTIBLE
        }
        call runVGMAP {
            input:
            in_left_read_pair_chunk_file=splitCramFastq.fastq1,
            in_right_read_pair_chunk_file=splitCramFastq.fastq2,
            in_vg_container=VG_CONTAINER,
            in_xg_file=XG_FILE,
            in_min_file=MIN_IDX_FILE,
            in_dist_file=DIST_IDX_FILE,
            in_gbwt_file=GBWT_FILE,
            in_sample_name=SAMPLE_NAME,
            in_chunk_id=chunk_id,
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
                previous_task_outputs = [splitCramFastq.fastq1, splitCramFastq.fastq2],
                current_task_output = vg_map_algorithm_chunk_gam_output,
                in_preemptible=PREEMPTIBLE
            }
        }
        if (!GOOGLE_CLEANUP_MODE && CLEANUP_MODE) {
            call cleanUpUnixFilesystem as cleanUpVGMapperInputsUnix {
                input:
                previous_task_outputs = [splitCramFastq.fastq1, splitCramFastq.fastq2],
                current_task_output = vg_map_algorithm_chunk_gam_output,
                in_preemptible=PREEMPTIBLE
            }
        }
    }
    
    # Merge chunked graph alignments (GAM files)
    call mergeAlignmentGAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_vg_container=VG_CONTAINER,
        in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
        in_merge_gam_disk=MERGE_GAM_DISK,
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

    # Genotype variants using the packed graph approach
    call runVGPackCaller {
        input: 
        in_sample_name=SAMPLE_NAME,
        in_vcf_file=VCF_FILE,
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

    # Save the VCF and potentially the packed graph and graph alignment files
    output {
        File output_vcf = runVGPackCaller.output_vcf
        File output_vcf_index = runVGPackCaller.output_vcf_index
        File? output_pack = runVGPackCaller.output_pack
        File? output_gam = mergeAlignmentGAMChunks.merged_gam_file
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

# Retrieve a chunk from a CRAM file and convert it to FASTQ
# CRAM + (i,N) -> FASTQ_i_of_N
task splitCramFastq {
    input {
        File in_cram_file
        File in_ref_file
        File in_ref_index_file
        Int in_cram_convert_cores
        Int in_cram_convert_disk
        Int in_chunk_id
        Int in_nb_chunks
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

        samtools collate $VIEW_CORES -k ~{in_chunk_id} -K ~{in_nb_chunks} --reference ~{in_ref_file} -Ouf ~{in_cram_file} . | samtools fastq -1 R1.fastq.gz -2 R2.fastq.gz -0 o.fq.gz -s s.fq.gz -c 1 -N -
    >>>
    output {
        File fastq1='R1.fastq.gz'
        File fastq2='R2.fastq.gz'
    }
    runtime {
        cpu: in_cram_convert_cores
        memory: "50 GB"
        disks: "local-disk " + in_cram_convert_disk + " SSD"
        docker: "jmonlong/samtools-jm:release-1.19jm0.2.1"
        preemptible: in_preemptible
    }
}

# Map reads to a variation graph.
# FASTQ + GRAPH INDEXES -> GAM
task runVGMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_min_file
        File in_dist_file
        File in_gbwt_file
        String in_vg_container
        String in_sample_name
        Int in_chunk_id
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

        vg gaffe \
          -x ~{in_xg_file} \
          -m ~{in_min_file} \
          -d ~{in_dist_file} \
          -p \
          -N ~{in_sample_name} \
          --gbwt-name ~{in_gbwt_file} \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.~{in_chunk_id}.gam
    >>>
    output {
        File chunk_gam_file = in_sample_name + '.' + in_chunk_id + '.gam'
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
        preemptible: in_preemptible
    }
}

# Merge GAMS
# GAM_1_of_N + ... + GAM_N_of_N -> GAM
task mergeAlignmentGAMChunks {
    input {
        String in_sample_name
        String in_vg_container
        Array[File] in_alignment_gam_chunk_files
        Int in_merge_gam_disk
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

    cat ~{sep=" " in_alignment_gam_chunk_files} > ~{in_sample_name}_merged.gam
    >>>
    output {
        File merged_gam_file = "~{in_sample_name}_merged.gam"
    }
    runtime {
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + in_merge_gam_disk  + " SSD"
        docker: in_vg_container
        preemptible: in_preemptible
    }
}

# Create packed graph and genotype VCF
# GAM + GRAPH INDEXES + VCF_original -> PACK + VCF_genotyped
task runVGPackCaller {
    input {
        String in_sample_name
        File? in_vcf_file
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
    Boolean vcf_options = defined(in_vcf_file)
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
           -Q 5 \
           -t ~{in_vgcall_cores} \
           -o ~{graph_tag}.~{in_sample_name}.pack

        SNARL_OPTION_STRING=""
        if [ ~{snarl_options} == true ]; then
          SNARL_OPTION_STRING="--snarls ~{in_snarl_file}"
        fi

        VCF_OPTION_STRING=""
        if [ ~{vcf_options} == true ]; then
            tabix -f ~{in_vcf_file}
          VCF_OPTION_STRING="-v ~{in_vcf_file}"
        fi

        vg call \
           -k ~{graph_tag}.~{in_sample_name}.pack \
           -t ~{in_vgcall_cores} \
           -s ~{in_sample_name} ${SNARL_OPTION_STRING} ${VCF_OPTION_STRING} \
           ~{in_xg_file} > ~{graph_tag}.unsorted.vcf

        head -10000 ~{graph_tag}.unsorted.vcf | grep "^#" >> ~{graph_tag}.~{in_sample_name}.vcf
        if [ "$(cat ~{graph_tag}.unsorted.vcf | grep -v '^#')" ]; then
            cat ~{graph_tag}.unsorted.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ~{graph_tag}.~{in_sample_name}.vcf
        fi
        bgzip ~{graph_tag}.~{in_sample_name}.vcf && \
            tabix -f -p vcf ~{graph_tag}.~{in_sample_name}.vcf.gz
    >>>
    output {
        File output_vcf = "~{graph_tag}.~{in_sample_name}.vcf.gz"
        File output_vcf_index = "~{graph_tag}.~{in_sample_name}.vcf.gz.tbi"
        File output_pack = "~{graph_tag}.~{in_sample_name}.pack"
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
