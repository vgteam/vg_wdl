version 1.0

workflow vgMapCallSV {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Read mapping and SV genotyping using vg. This workflow uses the new and faster mapper (giraffe). It takes a CRAM file and indexes for a graph containing the structural variants to genotype. The XG, minimizer and distance graph indexes as required. Including the snarls index is optional but speeds up the computation. \n\nIt outputs a VCF file. By default, the VCF contains SV calls in the context of the snarls in the graph (i.e. variant site/bubbles). This helps merging the calls from multiple samples that used the same input graph. However it contains additional flanking sequences in the REF/ALT fields so should be normalized (e.g. using 'bcftools norm -f REF.fa VCF.vcf') before comparing with other SV catalogs or annotating the variants. If you are sure you don't want this 'snarl calling' behavior, set SNARL_CALLING to false. Finally, if the graph was built using a VCF and you want to genotype SVs using that same VCF, pass it to VCF_FILE and switch off the SNARL_CALLING which is not useful in that case.\n\nThe CRAM files are converted to FASTQ in chunks. This is used to parallelize the CRAM conversion jobs and the mapping jobs. It also makes them short enough that they can be run on cheaper instances (e.g. pre-emptible). The NB_CRAM_CHUNKS parameter controls the number of chunks (~15 recommended for 20x Illumina WGS). The MAX_CRAM_CHUNKS parameter can be used to down-sample the reads by using only chunks up to MAX_CRAM_CHUNKS. For example, with NB_CRAM_CHUNKS=15 and MAX_CRAM_CHUNKS=5 only the first 5 chunks will be used, leading to downsampling 1/3 of the reads. To allow for preemptible instance, increase the PREEMPTIBLE attempts number parameter, e.g. PREEMPTIBLE=3." 
    }
    input {
        String SAMPLE_NAME                                  # The sample name
        File INPUT_CRAM_FILE                                # Input CRAM file
        File CRAM_REF                                       # Genome fasta file associated with the CRAM file
        File CRAM_REF_INDEX                                 # Index of the fasta file associated with the CRAM file
        File XG_FILE                                        # Path to .xg index file
        File GBWT_FILE                                      # Path to .gbwt index file
        File DIST_IDX_FILE                                  # Path to .dist index file
        File MIN_IDX_FILE                                   # Path to .min index file
        File? VCF_FILE                                      # (OPTIONAL) Path to the VCF used to create the graph
        File? SNARL_FILE                                    # (OPTIONAL) Path to snarl index file
        File? REFPATH_FILE                                  # (OPTIONAL) Path to a file listing the reference paths
        Int NB_CRAM_CHUNKS = 8                              # Number of chunks to split the reads in
        Int MAX_CRAM_CHUNKS = 4                             # Number of chunks to actually analyze
        Boolean GIRAFFE_FAST_MODE = true                    # Set to 'false' to not use the fast mode for the read mapping with giraffe
        String GIRAFFE_RESCUE_MODE = "dozeu"                # Rescue mode for the giraffe mapper. Either "dozeu" (default) or "gssw" (slower but currently a bit more accurate)
        Boolean SNARL_CALLING = true                        # Return SV calls in the snarl context. Easier to merge results across samples, although requires 'bcftools norm -f REF.fa' to clean it up.
        Int CRAM_CONVERT_CORES = 4                          # Resources for the different tasks
        Int CRAM_CONVERT_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAF_DISK = 400
        Int VGCALL_CORES = 16
        Int VGCALL_DISK = 200
        Int VGCALL_MEM = 100
        Int PREEMPTIBLE = 0                                # Number of attempts to use pre-emptible instances
        Boolean OUT_GAF = false                            # Should the chunked GAFs (aligned reads) be merged and part of the final output. Default if false.
    }

    # Split the reads (CRAM file) in chunks (FASTQ)
    call splitCramFastq {
        input:
        in_cram_file=INPUT_CRAM_FILE,
        in_ref_file=CRAM_REF,
        in_ref_index_file=CRAM_REF_INDEX,
        in_cram_convert_cores=CRAM_CONVERT_CORES,
        in_cram_convert_disk=CRAM_CONVERT_DISK,
        in_nb_chunks=NB_CRAM_CHUNKS,
        in_max_chunks=MAX_CRAM_CHUNKS,
        in_preemptible=PREEMPTIBLE
    }

    # Map each chunk separately.
    scatter (chunk_id in range(length(splitCramFastq.output_read_chunks_1))) {
        call runVGMAP {
            input:
            in_left_read_pair_chunk_file=splitCramFastq.output_read_chunks_1[chunk_id],
            in_right_read_pair_chunk_file=splitCramFastq.output_read_chunks_2[chunk_id],
            in_xg_file=XG_FILE,
            in_min_file=MIN_IDX_FILE,
            in_dist_file=DIST_IDX_FILE,
            in_gbwt_file=GBWT_FILE,
            in_sample_name=SAMPLE_NAME,
            in_chunk_id=chunk_id,
            in_fast_mode=GIRAFFE_FAST_MODE,
            in_rescue_mode=GIRAFFE_RESCUE_MODE,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM,
            in_preemptible=PREEMPTIBLE
        }        
        File vg_map_algorithm_chunk_gaf_output = runVGMAP.chunk_gaf_file
    }
    
    # Genotype variants using the packed graph approach
    call runVGPackCaller {
        input: 
        in_sample_name=SAMPLE_NAME,
        in_vcf_file=VCF_FILE,
        in_xg_file=XG_FILE, 
        in_alignment_gaf_chunk_files=vg_map_algorithm_chunk_gaf_output,
        in_snarl_calling=SNARL_CALLING,
        in_snarl_file=SNARL_FILE,
        in_refpath_file=REFPATH_FILE,
        in_vgcall_cores=VGCALL_CORES,
        in_vgcall_disk=VGCALL_DISK,
        in_vgcall_mem=VGCALL_MEM,
        in_preemptible=PREEMPTIBLE
    }

    if(OUT_GAF){
        # Merge chunked graph alignments (GAF files)
        call mergeAlignmentGAFChunks {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_gaf_chunk_files=vg_map_algorithm_chunk_gaf_output,
            in_merge_gaf_disk=MERGE_GAF_DISK,
            in_preemptible=PREEMPTIBLE
        }
    }
    
    # Save the VCF and potentially the packed graph and graph alignment files
    output {
        File output_vcf = runVGPackCaller.output_vcf
        File output_vcf_index = runVGPackCaller.output_vcf_index
        File? output_gaf = mergeAlignmentGAFChunks.merged_gaf_file
    }
}

########################
### TASK DEFINITIONS ###
########################

# Retrieve chunks from a CRAM file and convert them to FASTQ
# CRAM + (k,N) -> FASTQ_1_of_N, FASTQ_2_of_N, ..., FASTQ_k_of_N
task splitCramFastq {
    input {
        File in_cram_file
        File in_ref_file
        File in_ref_index_file
        Int in_cram_convert_cores
        Int in_cram_convert_disk
        Int in_nb_chunks
        Int in_max_chunks
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

        seq 0 ~{in_nb_chunks} | head -n ~{in_max_chunks} | parallel -j ~{in_cram_convert_cores} "samtools collate -k {} -K ~{in_nb_chunks} --reference ~{in_ref_file} -Ouf ~{in_cram_file} {} | samtools fastq -1 reads.{}.R1.fastq.gz -2 reads.{}.R2.fastq.gz -0 reads.{}.o.fq.gz -s reads.{}.s.fq.gz -c 1 -N -"
    >>>
    output {
        Array[File] output_read_chunks_1 = glob("reads.*.R1.fastq.gz")
        Array[File] output_read_chunks_2 = glob("reads.*.R2.fastq.gz")
    }
    runtime {
        cpu: in_cram_convert_cores
        memory: "50 GB"
        disks: "local-disk " + in_cram_convert_disk + " SSD"
        docker: "jmonlong/samtools-jm:release-1.19jm0.2.2"
        preemptible: in_preemptible
    }
}

# Map reads to a variation graph.
# FASTQ + GRAPH INDEXES -> GAF
task runVGMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_min_file
        File in_dist_file
        File in_gbwt_file
        String in_sample_name
        Int in_chunk_id
        Boolean in_fast_mode
        String in_rescue_mode
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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

        MODE_ARG=""
        if [ ~{in_fast_mode} == true ]; then
          MODE_ARG="-b fast"
        fi
        
        vg giraffe \
          -x ~{in_xg_file} \
          -m ~{in_min_file} \
          -d ~{in_dist_file} \
          -p $MODE_ARG \
          --rescue-algorithm ~{in_rescue_mode} \
          -N ~{in_sample_name} \
          --gbwt-name ~{in_gbwt_file} \
          -C 500 \
          -o gaf \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -t ~{in_map_cores} | gzip > ~{in_sample_name}.~{in_chunk_id}.gaf.gz
    >>>
    output {
        File chunk_gaf_file = in_sample_name + '.' + in_chunk_id + '.gaf.gz'
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.28.0"
        preemptible: in_preemptible
    }
}

# Merge GAFS
# GAF_1_of_N + ... + GAF_N_of_N -> GAF
task mergeAlignmentGAFChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_gaf_chunk_files
        Int in_merge_gaf_disk
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

    cat ~{sep=" " in_alignment_gaf_chunk_files} > ~{in_sample_name}_merged.gaf.gz
    >>>
    output {
        File merged_gaf_file = "~{in_sample_name}_merged.gaf.gz"
    }
    runtime {
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + in_merge_gaf_disk  + " SSD"
        docker: "quay.io/vgteam/vg:v1.28.0"
        preemptible: in_preemptible
    }
}

# Create packed graph and genotype VCF
# GAF + GRAPH INDEXES + VCF_original -> PACK + VCF_genotyped
task runVGPackCaller {
    input {
        String in_sample_name
        File? in_vcf_file
        File in_xg_file
        Array[File] in_alignment_gaf_chunk_files
        Boolean in_snarl_calling
        File? in_snarl_file
        File? in_refpath_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
        Int in_preemptible
    }

    Boolean snarl_options = defined(in_snarl_file)
    Boolean vcf_options = defined(in_vcf_file)
    Boolean refpath_options = defined(in_refpath_file)
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

        cat ~{sep=" " in_alignment_gaf_chunk_files} | vg pack \
                                                         -x ~{in_xg_file} \
                                                         -a - \
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

        REFPATH_OPTION_STRING=""
        if [ ~{refpath_options} == true ]; then
            for RP in `cat ~{in_refpath_file}`
            do
                REFPATH_OPTION_STRING="-p $RP $REFPATH_OPTION_STRING"
            done
        fi

        SNARL_CALLING_OPTION_STRING=""
        if [ ~{in_snarl_calling} == true ]; then
            SNARL_CALLING_OPTION_STRING="-a"
        fi

        vg call \
           -k ~{graph_tag}.~{in_sample_name}.pack \
           -t ~{in_vgcall_cores} \
           -s ~{in_sample_name} ${SNARL_CALLING_OPTION_STRING} ${SNARL_OPTION_STRING} ${VCF_OPTION_STRING} ${REFPATH_OPTION_STRING} \
           ~{in_xg_file} > ~{graph_tag}.unsorted.vcf

        head -10000 ~{graph_tag}.unsorted.vcf | grep "^#" >> ~{graph_tag}.~{in_sample_name}.vcf
        if [ "$(grep -c -v '^#' ~{graph_tag}.unsorted.vcf)" -gt 0 ]; then
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
        docker: "quay.io/vgteam/vg:v1.28.0"
        preemptible: in_preemptible
    }
}
