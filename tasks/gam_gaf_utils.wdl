version 1.0

task mergeGAMandSort {
    input {
        Array[File] in_gam_files
        String in_sample_name = "sample"
        Int in_cores = 16
        Int in_mem = 120
    }
    Int disk_size = round(4 * size(in_gam_files, 'G')) + 20
    command <<<
        set -eux -o pipefail

        cat ~{sep=" " in_gam_files} | vg gamsort -p -i ~{in_sample_name}.sorted.gam.gai \
                                               -t ~{in_cores} - > ~{in_sample_name}.sorted.gam
    >>>
    output {
        File gam = "~{in_sample_name}.sorted.gam"
        File gam_index = "~{in_sample_name}.sorted.gam.gai"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.43.0"
    }
}

task mergeGAFandSort {
    input {
        Array[File] in_gaf_files
        File in_xg_file
        String in_sample_name = "sample"
        Int in_cores = 16
        Int in_mem = 120
    }
    Int disk_size = 4 * round(size(in_xg_file, 'G') + size(in_gaf_files, 'G')) + 20
    Int half_cores = if in_cores > 1 then floor(in_cores/2) else 1
    command <<<
        set -eux -o pipefail

        zcat ~{sep=" " in_gaf_files} | vg convert -t ~{half_cores} -F - ~{in_xg_file} | \
            vg gamsort -p -i ~{in_sample_name}.sorted.gam.gai \
               -t ~{half_cores} - > ~{in_sample_name}.sorted.gam
    >>>
    output {
        File gam = "~{in_sample_name}.sorted.gam"
        File gam_index = "~{in_sample_name}.sorted.gam.gai"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.43.0"
    }
}

task splitGAM {
    input {
        File in_gam_file
	    Int in_read_per_chunk
        Int in_mem = 30
        Int in_cores = 6
    }
    Int disk_size = 3 * round(size(in_gam_file, 'G')) + 20

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

        vg chunk -t ~{in_cores} -a ~{in_gam_file} -b chunk -m ~{in_read_per_chunk}
    >>>
    output {
        Array[File] gam_chunk_files = glob("chunk*.gam")
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.43.0"
    }
}

task splitGAF {
    input {
        File in_gaf_file
	    Int in_read_per_chunk
        Int in_mem = 8
        Int in_cores = 2
    }
    Int disk_size = 3 * round(size(in_gaf_file, 'G')) + 20

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

        zcat ~{in_gaf_file} | split -l ~{in_read_per_chunk} --filter='gzip > $FILE.gaf.gz' - chunked_gaf
    >>>
    output {
        Array[File] gaf_chunk_files = glob("chunk*.gaf.gz")
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.43.0"
    }
}

task surjectGAMtoSortedBAM {
    input {
        File in_gam_file
        File in_xg_file
        File in_path_list_file
        String in_sample_name
        Boolean in_is_paired_end = true
        Int in_max_fragment_length = 3000
        Boolean make_bam_index = false
        Boolean input_is_gaf = false
        Int in_map_cores = 16
        String in_map_mem = 120
    }
    String out_prefix = basename(in_gam_file, ".gam")
    Int half_cores = in_map_cores / 2
    Int disk_size = 3 * round(size(in_xg_file, 'G') + size(in_gam_file, 'G')) + 20
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

        PAIR_ARGS=""
        if [ ~{in_is_paired_end} == true ]
        then
            PAIR_ARGS="--interleaved --max-frag-len ~{in_max_fragment_length}"
        fi
        
        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_xg_file} \
          -t ~{half_cores} \
          --bam-output ~{true="--gaf-input" false="" input_is_gaf} \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx $PAIR_ARGS \
          ~{in_gam_file} | samtools sort --threads ~{half_cores} \
                                    -O BAM > ~{out_prefix}.bam

        if [ ~{make_bam_index} == true ]; then
            samtools index ~{out_prefix}.bam ~{out_prefix}.bam.bai
        fi
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
        File? output_bam_index_file = "~{out_prefix}.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.43.0"
    }
}

