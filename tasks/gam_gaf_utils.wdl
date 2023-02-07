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
        docker: "quay.io/jmonlong/vg:d213c29"
    }
}

task mergeGAFandSort {
    input {
        File in_gaf_file
        File in_gbz_file
        String in_sample_name = "sample"
        Int in_cores = 16
        Int in_mem = 120
        Int disk_size = 4 * round(size(in_gbz_file, 'G') + size(in_gaf_file, 'G')) + 20
    }

    Int half_cores = if in_cores > 1 then floor(in_cores/2) else 1
    command <<<
        set -eux -o pipefail

        vg convert -t ~{half_cores} -F ~{in_gaf_file} ~{in_gbz_file} | \
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
        docker: "quay.io/jmonlong/vg:d213c29"
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
        docker: "quay.io/jmonlong/vg:d213c29"
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
        docker: "quay.io/jmonlong/vg:d213c29"
    }
}

task mergeGAF {
    input {
        String in_sample_name
        Array[File] in_gaf_chunk_files
        Int in_disk = round(3*size(in_gaf_chunk_files, 'G')) + 20
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

        cat ~{sep=" " in_gaf_chunk_files} > ~{in_sample_name}.gaf.gz
    >>>
    output {
        File output_merged_gaf = "~{in_sample_name}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: "6GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/jmonlong/vg:d213c29"
    }
}

## task surjectGAMtoSortedBAM
task surjectGAFtoSortedBAM {
    input {
        File in_gaf_file
        File in_gbz_file
        File in_path_list_file
        String in_sample_name
        Boolean in_paired_reads = true
        Int in_max_fragment_length = 3000
        Boolean make_bam_index = false
        Boolean input_is_gam = false
        Int nb_cores = 16
        String mem_gb = 120
        Int disk_size = 5 * round(size(in_gbz_file, 'G') + size(in_gaf_file, 'G')) + 50
    }
    String out_prefix = sub(sub(sub(basename(in_gaf_file), "\\.gz$", ""), "\\.gaf$", ""), "\\.gam$", "")
    Int half_cores = nb_cores / 2
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
        if [ ~{in_paired_reads} == true ]
        then
            PAIR_ARGS="--interleaved --max-frag-len ~{in_max_fragment_length}"
        fi
        
        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_gbz_file} \
          -t ~{half_cores} \
          --bam-output ~{true="" false="--gaf-input" input_is_gam} \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx $PAIR_ARGS \
          ~{in_gaf_file} | samtools sort --threads ~{half_cores} \
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
        memory: mem_gb + " GB"
        cpu: nb_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/vg:d213c29"
    }
}

task surjectGAFtoBAM {
    input {
        File in_gaf_file
        File in_gbz_file
        File in_path_list_file
        String in_sample_name
        Boolean in_paired_reads = true
        Int in_max_fragment_length = 3000
        Boolean input_is_gam = false
        Int nb_cores = 16
        String mem_gb = 120
        Int disk_size = 5 * round(size(in_gbz_file, 'G') + size(in_gaf_file, 'G')) + 50
    }
    String out_prefix = sub(sub(sub(basename(in_gaf_file), "\\.gz$", ""), "\\.gaf$", ""), "\\.gam$", "")
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
        if [ ~{in_paired_reads} == true ]
        then
            PAIR_ARGS="--interleaved --max-frag-len ~{in_max_fragment_length}"
        fi
        
        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_gbz_file} \
          -t ~{nb_cores} \
          --bam-output ~{true="" false="--gaf-input" input_is_gam} \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx $PAIR_ARGS \
          ~{in_gaf_file} > ~{out_prefix}.bam
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: mem_gb + " GB"
        cpu: nb_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/vg:d213c29"
    }
}

