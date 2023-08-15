version 1.0


task createDistanceIndex {
    input {
        File in_gbz_file
        Int nb_cores = 16
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
    }
    String output_prefix = sub(basename(in_gbz_file), "\\.gbz$", "")

    command {
        set -eux -o pipefail

        vg index -t ~{nb_cores} -j "~{output_prefix}.dist" ~{in_gbz_file}
    }

    output {
        File output_dist_index = "~{output_prefix}.dist"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"

    }
}

task createRIndex {
    input {
        File in_gbz_file
        Int nb_cores = 16
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
    }

    String out_prefix_name = sub( basename(in_gbz_file), "\\.gbz$", "")

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

        vg gbwt -p --num-threads ~{nb_cores} -r ~{out_prefix_name}.ri -Z ~{in_gbz_file}

    }

    output {
        File output_R_index = "~{out_prefix_name}.ri"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"

    }

}

task createHaplotypeIndex {
    input {
        File in_gbz_file
        File in_dist_index
        File in_R_index
        Int kmer_length
        Int window_length
        Int subchain_length
        Int nb_cores = 16
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G") + size(in_dist_index, "G") + size(in_R_index, "G")) + 20
    }

    String out_prefix_name = sub( basename(in_gbz_file), "\\.gbz$", "")

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


        vg haplotypes -v 2 --kmer-length ~{kmer_length} \
        --window-length ~{window_length} \
        --subchain-length ~{subchain_length} \
        -t ~{nb_cores} -d ~{in_dist_index} \
        -r ~{in_R_index} -H ~{out_prefix_name}.hapl ~{in_gbz_file}

    }

    output {
        File output_hap_index = "~{out_prefix_name}.hapl"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"

    }

}


task createMinimizerIndex {
    input {
        File in_gbz_file
        File in_dist_index
        String out_name
        Int nb_cores = 16
        Int in_extract_mem = 120
        Int in_extract_disk = 4 * round(size(in_gbz_file, "G") + size(in_dist_index, "G")) + 20
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

        vg minimizer -p -t ~{nb_cores} -o ~{out_name}.min -d ~{in_dist_index} ~{in_gbz_file}

    }

    output {
        File output_minimizer = "~{out_name}.min"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"

    }

}

