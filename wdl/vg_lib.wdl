

task vg_map_shard {
    File fastq1
    File fastq2
    Int shard_number
    
    File vg_index_tar
    String? map_options
    
    String vg_docker_image

    String output_name

    command <<<
        set -ex -o pipefail
        tar xvf "${vg_index_tar}" && rm "${vg_index_tar}"
        vg map -f "${fastq1}" -f "${fastq2}" -d vg/index -t $(nproc) ${map_options} > "${output_name}.${shard_number}.gam"
    >>>

    runtime {
        cpu: "32"
        memory: "160000 MB"
        disks: "local-disk 200 HDD"
        docker: "${vg_docker_image}"
        preemptible : 3
    }

    output {
        File gam = "${output_name}.${shard_number}.gam"
    }
}

task vg_map {
    File fastq1
    File fastq2
    Int diskGB
    Int? threads

    File gcsa_index
    File xg_index
    File lcp_index
    String? map_options

    String output_name

    command <<<
        set -ex -o pipefail
        vg map -f ${fastq1} -f ${fastq2} -g ${gcsa_index} -x ${xg_index} -t ${threads} ${map_options} > ${output_name}.gam
    >>>

    runtime {
        cpu: "${threads}}"
        memory: "160000 MB"
        disks: "local-disk " + diskGB + "  HDD"
        docker: "variationgraphs/vg:latest"
        preemptible : 3
    }

    output {
        File gam = "${output_name}.gam"
    }
}

task vg_index {
    File vgFile

    Int diskGB
    Int? memory
    Int? threads
    String? index_options
    String? vg_docker

    String output_name

    command <<<
        set -ex -o pipefail
        vg index -t ${threads} -x ${output_name}.xg -g ${output_name}.gcsa ${index_options} ${vgFile}
    >>>

    runtime {
        cpu: "${threads}"
        memory: "${memory} GB"
        disks: "local-disk " + diskGB + "  HDD"
        docker: "${vg_docker}"
        preemptible : 0
    }

    output {
        File xg = "${output_name}.xg"
        File gcsa = "${output_name}.gcsa"
        File lcp = "${output_name}.lcp"
    }
}

task vg_sift{

}

task vg_call{

}

task vg_recall{

}
