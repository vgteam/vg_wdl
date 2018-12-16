task shardFastq{
    File fastq1
    File fastq2
    Int diskGB
    
    command <<<
        fastq_splitter.sh ${fastq1} ${fastq2}}
    >>>

    runtime{
        cpu: "8"
        memory : "4 GB"
        disks : "local-disk " + diskGB + " HDD"
        docker: "erictdawson/fastq-splitter"
    }

    output{
        Array[File] firstSplits = glob("*")
        Array[File] secondSplits = glob("*")
    }
}

task vg_map {
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

task concat {
    Array[File] shard_gam
    String output_name

    command <<<
        cat ${sep=" " shard_gam} > "${output_name}.gam"
    >>>

    runtime {
        disks: "local-disk 400 HDD"
    }

    output {
        File gam = "${output_name}.gam"
    }
}

workflow vg_map_shards {
    File fastq1
    File fastq2
    String output_name

    File vg_index_tar

    Int fqDiskSZ = ceil(size(fastq1, "GB") + size(fastq2, "GB") + 20)

    call shardFastq{
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            diskGB=fqDiskSZ
    }

}
