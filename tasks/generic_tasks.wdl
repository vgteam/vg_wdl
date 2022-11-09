version 1.0

task shard_fastq {
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
        preemptible : 3
    }

    output{
        Array[File] firstSplits = glob("*_1.fastq.part*")
        Array[File] secondSplits = glob("*_2.fastq.part*")
    }
}


task concat {
    Array[File] shard_gam
    String output_name
    Int diskGB

    command <<<
        cat ${sep=" " shard_gam} > "${output_name}.gam"
    >>>

    runtime {
        docker : "erictdawson/base"
        cpu : 1
        memory : "1 GB"
        disks: "local-disk " + diskGB + "  HDD"
        preemptible : 4
    }

    output {
        File gam = "${output_name}.gam"
    }
}

task falseTouchFile{
    File input
    Int diskGB

    command <<<
        touch ${input}
    >>>

    runtime{
        docker : "erictdawson/base"
        cpu : 1
        memory : "1 GB"
        disks : "local-disk " + diskGB + " HDD
        preemptible : 7
    }

    output{
        File output = "${input}
    }

}

workflow strawman_flow{

}
