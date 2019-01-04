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
    }

    output{
        Array[File] firstSplits = glob("*")
        Array[File] secondSplits = glob("*")
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

