task vgCallTask{
    File graph
    File augmented_graph
    File support_file
    File trans_file
    File ref_path

    String outbase

    String? vg_docker
    Int? threads
    Int? memory
    Int? diskGB

    command <<<
        vg call ${augmented_graph} -b ${graph} -s ${support_file} -z ${trans_file} -r ${ref_path} > ${outbase}.vcf
    >>>

    runtime {
        docker : "${vg_docker}"
        memory : ${memory} + " GB"
        disks : "local-disk " + ${diskGB} + " HDD"
        preemptible : 2
    }

    output {
        File callsVCF = "${outbase}.vcf"
    }
}

worfklow vgCallWorkflow{
    File graph
    File augmented_graph
    File support_file
    File trans_file
    String ref_path

    String? vg_docker
    Int? threads
    Int? memory
    Int? diskBaseGB

    Int diskGB = ceil(size(graph, "GB") + size(augmented_graph, "GB") + size(support_file) + size(trans_file)) + 20;

    call vgCallTask{
        graph=graph,
        augmented_graph=augmented_graph,
        support_file=support_file,
        trans_file=trans_file,
        ref_path=ref_path,
        outbase=outbase,
        vg_docker=vg_docker,
        threads=threads,
        memory=memory,
        diskGB=diskGB

    }

}
