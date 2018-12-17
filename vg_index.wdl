task vgIndexTask {
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
        disks: "local-disk " + diskGB "  HDD"
        docker: "${vg_docker}"
        preemptible : 0
    }

    output {
        File xg = "${output_name}.xg"
        File gcsa = "${output_name}.gcsa"
        File lcp = "${output_name}.lcp"
    }
}
workflow vgIndex {
    File vgFile
    Int? threads
    String? index_options
    String? vg_docker
    
    String output_name

    Int vgDiskSZ = ceil(size(vgFile, "GB") + 1000;

    call vgMapTask{
        input:
            vgFile=vgFile,
            diskGB=vgDiskSZ,
            threads=threads,
            memory=memory,
            output_name=output_name,
            index_options=index_options,
            vg_docker=vg_docker
    }

}
