task vgMapTask {
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
workflow vg_map {
    File fastq1
    File fastq2
    Int? threads
    String? map_options
    
    String output_name

    File xg_index
    File gcsa_index
    File lcp_index

    Int fqDiskSZ = ceil(size(fastq1, "GB") + size(fastq2, "GB") + 20) * 2

    call vgMapTask{
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            diskGB=fqDiskSZ,
            threads=threads,
            output_name=output_name,
            map_options=map_options,
            xg_index=xg_index,
            gcsa_index=gcsa_index,
            lcp_index=lcp_index
    }

}
