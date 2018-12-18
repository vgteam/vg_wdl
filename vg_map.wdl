task vgMapTask {
    File fastq1
    File fastq2
    Int diskGB
    Int? threads
    
    File gcsa_ind
    File xg_ind
    String? map_options
    
    String output_name

    command <<<
        set -ex -o pipefail
        tar xvf ${vg_index_tar} && rm ${vg_index_tar} && \
        vg map -f ${fastq1} -f ${fastq2} -g ${gcsa_ind} -x ${xg_ind} -t ${threads} ${map_options} > ${output_name}.gam
    >>>

    runtime {
        cpu: "${threads}}"
        memory: "160000 MB"
        disks: "local-disk " + diskGB "  HDD"
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

    File vg_index_tar

    Int fqDiskSZ = ceil(size(fastq1, "GB") + size(fastq2, "GB") + 20) * 2;

    call vgMapTask{
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            diskGB=fqDiskSZ,
            threads=threads,
            output_name=output_name,
            map_options=map_options,
            vg_index_tar=vg_index_tar
    }

}
