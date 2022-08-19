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
        docker: "quay.io/vgteam/vg:v1.37.0"
    }
}

task mergeGAFandSort {
    input {
        Array[File] in_gaf_files
        File in_xg_file
        String in_sample_name = "sample"
        Int in_cores = 16
        Int in_mem = 120
    }
    Int disk_size = 4 * round(size(in_xg_file, 'G') + size(in_gaf_files, 'G')) + 20
    Int half_cores = if in_cores > 1 then floor(in_cores/2) else 1
    command <<<
        set -eux -o pipefail

        zcat ~{sep=" " in_gaf_files} | vg convert -t ~{half_cores} -F - ~{in_xg_file} | \
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
        docker: "quay.io/vgteam/vg:v1.37.0"
    }
}
