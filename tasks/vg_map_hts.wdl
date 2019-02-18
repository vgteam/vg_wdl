version 1.0

# vg map the reads in a SAM/BAM/CRAM file
task vg_map_hts {
    input {
        File sam_bam_cram
        File xg
        File gcsa
        File gcsa_lcp
        String vg_map_options = ""
        String vg_docker
    }

    command <<<
        set -ex -o pipefail
        ofn=$(basename "~{sam_bam_cram}" .sam)
        ofn=$(basename "$ofn" .bam)
        ofn=$(basename "$ofn" .cram)
        vg map -t "$(nproc)" ~{vg_map_options} -x "~{xg}" -g "~{gcsa}" --hts-input "~{sam_bam_cram}" > "${ofn}.gam"
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        File gam = glob("*.gam")[0]
    }
}
