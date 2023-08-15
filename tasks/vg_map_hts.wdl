version 1.0

# vg map the reads in a SAM/BAM/CRAM file

task vg_map_hts {
    input {
        File sam_bam_cram
        File xg
        File gcsa
        File gcsa_lcp
        File? gbwt
        String vg_map_options = ""
        String vg_docker
    }

    command <<<
        set -ex -o pipefail
        ofn=$(basename "~{sam_bam_cram}" .sam)
        ofn=$(basename "$ofn" .bam)
        ofn=$(basename "$ofn" .cram)
        vg map -t "$(nproc)" ~{vg_map_options} --xg-name "~{xg}" --gcsa-name "~{gcsa}" ~{"--gbwt-name " + gbwt} --hts-input "~{sam_bam_cram}" > "${ofn}.gam"
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        File gam = glob("*.gam")[0]
    }
}

task runVGGIRAFFE {
    input {
        File fastq_file_1
        File? fastq_file_2
        File in_gbz_file
        File in_dist_file
        File in_min_file
        String in_giraffe_options
        String in_sample_name
        Int nb_cores = 16
        String mem_gb = 120
        Int disk_size = 3 * round(size(fastq_file_1, 'G') + size(fastq_file_2, 'G') + size(in_gbz_file, 'G') + size(in_dist_file, 'G') + size(in_min_file, 'G')) + 50
    }

    String out_prefix = sub(sub(sub(basename(fastq_file_1), "\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
    Boolean paired_reads = defined(fastq_file_2)
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        PAIR_ARGS=""
        if [ ~{paired_reads} == true ]
        then
            PAIR_ARGS="-f ~{fastq_file_2}"
        fi

        vg giraffe \
        --progress \
        --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
        --sample "~{in_sample_name}" \
        ~{in_giraffe_options} \
        --output-format gaf \
        -f ~{fastq_file_1} ${PAIR_ARGS} \
        -Z ~{in_gbz_file} \
        -d ~{in_dist_file} \
        -m ~{in_min_file} \
        -t ~{nb_cores} | gzip > ~{out_prefix}.gaf.gz
    >>>
    output {
        File chunk_gaf_file = "~{out_prefix}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: mem_gb + " GB"
        cpu: nb_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"
    }
}

task extractSubsetPathNames {
    input {
        File in_gbz_file
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
        Int in_extract_mem = 120
    }

    command {
        set -eux -o pipefail

        vg gbwt -CL -Z ${in_gbz_file} | sort > path_list.txt

        grep -v _decoy path_list.txt | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM | grep -v chain_ > path_list.sub.txt
    }
    output {
        File output_path_list_file = "path_list.sub.txt"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"
    }
}

task extractReference {
    input {
        File in_gbz_file
        File in_path_list_file
        String in_prefix_to_strip = ""
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
    }

    command {
        set -eux -o pipefail

        # Subset to just the paths we care about (may be the whole file) so we
        # get a good dict with just those paths later
        vg paths \
        --extract-fasta \
        -p ${in_path_list_file} \
        --xg ${in_gbz_file} > ref.fa

        if [ ~{in_prefix_to_strip} != "" ]
        then
            mv ref.fa ref.prefix.fa
            sed -e "s/>~{in_prefix_to_strip}/>/g" ref.prefix.fa > ref.fa
        fi
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"
    }
}


task samplingHaplotypes {
    input {
        File in_gbz_file
        File in_hap_index
        File in_kmer_info
        String output_file_name
        String working_directory
        Int haplotype_number
        Float present_discount
        Float het_adjust
        Float absent_score
        Boolean include_reference
        Int nb_cores = 16
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G") + size(in_hap_index, "G") + size(in_kmer_info, "G")) + 20
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        INCLUDE_REF=""
        if [ ~{include_reference} == true ]
        then
            INCLUDE_REF="--include-reference"
        fi

        vg haplotypes -v 2 -t ~{nb_cores} \
        --num-haplotypes ~{haplotype_number} \
        --present-discount ~{present_discount} \
        --het-adjustment ~{het_adjust} \
        --absent-score ~{absent_score} \
        ${INCLUDE_REF} \
        -i ~{in_hap_index} \
        -k ~{in_kmer_info} \
        -g ~{working_directory}/~{output_file_name}.gbz ~{in_gbz_file}
    >>>

    output {
        File output_graph = working_directory+ "/" + output_file_name+".gbz"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.50.1"

    }

}


