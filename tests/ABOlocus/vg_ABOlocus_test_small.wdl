version 1.0

# This example uses vg_construct_and_index.wdl to build a graph including the
# ABO locus from GRCh38 (a 50Kbp region including the gene), with the 1000
# Genomes small variants (SNVs & short indels) and individual haplotypes.

import "../../workflows/vg_construct_and_index.wdl"

workflow vg_ABOlocus_test {
    input {
        File ABOlocus_fa_gz
        File ABOlocus_small_vcf_gz
        File reads_bam
        String vg_docker = "quay.io/vgteam/vg:v1.14.0"
    }

    # build & check the ABOlocus graph
    call vg_construct_and_index.vg_construct_and_index as cons { input:
        graph_name = "ABOlocus_small",
        ref_fasta_gz = ABOlocus_fa_gz,
        contigs = ["ABOlocus"],
        contigs_vcf_gz = [ABOlocus_small_vcf_gz],
        use_haplotypes = true,
        vg_docker = vg_docker
    }

    # extract FASTQs from BAM for use with vg mpmap
    call bam_to_paired_fastq { input:
        bam = reads_bam
    }

    # map reads to the ABOlocus graph & check the mappings
    call vg_mpmap { input:
        fastq_1 = bam_to_paired_fastq.fastq_1_gz,
        fastq_2 = bam_to_paired_fastq.fastq_2_gz,
        xg = cons.xg,
        gcsa = cons.gcsa,
        gcsa_lcp = cons.gcsa_lcp,
        gbwt = cons.gbwt,
        vg_docker = vg_docker
    }

    call check_gam { input:
        fastq_1_gz = bam_to_paired_fastq.fastq_1_gz,
        gam = vg_mpmap.gam,
        vg_docker = vg_docker
    }

    output {
        Int reads = check_gam.reads
        Int reads_aligned_identically = check_gam.reads_aligned_identically
    }
}

task bam_to_paired_fastq {
    input {
        File bam
    }

    command <<<
        set -eux -o pipefail
        apt-get update
        apt-get install -y openjdk-11-jre-headless wget samtools pigz
        wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar

        samtools sort -n -@ $(nproc) -O BAM -o namesorted.bam "~{bam}"
        nm=$(basename "~{bam}" .bam)
        java -jar picard.jar SamToFastq I=namesorted.bam RE_REVERSE=true INCLUDE_NON_PF_READS=true "FASTQ=${nm}_1.fastq" "SECOND_END_FASTQ=${nm}_2.fastq" "UNPAIRED_FASTQ=${nm}_unpaired.fastq" VALIDATION_STRINGENCY=LENIENT
        pigz *.fastq
    >>>

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File fastq_1_gz = glob("*_1.fastq.gz")[0]
        File fastq_2_gz = glob("*_2.fastq.gz")[0]
        File fastq_unpaired_gz = glob("*_unpaired.fastq.gz")[0]
    }
}

task vg_mpmap {
    input {
        File fastq_1
        File fastq_2
        File xg
        File gcsa
        File gcsa_lcp
        File? gbwt
        String vg_mpmap_options = ""
        String vg_docker
    }

    command <<<
        set -ex -o pipefail
        ofn=$(basename "~{fastq_1}" .gz)
        ofn=$(basename "$ofn" .fastq)
        ofn=$(basename "$ofn" .fq)
        ofn=$(basename "$ofn" _1)
        vg mpmap -t "$(nproc)" ~{vg_mpmap_options} --xg-name "~{xg}" --gcsa-name "~{gcsa}" ~{"--gbwt-name " + gbwt} --single-path-mode -f "~{fastq_1}" -f "~{fastq_2}" > "${ofn}.gam"
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        File gam = glob("*.gam")[0]
    }
}


task check_gam {
    # checks the gam has the expected # of mappings, and finds the proportion aligned at 100% identity
    input {
        File fastq_1_gz
        File gam
        String vg_docker
    }

    command <<<
        n_reads=$(zcat ~{fastq_1_gz} | wc -l)
        n_reads=$(expr "$n_reads" / 2)
        vg view -a "~{gam}" -j > mappings.json
        if [ "$n_reads" -ne "$(wc -l < mappings.json)" ]; then
            echo "wrong read count" >&2
            exit 1
        fi
        echo "$n_reads" > n_reads
        jq ".identity>0.99" mappings.json | grep true | wc -l > n_identical
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        Int reads = read_int("n_reads")
        Int reads_aligned_identically = read_int("n_identical")
    }
}
