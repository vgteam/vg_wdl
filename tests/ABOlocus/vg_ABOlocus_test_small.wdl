version 1.0

# This example uses vg_construct_and_index.wdl to build a graph including the
# ABO locus from GRCh38 (a 50Kbp region including the gene), with the 1000
# Genomes small variants (SNVs & short indels) and individual haplotypes.

import "../../workflows/vg_construct_and_index.wdl"
import "../../tasks/vg_map_hts.wdl"

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

    output {
    }
}
