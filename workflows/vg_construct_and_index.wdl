version 1.0

# Construct variation graph from reference genome FASTA and VCF, then generate
# xg and GCSA+lcp indices.
# Based on https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph
workflow vg_construct_and_index {
    input {
        String graph_name
        File ref_fasta_gz
        File vcf_gz
        Array[String]+ contigs = [
             "1",  "2",  "3",  "4",  "5",  "6",
             "7",  "8",  "9", "10", "11", "12",
            "13", "14", "15", "16", "17", "18",
            "19", "20", "21", "22",  "X",  "Y"
        ]
        String vg_docker = "quay.io/vgteam/vg:v1.13.0"
    }

    # make graph for each reference contig
    scatter (contig in contigs) {
        call construct_graph { input:
            ref_fasta_gz = ref_fasta_gz,
            vcf_gz = vcf_gz,
            contig = contig,
            vg_docker = vg_docker
        }
    }

    # combine them into a single graph with unique node IDs
    call combine_graphs { input:
        graph_name = graph_name,
        contigs_vg = construct_graph.contig_vg,
        vg_docker = vg_docker
    }

    # make xg index
    call xg_index { input:
        graph_name = graph_name,
        vg = combine_graphs.vg,
        vg_docker = vg_docker
    }

    # prune the graph of repetitive sequences in preparation for GCSA indexing
    scatter (contig_vg in combine_graphs.contigs_uid_vg) {
        call prune_graph { input:
            contig_vg = contig_vg,
            vg_docker = vg_docker
        }
    }

    # make GCSA index
    call gcsa_index { input:
        graph_name = graph_name,
        contigs_pruned_vg = prune_graph.contig_pruned_vg,
        vg_docker = vg_docker
    }

    output {
        File vg = combine_graphs.vg
        File xg = xg_index.xg
        File gcsa = gcsa_index.gcsa
        File gcsa_lcp = gcsa_index.lcp
    }
}

# construct the graph for one reference contig
task construct_graph {
    input {
        File ref_fasta_gz
        File vcf_gz
        String contig
        String vg_construct_options="--node-max 32 --handle-sv"
        String vg_docker
    }

    command {
        set -ex -o pipefail
        pigz -dc ${ref_fasta_gz} > ref.fa
        tabix "${vcf_gz}"
        vg construct -R "${contig}" -C -r ref.fa -v "${vcf_gz}" --region-is-chrom ${vg_construct_options} > "${contig}.vg"
    }

    output {
        File contig_vg = "${contig}.vg"
    }

    runtime {
        docker: vg_docker
    }
}

# combine graphs from several contigs into one, with unique node IDs
task combine_graphs {
    input {
        String graph_name
        Array[File]+ contigs_vg
        String vg_docker
    }

    command {
        set -ex -o pipefail
        mkdir vg/
        cp ${sep=" " contigs_vg} vg/
        vg ids -j vg/*.vg
        mkdir concat
        cat vg/*.vg > "concat/${graph_name}.vg"
    }

    output {
        File vg = "concat/${graph_name}.vg"
        Array[File]+ contigs_uid_vg = glob("vg/*.vg")
    }

    runtime {
        docker: vg_docker
    }
}

task xg_index {
    input {
        String graph_name
        File vg
        String xg_options = ""
        String vg_docker
    }

    command {
        set -ex -o pipefail
        vg index -x "${graph_name}.xg" ~{xg_options} "~{vg}"
    }

    output {
        File xg ="${graph_name}.xg"
    }

    runtime {
        docker: vg_docker
    }
}

task prune_graph {
    input {
        File contig_vg
        String prune_options = ""
        String vg_docker
    }

    command {
        set -ex -o pipefail
        vg prune -r "${contig_vg}" ~{prune_options} > "$(basename '${contig_vg}' .vg).pruned.vg"
    }

    output {
        File contig_pruned_vg = glob("*.pruned.vg")[0]
    }

    runtime {
        docker: vg_docker
    }
}

task gcsa_index {
    input {
        String graph_name
        Array[File]+ contigs_pruned_vg
        String gcsa_options = ""
        String vg_docker
    }

    command {
        set -ex -o pipefail
        vg index -g "${graph_name}.gcsa" ${gcsa_options} ${sep=" " contigs_pruned_vg}
    }

    output {
        File gcsa = "${graph_name}.gcsa"
        File lcp = "${graph_name}.gcsa.lcp"
    }

    runtime {
        docker: vg_docker
    }
}
