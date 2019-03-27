version 1.0

# Construct variation graph from reference genome FASTA and VCF, then generate
# xg and GCSA+lcp indices.
# Based on https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph
workflow vg_construct_and_index {
    input {
        String graph_name
        File ref_fasta_gz
        Array[String]+ contigs = [
             "1",  "2",  "3",  "4",  "5",  "6",
             "7",  "8",  "9", "10", "11", "12",
            "13", "14", "15", "16", "17", "18",
            "19", "20", "21", "22",  "X",  "Y"
        ]
        Array[File]+ contigs_vcf_gz
        Boolean use_haplotypes = false
        String vg_docker = "quay.io/vgteam/vg:v1.14.0"
    }

    # make graph for each reference contig
    scatter (i in range(length(contigs))) {
        call construct_graph { input:
            ref_fasta_gz = ref_fasta_gz,
            vcf_gz = contigs_vcf_gz[i],
            contig = contigs[i],
            use_haplotypes = use_haplotypes,
            vg_docker = vg_docker
        }

        if (use_haplotypes) {
            call gbwt_index { input:
                vg = construct_graph.contig_vg,
                vcf_gz = contigs_vcf_gz[i],
                vg_docker = vg_docker
            }
        }
    }

    if (use_haplotypes) {
        Array[File]+ gbwt_threads = select_all(gbwt_index.threads)
        call gbwt_merge { input:
            gbwts = select_all(gbwt_index.gbwt),
            graph_name = graph_name,
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
        threads = gbwt_threads,
        vg_docker = vg_docker
    }

    # prune the graph of repetitive sequences in preparation for GCSA indexing
    if (!use_haplotypes) {
        scatter (contig_vg in combine_graphs.contigs_uid_vg) {
            call prune_graph { input:
                contig_vg = contig_vg,
                vg_docker = vg_docker
            }
        }
    }
    if (use_haplotypes) {
        call prune_graph_with_haplotypes { input:
            contigs_vg = combine_graphs.contigs_uid_vg,
            contigs_gbwt = select_all(gbwt_index.gbwt),
            empty_id_map = combine_graphs.empty_id_map,
            vg_docker = vg_docker
        }
    }

    # make GCSA index
    call gcsa_index { input:
        graph_name = graph_name,
        contigs_pruned_vg = select_first([prune_graph.contig_pruned_vg, prune_graph_with_haplotypes.contigs_pruned_vg]),
        empty_id_map = combine_graphs.empty_id_map,
        vg_docker = vg_docker
    }

    output {
        File vg = combine_graphs.vg
        File xg = xg_index.xg
        File? gbwt = gbwt_merge.gbwt
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
        Boolean use_haplotypes
        String vg_construct_options="--node-max 32 --handle-sv"
        String vg_docker
    }

    command {
        set -ex -o pipefail
        pigz -dc ${ref_fasta_gz} > ref.fa
        tabix "${vcf_gz}"

        vg construct -R "${contig}" -C -r ref.fa -v "${vcf_gz}" --region-is-chrom ${vg_construct_options} ${if use_haplotypes then "-a" else ""} > "${contig}.vg"
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
        # we approach this in a particular way to ensure the output array contigs_uid_vg has the
        # same order as the input array contigs_vg (so we can't rely on glob patterns)
        mkdir vg/
        while read -r contig_vg; do
            nm=$(basename "$contig_vg")
            cp "$contig_vg" "vg/$nm"
            echo "vg/$nm" >> contigs_uid_vg
        done < "~{write_lines(contigs_vg)}"
        xargs -n 999999 vg ids -j -m empty.id_map < contigs_uid_vg
        mkdir concat
        xargs -n 999999 cat < contigs_uid_vg > "concat/${graph_name}.vg"
    }

    output {
        File vg = "concat/${graph_name}.vg"
        File empty_id_map = "empty.id_map"
        Array[File]+ contigs_uid_vg = read_lines("contigs_uid_vg")
    }

    runtime {
        docker: vg_docker
    }
}

task gbwt_index {
    input {
        File vg
        File vcf_gz
        String vg_docker
    }

    command {
        set -ex -o pipefail
        nm=$(basename "~{vg}" .vg)
        vg index -G "$nm.gbwt" -F "$nm.threads" -v "~{vcf_gz}" "~{vg}"
    }

    output {
        File gbwt = glob("*.gbwt")[0]
        File threads = glob("*.threads")[0]
    }

    runtime {
        docker: vg_docker
    }
}

task gbwt_merge {
    input {
        Array[File]+ gbwts
        String graph_name
        String vg_docker
    }

    command {
        set -ex -o pipefail
        vg gbwt -m -f -o "~{graph_name}.gbwt" ~{sep=" " gbwts}
    }

    output {
        File gbwt = "~{graph_name}.gbwt"
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
        Array[File]+? threads
        String vg_docker
    }
    Array[File] threads_args = prefix("-F ", select_first([threads, []]))

    command {
        set -ex -o pipefail
        vg index -x "${graph_name}.xg" ~{xg_options} ~{sep=" " threads_args} "~{vg}"
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
        nm=$(basename "${contig_vg}" .vg)
        vg prune -r "${contig_vg}" ~{prune_options} > "$nm.pruned.vg"
    }

    output {
        File contig_pruned_vg = glob("*.pruned.vg")[0]
    }

    runtime {
        docker: vg_docker
    }
}

task prune_graph_with_haplotypes {
    input {
        Array[File]+ contigs_vg
        Array[File]+ contigs_gbwt
        File empty_id_map
        String prune_options = ""
        String vg_docker
    }

    command <<<
        set -ex -o pipefail
        cp "~{empty_id_map}" mapping
        paste -d ";" "~{write_lines(contigs_vg)}" "~{write_lines(contigs_gbwt)}" > inputs
        while IFS=';' read -ra p; do
            contig_vg="${p[0]}"
            contig_gbwt="${p[1]}"
            nm=$(basename "${contig_vg}" .vg)
            contig_pruned_vg="${nm}.pruned.vg"
            vg prune -u -g "$contig_gbwt" -a -m mapping "$contig_vg" ~{prune_options} > "$contig_pruned_vg"
            echo "$contig_pruned_vg" >> contigs_pruned_vg
        done <<< "inputs"
    >>>

    output {
        Array[File]+ contigs_pruned_vg=read_lines("contigs_pruned_vg")
    }

    runtime {
        docker: vg_docker
    }
}

task gcsa_index {
    input {
        String graph_name
        Array[File]+ contigs_pruned_vg
        File empty_id_map
        String gcsa_options = ""
        String vg_docker
    }

    command {
        set -ex -o pipefail
        vg index -g "${graph_name}.gcsa" -f "${empty_id_map}" ${gcsa_options} ${sep=" " contigs_pruned_vg}
    }

    output {
        File gcsa = "${graph_name}.gcsa"
        File lcp = "${graph_name}.gcsa.lcp"
    }

    runtime {
        docker: vg_docker
    }
}
