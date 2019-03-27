version 1.0

# Construct variation graph from reference genome FASTA and VCF, then generate
# xg and GCSA+lcp indices. Optionally, create GBWT index for phased haplotypes
# in the VCF.
# References:
#     https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph
#     https://github.com/vgteam/vg/wiki/Index-Construction
workflow vg_construct_and_index {
    input {
        # Overall name for the final graph (used in filenames)
        String graph_name

        # Reference genome FASTA
        File ref_fasta_gz

        # Desired reference genome contigs
        Array[String]+ contigs = [
             "1",  "2",  "3",  "4",  "5",  "6",
             "7",  "8",  "9", "10", "11", "12",
            "13", "14", "15", "16", "17", "18",
            "19", "20", "21", "22",  "X",  "Y"
        ]

        # VCF for each contig
        Array[File]+ contigs_vcf_gz

        # set true to GBWT index the VCF phased haplotypes
        Boolean use_haplotypes = false

        # vg docker image tag
        String vg_docker = "quay.io/vgteam/vg:v1.14.0"
    }

    # construct graph for each reference contig
    scatter (p in zip(contigs, contigs_vcf_gz)) {
        call construct_graph { input:
            ref_fasta_gz = ref_fasta_gz,
            contig = p.left,
            vcf_gz = p.right,
            use_haplotypes = use_haplotypes,
            vg_docker = vg_docker
        }
    }

    # combine them into a single graph with unique node IDs
    call combine_graphs { input:
        graph_name = graph_name,
        contigs_vg = construct_graph.contig_vg,
        vg_docker = vg_docker
    }

    # make GBWT index, if so configured
    if (use_haplotypes) {
        scatter (p in zip(combine_graphs.contigs_uid_vg, contigs_vcf_gz)) {
            call gbwt_index { input:
                vg = p.left,
                vcf_gz = p.right,
                vg_docker = vg_docker
            }
        }
        if (length(gbwt_index.gbwt) > 1) {
            call gbwt_merge { input:
                gbwts = gbwt_index.gbwt,
                graph_name = graph_name,
                vg_docker = vg_docker
            }
        }
        File? final_gbwt = if (defined(gbwt_merge.gbwt)) then gbwt_merge.gbwt else gbwt_index.gbwt[0]
    }

    # make xg index
    call xg_index { input:
        graph_name = graph_name,
        vg = combine_graphs.vg,
        vg_docker = vg_docker
    }

    # Prune the graph of repetitive sequences in preparation for GCSA indexing.
    # The workflow bifurcates for this, as the necessary invocations differ
    # significantly if we're using haplotypes or not.
    # If setting prune_options, make sure to set prune_graph.prune_options or
    # prune_graph_with_haplotypes.prune_options according to the desired path.
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
            contigs_gbwt = select_first([gbwt_index.gbwt]),
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
        File gcsa = gcsa_index.gcsa
        File gcsa_lcp = gcsa_index.lcp
        File? gbwt = final_gbwt
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
        set -exu -o pipefail
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
        set -exu -o pipefail
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

# construct the GBWT index of phased haplotypes
task gbwt_index {
    input {
        File vg
        File vcf_gz
        String vg_docker
    }

    command {
        set -exu -o pipefail
        nm=$(basename "~{vg}" .vg)
        tabix "~{vcf_gz}"

        vg index -G "$nm.gbwt" -v "~{vcf_gz}" "~{vg}"
    }

    output {
        File gbwt = glob("*.gbwt")[0]
    }

    runtime {
        docker: vg_docker
    }
}

# merge multiple GBWT indices (from disjoint contigs)
task gbwt_merge {
    input {
        Array[File]+ gbwts
        String graph_name
        String vg_docker
    }

    command {
        set -exu -o pipefail
        vg gbwt -m -f -o "~{graph_name}.gbwt" ~{sep=" " gbwts}
    }

    output {
        File gbwt = "~{graph_name}.gbwt"
    }

    runtime {
        docker: vg_docker
    }
}

# make xg index
task xg_index {
    input {
        String graph_name
        File vg
        String xg_options = ""
        String vg_docker
    }

    command {
        set -exu -o pipefail
        vg index -x "${graph_name}.xg" ~{xg_options} "~{vg}"
    }

    output {
        File xg ="${graph_name}.xg"
    }

    runtime {
        docker: vg_docker
    }
}

# prune repetitive/complex regions from graph (non-haplotype version)
task prune_graph {
    input {
        File contig_vg
        String prune_options = ""
        String vg_docker
    }

    command {
        set -exu -o pipefail
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

# prune repetitive/comlpex regions from graph (with haplotypes)
task prune_graph_with_haplotypes {
    input {
        Array[File]+ contigs_vg
        Array[File]+ contigs_gbwt
        File empty_id_map
        String prune_options = ""
        String vg_docker
    }

    command <<<
        set -exu -o pipefail
        cp "~{empty_id_map}" mapping
        paste -d ";" "~{write_lines(contigs_vg)}" "~{write_lines(contigs_gbwt)}" > inputs
        while IFS=';' read -ra p; do
            contig_vg="${p[0]}"
            contig_gbwt="${p[1]}"
            nm=$(basename "${contig_vg}" .vg)
            contig_pruned_vg="${nm}.pruned.vg"
            vg prune -u -g "$contig_gbwt" -a -m mapping "$contig_vg" ~{prune_options} > "$contig_pruned_vg"
            echo "$contig_pruned_vg" >> contigs_pruned_vg
        done < "inputs"
    >>>

    output {
        Array[File]+ contigs_pruned_vg=read_lines("contigs_pruned_vg")
    }

    runtime {
        docker: vg_docker
    }
}

# make GCSA index
task gcsa_index {
    input {
        String graph_name
        Array[File]+ contigs_pruned_vg
        File empty_id_map
        String gcsa_options = ""
        String vg_docker
    }

    command {
        set -exu -o pipefail
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
