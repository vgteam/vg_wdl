version 1.0

# This example uses vg_construct_and_index.wdl to build a graph including the
# ABO locus from GRCh38 (a 50Kbp region including the gene), and a known 5.8kbp
# deletion segregating in human populations. It then maps 1000 Genomes Project
# reads from a carrier of this deletion, verifying that vg maps some reads
# along the deletion edge.
#
#
# Recent study on this deletion allele:
#
# Moller M, Hellberg A, Olsson ML (2018). Thorough analysis of unorthodox ABO
#   deletions called by the 1000 Genomes Project. Vox Sang 113(2):185-197
#   doi:10.1111/vox.12613
#   https://www.ncbi.nlm.nih.gov/pubmed/29214632

import "../../workflows/vg_construct_and_index.wdl"
import "../../tasks/vg_map_hts.wdl"

workflow vg_ABOlocus_test {
    input {
        File ABOlocus_fa_gz
        File ABOlocus_SV_vcf_gz
        File reads_bam
        String vg_docker = "quay.io/vgteam/vg:v1.14.0"
    }

    # build & check the ABOlocus graph
    call vg_construct_and_index.vg_construct_and_index as cons { input:
        graph_name = "ABOlocus_SV",
        ref_fasta_gz = ABOlocus_fa_gz,
        contigs = ["ABOlocus"],
        contigs_vcf_gz = [ABOlocus_SV_vcf_gz],
        vg_docker = vg_docker
    }

    call check_graph { input:
        vg = cons.vg,
        vg_docker = vg_docker
    }

    # map reads to the ABOlocus graph & check the mappings
    call vg_map_hts.vg_map_hts as map { input:
        sam_bam_cram = reads_bam,
        xg = cons.xg,
        gcsa = cons.gcsa,
        gcsa_lcp = cons.gcsa_lcp,
        vg_docker = vg_docker
    }

    call check_gam { input:
        gam = map.gam,
        deletion_nodes = check_graph.deletion_nodes,
        vg_docker = vg_docker
    }

    output {
        Int deletion_edge_count = check_graph.deletion_edge_count
        Array[String] deletion_nodes = check_graph.deletion_nodes
        Int deletion_read_count = check_gam.deletion_read_count
    }
}

task check_graph {
    # checks the graph has the expected deletion edge
    input {
        File vg
        String vg_docker
    }

    command <<<
        set -ex -o pipefail
        apt-get update && apt-get install -y python3
        vg view "~{vg}" > graph.gfa
        python3 - <<EOF
        import sys
        node_sequences = {}
        deletion_edges = []
        with open("graph.gfa") as gfa:
            for line in gfa:
                line = line.strip().split("\t")
                if line[0] == "S":
                    node_sequences[line[1]] = line[2]
                if line[0] == "L" and int(line[3]) - int(line[1]) > 1:
                    deletion_edges.append(line)
        assert deletion_edges, "Graph has no deletion edges"
        for line in deletion_edges:
            assert node_sequences[line[1]].endswith("AAAT") and node_sequences[line[3]].startswith("CCCT")
        with open("deletion_edge_count","w") as outfile:
            print(str(len(deletion_edges)), file=outfile)
        with open("deletion_nodes","w") as outfile:
            for l in deletion_edges:
                print(l[1], file=outfile)
                print(l[3], file=outfile)
        EOF
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        Int deletion_edge_count = read_int("deletion_edge_count")
        Array[String] deletion_nodes = read_lines("deletion_nodes")
    }
}

task check_gam {
    # checks the gam has reads mapped along the deletion edge
    input {
        File gam
        Array[String] deletion_nodes
        String vg_docker
    }

    command <<<
        apt-get update && apt-get install -y python3
        vg view -a "~{gam}" -j > mappings.json
        python3 - <<EOF
        import json
        deletion_nodes = set()
        with open("~{write_lines(deletion_nodes)}") as infile:
            for line in infile:
                deletion_nodes.add(line.strip())
        deletion_reads = set()
        with open("mappings.json") as infile:
            for line in infile:
                line = json.loads(line)
                if "path" in line and "mapping" in line["path"]:
                    nodes = set()
                    for m in line["path"]["mapping"]:
                        nodes.add(m["position"]["node_id"])
                    if len(nodes.intersection(deletion_nodes)) > 1:
                        deletion_reads.add(line["name"])
        assert deletion_reads, "No reads mapped to deletion"
        with open("deletion_read_count", "w") as outfile:
            print(str(len(deletion_reads)), file=outfile)
        EOF
    >>>

    runtime {
        docker: vg_docker
    }

    output {
        Int deletion_read_count = read_int("deletion_read_count")
    }
}
