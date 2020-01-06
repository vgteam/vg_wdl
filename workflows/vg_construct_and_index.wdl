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
            "19", "20", "21", "22",  "X",  "Y", "MT"
        ]

        # VCF for each contig
        Array[File]+ contigs_vcf_gz

        # set true to GBWT index the VCF phased haplotypes
        Boolean use_haplotypes = false
        
        # set true to generate SNARLS index of VG graph
        Boolean make_snarls = false
        
        # set true to include decoy sequences from reference genome FASTA into VG graph construction
        Boolean use_decoys = false
        
        # regex to use in grep for extracting decoy contig names from reference FASTA
        String decoy_regex = ">GL\|>NC_007605\|>hs37d5"
        
        # vg docker image tag
        String vg_docker = "quay.io/vgteam/vg:v1.19.0"
    }

    # construct graph for each reference contig
    scatter (p in zip(contigs, contigs_vcf_gz)) {
        call construct_graph as construct_chromosome_graph { input:
            ref_fasta_gz = ref_fasta_gz,
            contig = p.left,
            vcf_gz = p.right,
            use_haplotypes = use_haplotypes,
            vg_docker = vg_docker,
            construct_cores = 2
        }
    }
    
    # extract decoy sequences from fasta and construct graph fro each decoy contig, if so configured
    if (use_decoys) {
        call extract_decoys { input:
            ref_fasta_gz = ref_fasta_gz,
            decoy_regex = decoy_regex,
            vg_docker = vg_docker
        }
        scatter (contig in extract_decoys.decoy_contig_ids) {
            call construct_graph as construct_decoy_graph { input:
                ref_fasta_gz = ref_fasta_gz,
                contig = contig,
                use_haplotypes = false,
                vg_docker = vg_docker,
                construct_cores = 2
            }
        }
        call concat as concat_vg_graph_lists { input:
            array_1 = construct_chromosome_graph.contig_vg,
            array_2 = construct_decoy_graph.contig_vg,
            vg_docker = vg_docker
        }
    }
    
    Array[File]? combined_contig_vg_list = if (use_decoys) then concat_vg_graph_lists.out else construct_chromosome_graph.contig_vg

    # combine them into a single graph with unique node IDs
    call combine_graphs { input:
        graph_name = graph_name,
        contigs_vg = select_first([combined_contig_vg_list]),
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
    
    # make snarls index
    if (make_snarls) {
        scatter (contig_vg in combine_graphs.all_contigs_uid_vg) {
            call snarls_index { input:
                vg = contig_vg,
                graph_name = graph_name,
                vg_docker = vg_docker
            }
        }
        call snarls_merge { input:
            snarls = snarls_index.snarls,
            graph_name = graph_name,
            vg_docker = vg_docker
        }
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
        scatter (contig_vg in combine_graphs.all_contigs_uid_vg) {
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
        # If decoys are to be included, prune those graphs and merge
        #  them with the haplotype pruned graphs
        if (use_decoys) {
            scatter (contig_vg in select_first([combine_graphs.decoy_contigs_uid_vg])) {
                call prune_graph as prune_decoy_graphs { input:
                    contig_vg = contig_vg,
                    vg_docker = vg_docker
                }
            }
            call concat as concat_pruned_vg_graph_lists { input:
                array_1 = prune_graph_with_haplotypes.contigs_pruned_vg,
                array_2 = prune_decoy_graphs.contig_pruned_vg,
                vg_docker = vg_docker
            }
        }
    }
    
    # make GCSA index
    call gcsa_index { input:
        graph_name = graph_name,
        contigs_pruned_vg = select_first([concat_pruned_vg_graph_lists.out, prune_graph.contig_pruned_vg, prune_graph_with_haplotypes.contigs_pruned_vg]),
        id_map = select_first([prune_graph_with_haplotypes.pruned_id_map,combine_graphs.empty_id_map]),
        vg_docker = vg_docker
    }

    output {
        File vg = combine_graphs.vg
        File xg = xg_index.xg
        File gcsa = gcsa_index.gcsa
        File gcsa_lcp = gcsa_index.lcp
        File? gbwt = final_gbwt
        File? snarls = snarls_merge.merged_snarls
    }
}

# extract decoy contigs from reference FASTA file
task extract_decoys {
    input {
        File ref_fasta_gz
        String decoy_regex
        String vg_docker
    }
    
    command <<<
        set -exu -o pipefail
        GREP_REGEX="~{decoy_regex}"
        zcat ~{ref_fasta_gz} | grep -E "${GREP_REGEX}" | cut -f 1 -d ' ' | cut -f 2 -d '>' >> decoy_contig_ids.txt
    >>>
    output {
        Array[String] decoy_contig_ids = read_lines("decoy_contig_ids.txt")
    }
    runtime {
        time: 10
        memory: 5
        docker: vg_docker
    }
}

# construct the graph for one reference contig
task construct_graph {
    input {
        File ref_fasta_gz
        File? vcf_gz
        String contig
        Boolean use_haplotypes
        String vg_construct_options="--node-max 32 --handle-sv"
        String vg_docker
        Int construct_cores
    }
    
    Boolean use_vcf = defined(vcf_gz)
    
    command <<<
        set -exu -o pipefail
        pigz -dc ~{ref_fasta_gz} > ref.fa
        VCF_OPTION_STRING=""
        if [ ~{use_vcf} == true ]; then
            tabix "~{vcf_gz}"
            VCF_OPTION_STRING="-v ~{vcf_gz} --region-is-chrom"
        fi

        vg construct --threads ~{construct_cores} -R "~{contig}" -C -r ref.fa ${VCF_OPTION_STRING} ~{vg_construct_options} ~{if use_haplotypes then "-a" else ""} > "~{contig}.vg" \
        && rm -f ref.fa
    >>>

    output {
        File contig_vg = "${contig}.vg"
    }

    runtime {
        time: 20
        cpu: construct_cores
        memory: 5 + " GB"
        disks: "local-disk 10 SSD"
        docker: vg_docker
    }
}

# helper task concatenates two Array[File] into a single Array[File]
#   borrowed from https://gatkforums.broadinstitute.org/wdl/discussion/8511/concatenating-arrays
task concat {
    input {
        Array[File]+ array_1
        Array[File]+ array_2
        String vg_docker
    }
    command {
        cat "~{write_lines(array_1)}"
        cat "~{write_lines(array_2)}"
    }
    output {
        Array[File]+ out = read_lines(stdout())
    }
    runtime {
        time: 5
        memory: 2 + " GB"
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
        touch decoy_contigs_uid_vg
        while read -r contig_vg; do
            nm=$(basename "$contig_vg")
            cp "$contig_vg" "vg/$nm"
            if [[ $nm == *"GL"* || $nm == *"NC_007605"* || $nm == *"hs37d5"* ]]; then
                echo "vg/$nm" >> decoy_contigs_uid_vg
            else
                echo "vg/$nm" >> contigs_uid_vg
            fi
            echo "vg/$nm" >> all_contigs_uid_vg
        done < "~{write_lines(contigs_vg)}"
        xargs -n 999999 vg ids -j -m empty.id_map < all_contigs_uid_vg
        mkdir concat
        xargs -n 999999 cat < all_contigs_uid_vg > "concat/${graph_name}.vg"
    }

    output {
        File vg = "concat/${graph_name}.vg"
        File empty_id_map = "empty.id_map"
        Array[File]+ contigs_uid_vg = read_lines("contigs_uid_vg")
        Array[File]+ all_contigs_uid_vg = read_lines("all_contigs_uid_vg")
        Array[File]? decoy_contigs_uid_vg = read_lines("decoy_contigs_uid_vg")
    }

    runtime {
        time: 180
        memory: 20 + " GB"
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

    command <<<
        set -exu -o pipefail
        nm=$(basename "~{vg}" .vg)
        tabix "~{vcf_gz}"

        vg index --threads "$(nproc --all)" -G "$nm.gbwt" -v "~{vcf_gz}" "~{vg}"
    >>>

    output {
        File gbwt = glob("*.gbwt")[0]
    }

    runtime {
        time: 30
        cpu: 30
        memory: 15 + " GB"
        disks: "local-disk 10 SSD"
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
        time: 10
        memory: 15 + " GB"
        disks: "local-disk 10 SSD"
        docker: vg_docker
    }
}

# construct the Snarls index from the VG graph
task snarls_index {
    input {
        File vg
        String graph_name
        String vg_docker
    }
    
    command {
        set -exu -o pipefail
        nm=$(basename "~{vg}" .vg)
        vg snarls -t ~{vg} > "$nm.snarls"
    }
    
    output {
        File snarls = glob("*.snarls")[0]
    }
     
    runtime {
        memory: 40 + " GB"
        disks: "local-disk 100 SSD"
        docker: vg_docker
    }
}

# merge multiple SNARL indices (from disjoint contigs)
task snarls_merge {
    input {
        Array[File]+ snarls
        String graph_name
        String vg_docker
    }
    
    command {
        set -exu -o pipefail
        cat ${sep=" " snarls} > "~{graph_name}.snarls"
    }
    
    output {
        File merged_snarls = "~{graph_name}.snarls"
    }
     
    runtime {
        disks: "local-disk 20 SSD"
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

    command <<<
        set -exu -o pipefail
        vg index --threads "$(nproc --all)" -x "~{graph_name}.xg" ~{xg_options} "~{vg}"
    >>>

    output {
        File xg ="~{graph_name}.xg"
    }

    runtime {
        time: 240
        cpu: 32
        memory: 50 + " GB"
        disks: "local-disk 10 SSD"
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
        vg prune --threads 2 -r "${contig_vg}" ~{prune_options} > "$nm.pruned.vg"
    }

    output {
        File contig_pruned_vg = glob("*.pruned.vg")[0]
    }

    runtime {
        time: 180
        cpu: 2
        memory: 20 + " GB"
        disks: "local-disk 10 SSD"
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
            vg prune --threads 2 -u -g "$contig_gbwt" -a -m mapping "$contig_vg" ~{prune_options} > "$contig_pruned_vg"
            echo "$contig_pruned_vg" >> contigs_pruned_vg
        done < "inputs"
    >>>

    output {
        Array[File]+ contigs_pruned_vg=read_lines("contigs_pruned_vg")
        File pruned_id_map = "mapping"
    }

    runtime {
        time: 180
        cpu: 2
        memory: 20 + " GB"
        disks: "local-disk 10 SSD"
        docker: vg_docker
    }
}

# make GCSA index
task gcsa_index {
    input {
        String graph_name
        Array[File]+ contigs_pruned_vg
        File id_map
        String gcsa_options = ""
        String vg_docker
    }

    command <<<
        set -exu -o pipefail
        vg index --threads 32 -p -g "~{graph_name}.gcsa" -f "~{id_map}" ~{gcsa_options} ~{sep=" " contigs_pruned_vg}
    >>>

    output {
        File gcsa = "~{graph_name}.gcsa"
        File lcp = "~{graph_name}.gcsa.lcp"
    }

    runtime {
        time: 1200
        cpu: 32
        memory: 100 + " GB"
        disks: "local-disk 50 SSD"
        docker: vg_docker
    }
}
