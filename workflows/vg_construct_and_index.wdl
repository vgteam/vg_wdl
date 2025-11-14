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

        # set true to build giraffe-required indexes
        Boolean giraffe_indexes = false
        
        # set true to include decoy sequences from reference genome FASTA into VG graph construction
        Boolean use_decoys = false
        
        # set true to include structural variants from input vcfs into VG graph construction
        Boolean use_svs = false
        
        # Set to 'true' to use small resources for tiny test dataset
        Boolean in_small_resources = false
         
        # regex to use in grep for extracting decoy contig names from reference FASTA
        String decoy_regex = ">GL\|>NC_007605\|>hs37d5\|>hs38d1_decoys\|>chrEBV\|>chrUn\|>chr\([1-2][1-9]\|[1-9]\|Y\)_"
        
        # vg docker image tag
        String vg_docker = "quay.io/vgteam/vg:v1.64.0"
    }

    # construct graph for each reference contig
    scatter (p in zip(contigs, contigs_vcf_gz)) {
        call construct_graph as construct_chromosome_graph { input:
            ref_fasta_gz = ref_fasta_gz,
            contig = p.left,
            vcf_gz = p.right,
            use_haplotypes = giraffe_indexes,
            use_svs = use_svs,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
    }
    
    # extract decoy sequences from fasta and construct graph fro each decoy contig, if so configured
    if (use_decoys) {
        call extract_decoys { input:
            ref_fasta_gz = ref_fasta_gz,
            decoy_regex = decoy_regex,
            in_small_resources = in_small_resources
        }
        scatter (contig in extract_decoys.decoy_contig_ids) {
            call construct_graph as construct_decoy_graph { input:
                ref_fasta_gz = ref_fasta_gz,
                contig = contig,
                use_haplotypes = false,
                use_svs = false,
                vg_docker = vg_docker,
                in_small_resources = in_small_resources
            }
        }
    }
    

    # combine them into a single graph with unique node IDs
    call combine_graphs { input:
        graph_name = graph_name,
        contigs_vg = select_first([construct_chromosome_graph.contig_vg]),
        decoy_contigs_vg = construct_decoy_graph.contig_vg,
        vg_docker = vg_docker,
        in_small_resources = in_small_resources
    }

    # make GBWT index, if so configured
    if (giraffe_indexes) {
        scatter (p in zip(combine_graphs.contigs_uid_vg, contigs_vcf_gz)) {
            call gbwt_index { input:
                vg = p.left,
                vcf_gz = p.right,
                vg_docker = vg_docker,
                in_small_resources = in_small_resources
            }
        }
        if (length(gbwt_index.gbwt) > 1) {
            call gbwt_merge { input:
                gbwts = gbwt_index.gbwt,
                graph_name = graph_name,
                vg_docker = vg_docker,
                in_small_resources = in_small_resources
            }
        }
        File? final_gbwt = if (defined(gbwt_merge.gbwt)) then gbwt_merge.gbwt else gbwt_index.gbwt[0]
    }
    
    # make xg index
    call xg_index { input:
        graph_name = graph_name,
        vg = combine_graphs.vg,
        vg_docker = vg_docker,
        in_small_resources = in_small_resources
    }

    # Prune the graph of repetitive sequences in preparation for GCSA indexing.
    # The workflow bifurcates for this, as the necessary invocations differ
    # significantly if we're using haplotypes or not.
    # If setting prune_options, make sure to set prune_graph.prune_options or
    # prune_graph_with_haplotypes.prune_options according to the desired path.
    if (!giraffe_indexes) {
        scatter (contig_vg in combine_graphs.all_contigs_uid_vg) {
            call prune_graph { input:
                contig_vg = contig_vg,
                vg_docker = vg_docker,
                in_small_resources = in_small_resources
            }
        }
    }
    if (giraffe_indexes) {
        call prune_graph_with_haplotypes { input:
            contigs_vg = combine_graphs.contigs_uid_vg,
            contigs_gbwt = select_first([gbwt_index.gbwt]),
            empty_id_map = combine_graphs.empty_id_map,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
        # If decoys are to be included, prune those graphs and merge
        #  them with the haplotype pruned graphs
        if (use_decoys) {
            scatter (contig_vg in select_first([combine_graphs.decoy_contigs_uid_vg])) {
                call prune_graph as prune_decoy_graphs { input:
                    contig_vg = contig_vg,
                    vg_docker = vg_docker,
                    in_small_resources = in_small_resources
                }
            }
            Array[File] concat_pruned_vg_graph_lists = flatten([select_first([prune_graph_with_haplotypes.contigs_pruned_vg,[]]), select_first([prune_decoy_graphs.contig_pruned_vg,[]])]) 
        }
    }
    
    if (giraffe_indexes) {
        # make Trivial Snarls index
        scatter (contig_vg in combine_graphs.all_contigs_uid_vg) {
            call snarls_index { input:
                vg = contig_vg,
                graph_name = graph_name,
                vg_docker = vg_docker,
                in_small_resources = in_small_resources
            }
        }
        call snarls_merge { input:
            snarls = snarls_index.snarls,
            graph_name = graph_name,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
        # make Distance index
        call dist_index { input:
            graph_name = graph_name,
            xg = xg_index.xg,
            trivial_snarls = snarls_merge.merged_snarls,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
        # make Graph GBWT index
        call sampled_gbwt_index { input:
            graph_name = graph_name,
            gbwt = final_gbwt,
            xg = xg_index.xg,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
        # make Minimizer index
        call min_index { input:
            graph_name = graph_name,
            sampled_gbwt = sampled_gbwt_index.sampled_gbwt,
            sampled_gg = sampled_gbwt_index.sampled_gg,
            dist = dist_index.dist,
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
    }
    if (!giraffe_indexes) {
        # make GCSA index
        call gcsa_index { input:
            graph_name = graph_name,
            contigs_pruned_vg = select_first([concat_pruned_vg_graph_lists, prune_graph.contig_pruned_vg, prune_graph_with_haplotypes.contigs_pruned_vg]),
            id_map = select_first([prune_graph_with_haplotypes.pruned_id_map,combine_graphs.empty_id_map]),
            vg_docker = vg_docker,
            in_small_resources = in_small_resources
        }
    }

    output {
        File vg = combine_graphs.vg
        File xg = xg_index.xg
        File? gcsa = gcsa_index.gcsa
        File? gcsa_lcp = gcsa_index.lcp
        File? gbwt = select_first([sampled_gbwt_index.sampled_gbwt, final_gbwt])
        File? ggbwt = sampled_gbwt_index.sampled_gg
        File? snarls = snarls_merge.merged_snarls
        File? dist = dist_index.dist
        File? min = min_index.min
    }
}

# extract decoy contigs from reference FASTA file
task extract_decoys {
    input {
        File ref_fasta_gz
        String decoy_regex
        Boolean in_small_resources
    }

    String in_mem = if in_small_resources then "1" else "5"

    command <<<
        set -exu -o pipefail
        GREP_REGEX="~{decoy_regex}"
        zcat ~{ref_fasta_gz} | grep -E "${GREP_REGEX}" | cut -f 1 -d ' ' | cut -f 2 -d '>' >> decoy_contig_ids.txt
    >>>
    output {
        Array[String] decoy_contig_ids = read_lines("decoy_contig_ids.txt")
    }
    runtime {
        preemptible: 2
        time: 10
        memory: in_mem + " GB"
        docker: "ubuntu:24.04"
    }
}

# construct the graph for one reference contig
task construct_graph {
    input {
        File ref_fasta_gz
        File? vcf_gz
        String contig
        Boolean use_haplotypes
        Boolean use_svs
        String vg_construct_options="--node-max 32"
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 20 else 20
    String in_mem = if in_small_resources then "20" else "50"
 
    Boolean use_vcf = defined(vcf_gz)
    
    command <<<
        set -exu -o pipefail
        pigz -dc ~{ref_fasta_gz} > ref.fa
        VCF_OPTION_STRING=""
        if [ ~{use_vcf} == true ]; then
            tabix "~{vcf_gz}"
            VCF_OPTION_STRING="-v ~{vcf_gz} --region-is-chrom"
        fi
        
        VG_SV_OPTION=""
        if [ ~{use_svs} == true ]; then
            VG_SV_OPTION="--handle-sv"
        fi

        vg construct --threads ~{in_cores} -R "~{contig}" -C -r ref.fa ${VCF_OPTION_STRING} ~{vg_construct_options} ${VG_SV_OPTION} ~{if use_haplotypes then "-a" else ""} > "~{contig}.vg" \
        && rm -f ref.fa
    >>>

    output {
        File contig_vg = "~{contig}.vg"
    }

    runtime {
        preemptible: 2
        time: 200
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
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
        preemptible: 2
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
        Array[File]? decoy_contigs_vg
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 2 else 50
    String in_mem = if in_small_resources then "20" else "40"
    
    Boolean decoy_contigs_exist = defined(decoy_contigs_vg)
    Array[File] decoy_contigs_vg_resolved = select_first([decoy_contigs_vg, []])
    
    command {
        set -exu -o pipefail
        # we approach this in a particular way to ensure the output array contigs_uid_vg has the
        # same order as the input array contigs_vg (so we can't rely on glob patterns)
        mkdir vg_decoy_contigs/ vg_contigs/ vg_all_contigs/
        while read -r contig_vg; do
            nm=$(basename "$contig_vg")
            if [[ $nm == *"GL"* || $nm == *"NC_007605"* || $nm == *"hs37d5"* || $nm == *"KI"* || $nm == *"chrEBV"* || $nm == *"chrUn"* || $nm == *"hs38d1_decoys"* ]]; then
                cp "$contig_vg" "vg_decoy_contigs/$nm"
                echo "vg_decoy_contigs/$nm" >> all_contigs_uid_vg
            else
                cp "$contig_vg" "vg_contigs/$nm"
                echo "vg_contigs/$nm" >> all_contigs_uid_vg
            fi
        done < "~{write_lines(contigs_vg)}"
        
        if [ ~{decoy_contigs_exist} == true ]; then
            while read -r decoy_contig_vg; do
                nm=$(basename "$decoy_contig_vg")
                if [[ $nm == *"GL"* || $nm == *"NC_007605"* || $nm == *"hs37d5"* || $nm == *"KI"* || $nm == *"chrEBV"* || $nm == *"chrUn"* || $nm == *"hs38d1_decoys"* ]]; then
                    cp "$decoy_contig_vg" "vg_decoy_contigs/$nm"
                    echo "vg_decoy_contigs/$nm" >> all_contigs_uid_vg
                else
                    cp "$decoy_contig_vg" "vg_contigs/$nm"
                    echo "vg_contigs/$nm" >> all_contigs_uid_vg
                fi
            done < "~{write_lines(decoy_contigs_vg_resolved)}"
        fi
        xargs -n 999999 vg ids -j -m empty.id_map < all_contigs_uid_vg
        mkdir concat
        xargs -n 999999 vg combine < all_contigs_uid_vg > "concat/~{graph_name}.vg"
    }

    output {
        File vg = "concat/~{graph_name}.vg"
        File empty_id_map = "empty.id_map"
        Array[File]+ contigs_uid_vg = glob("vg_contigs/*.vg")
        Array[File]? decoy_contigs_uid_vg = glob("vg_decoy_contigs/*.vg")
        Array[File]+ all_contigs_uid_vg = flatten([contigs_uid_vg, select_first([decoy_contigs_uid_vg,[]])])
    }

    runtime {
        preemptible: 2
        time: 800
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# construct the GBWT index of phased haplotypes
task gbwt_index {
    input {
        File vg
        File vcf_gz
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 5 else 20
    String in_mem = if in_small_resources then "5" else "50"

    command <<<
        set -exu -o pipefail
        nm=$(basename "~{vg}" .vg)
        tabix "~{vcf_gz}"

        vg index --threads ~{in_cores} --force-phasing --discard-overlaps -G "$nm.gbwt" -v "~{vcf_gz}" "~{vg}"
    >>>

    output {
        File gbwt = glob("*.gbwt")[0]
    }

    runtime {
        preemptible: 2
        time: 800
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# merge multiple GBWT indices (from disjoint contigs)
task gbwt_merge {
    input {
        Array[File]+ gbwts
        String graph_name
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 2 else 10
    String in_mem = if in_small_resources then "5" else "50"

    command {
        set -exu -o pipefail
        vg gbwt -m -f -o "~{graph_name}.gbwt" ~{sep=" " gbwts}
    }

    output {
        File gbwt = "~{graph_name}.gbwt"
    }

    runtime {
        preemptible: 2
        time: 100
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# construct the Snarls index from the VG graph
task snarls_index {
    input {
        File vg
        String graph_name
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 1 else 1
    Int in_disk = if in_small_resources then 1 else 100
    String in_mem = if in_small_resources then "10" else "120"
    
    command {
        set -exu -o pipefail
        nm=$(basename "~{vg}" .vg)
        vg snarls -t ~{in_cores} --include-trivial ~{vg} > "$nm.snarls"
    }
    
    output {
        File snarls = glob("*.snarls")[0]
    }
     
    runtime {
        preemptible: 2
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# merge multiple SNARL indices (from disjoint contigs)
task snarls_merge {
    input {
        Array[File]+ snarls
        String graph_name
        String vg_docker
        Boolean in_small_resources
    }

    Int in_disk = if in_small_resources then 1 else 20
    
    command {
        set -exu -o pipefail
        cat ${sep=" " snarls} > "~{graph_name}.snarls"
    }
    
    output {
        File merged_snarls = "~{graph_name}.snarls"
    }
     
    runtime {
        preemptible: 2
        disks: "local-disk " + in_disk + " SSD"
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
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 8 else 16
    Int in_disk = if in_small_resources then 2 else 60
    String in_mem = if in_small_resources then "20" else "60"

    command <<<
        set -exu -o pipefail
        vg index --threads ~{in_cores} -x "~{graph_name}.xg" ~{xg_options} "~{vg}"
    >>>

    output {
        File xg ="~{graph_name}.xg"
    }

    runtime {
        preemptible: 2
        time: 240
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# Make the distance index
task dist_index { 
    input {
        String graph_name
        File xg
        File trivial_snarls
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 4
    Int in_disk = if in_small_resources then 2 else 50
    String in_mem = if in_small_resources then "10" else "100"
    
    command <<<
        set -exu -o pipefail
        vg index --threads ~{in_cores} -x ~{xg} -s ~{trivial_snarls} -j "~{graph_name}.dist"
    >>>
    output {
        File dist = "~{graph_name}.dist"
    }
    runtime {
        preemptible: 2
        time: 240
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# Make the sampled gbwt and graph gbwt index
task sampled_gbwt_index { 
    input {
        String graph_name
        File? gbwt
        File xg
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 20 else 80
    String in_mem = if in_small_resources then "20" else "80"
    
    command <<<
        set -exu -o pipefail
        vg gbwt -l -n 64 --num-threads ~{in_cores} ~{gbwt} -x ~{xg} -o "~{graph_name}.sampled.gbwt" -g "~{graph_name}.sampled.gg"
    >>>
    output {
        File sampled_gbwt = "~{graph_name}.sampled.gbwt"
        File sampled_gg = "~{graph_name}.sampled.gg"
    }
    runtime {
        preemptible: 2
        time: 240
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# Make the minimizer index
task min_index { 
    input {
        String graph_name
        File sampled_gbwt
        File sampled_gg
        File dist
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 4 else 100
    String in_mem = if in_small_resources then "5" else "80"
    
    command <<<
        set -exu -o pipefail
        vg minimizer --threads ~{in_cores} -g ~{sampled_gbwt} -G ~{sampled_gg} -d ~{dist} -i "~{graph_name}.min"
    >>>
    output {
        File min = "~{graph_name}.min"
    }
    runtime {
        preemptible: 2
        time: 240
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}

# prune repetitive/complex regions from graph (non-haplotype version)
task prune_graph {
    input {
        File contig_vg
        String prune_options = ""
        String vg_docker
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 1 else 10
    String in_mem = if in_small_resources then "10" else "20"

    command {
        set -exu -o pipefail
        nm=$(basename "${contig_vg}" .vg)
        vg prune --threads ~{in_cores} -r "${contig_vg}" ~{prune_options} > "$nm.pruned.vg"
    }

    output {
        File contig_pruned_vg = glob("*.pruned.vg")[0]
    }

    runtime {
        preemptible: 2
        time: 180
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
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
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 4 else 50
    String in_mem = if in_small_resources then "5" else "20"

    command <<<
        set -exu -o pipefail
        cp "~{empty_id_map}" mapping
        paste -d ";" "~{write_lines(contigs_vg)}" "~{write_lines(contigs_gbwt)}" > inputs
        mkdir outdir
        while IFS=';' read -ra p; do
            contig_vg="${p[0]}"
            contig_gbwt="${p[1]}"
            nm=$(basename "${contig_vg}" .vg)
            contig_pruned_vg="${nm}.pruned.vg"
            vg prune --threads ~{in_cores} -u -g "$contig_gbwt" -a -m mapping "$contig_vg" ~{prune_options} > "outdir/$contig_pruned_vg"
        done < "inputs"
    >>>

    output {
        Array[File]+ contigs_pruned_vg = glob("outdir/*.vg") 
        File pruned_id_map = "mapping"
    }

    runtime {
        preemptible: 2
        time: 180
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
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
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 8 else 32
    Int in_disk = if in_small_resources then 4 else 50
    String in_mem = if in_small_resources then "4" else "250"

    command <<<
        set -exu -o pipefail
        vg index --threads ~{in_cores} -p -g "~{graph_name}.gcsa" -f "~{id_map}" ~{gcsa_options} ~{sep=" " contigs_pruned_vg}
    >>>

    output {
        File gcsa = "~{graph_name}.gcsa"
        File lcp = "~{graph_name}.gcsa.lcp"
    }

    runtime {
        preemptible: 2
        time: 1200
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: vg_docker
    }
}
