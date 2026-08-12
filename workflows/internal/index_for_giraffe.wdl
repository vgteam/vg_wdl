version 1.0

import "../../tasks/vg_indexing.wdl" as index
import "../haplotype_sampling.wdl" as hapl

workflow IndexForGiraffe {
    meta {
        description: "## Giraffe indexes workflow \n Produce the set of indexes that `vg giraffe` needs for one graph: a GBZ and the matching distance, minimizer and zipcodes indexes. Indexes that are passed in are used as they are, and the rest are built with the given vg container. When haplotype sampling is on, the graph is replaced by a sample-specific sampled graph and all of the mapping indexes are rebuilt for it, so the reads to sample from are required."
    }

    parameter_meta {
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "(OPTIONAL) Path to .dist index file. Built if not provided."
        MIN_FILE: "(OPTIONAL) Path to .min index file. Built if not provided, along with a new .zipcodes file."
        ZIPCODES_FILE: "(OPTIONAL) Path to .zipcodes index file matching MIN_FILE. Ignored unless MIN_FILE is also provided."
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        R_INDEX_FILE: "(OPTIONAL) Path to .ri file used in haplotype sampling"
        KFF_FILE: "(OPTIONAL) Path to .kff file of sample kmer counts used in haplotype sampling"
        HAPLOTYPE_SAMPLING: "Whether to haplotype sample the graph down to the sample's haplotypes. Default is 'true'."
        INPUT_READ_FILE_FIRST: "Reads to count kmers in for haplotype sampling. Required if HAPLOTYPE_SAMPLING is true and KFF_FILE is not given."
        INPUT_READ_FILE_SECOND: "(OPTIONAL) Second read pair file to count kmers in for haplotype sampling."
        GIRAFFE_PRESET: "Name of the Giraffe mapper parameter preset the indexes are for (default, fast, hifi, or r10). This sets the minimizer size and window size."
        HAPLOTYPE_NUMBER: "Number of generated synthetic haplotypes used in haplotype sampling. (Default: 32)"
        DIPLOID: "Whether to use diploid sampling while doing haplotype sampling. (Default: true)"
        SET_REFERENCE: "(OPTIONAL) Name of the single reference to keep for haplotype sampling."
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for indexing tasks (distance index, r-index, haplotype index, and sampling). (Default: 120)"
        CORES: "Number of cores to use for indexing. (Default: 16)"
        VG_DOCKER: "Container image to use when running vg"
    }

    input {
        File GBZ_FILE
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        File? HAPL_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        Boolean HAPLOTYPE_SAMPLING = true
        File? INPUT_READ_FILE_FIRST
        File? INPUT_READ_FILE_SECOND
        String GIRAFFE_PRESET = "default"
        Int HAPLOTYPE_NUMBER = 32
        Boolean DIPLOID = true
        String? SET_REFERENCE
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int HAPLOTYPE_INDEXING_MEM = 120
        Int CORES = 16
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
    }

    # Short reads want short minimizers in tight windows; the long read presets
    # want the opposite. These have to match what the mapper expects.
    Int minimizer_k = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 29 else 31
    Int minimizer_w = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 11 else 50

    if (HAPLOTYPE_SAMPLING) {
        call hapl.HaplotypeSampling {
            input:
            GBZ_FILE=GBZ_FILE,
            INPUT_READ_FILE_FIRST=select_first([INPUT_READ_FILE_FIRST]),
            INPUT_READ_FILE_SECOND=INPUT_READ_FILE_SECOND,
            HAPL_FILE=HAPL_FILE,
            DIST_FILE=DIST_FILE,
            R_INDEX_FILE=R_INDEX_FILE,
            KFF_FILE=KFF_FILE,
            HAPLOTYPE_NUMBER=HAPLOTYPE_NUMBER,
            DIPLOID=DIPLOID,
            SET_REFERENCE=SET_REFERENCE,
            CORES=CORES,
            KMER_COUNTING_MEM=KMER_COUNTING_MEM,
            HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
            VG_DOCKER=VG_DOCKER
        }
    }

    File file_gbz = select_first([HaplotypeSampling.sampled_graph, GBZ_FILE])

    # DIST_FILE and MIN_FILE (if given) are indexes of GBZ_FILE. If haplotype
    # sampling ran, file_gbz is a different graph, so they can't be reused
    # and have to be rebuilt.
    if (HAPLOTYPE_SAMPLING || !defined(DIST_FILE)) {
        call index.createDistanceIndex {
            input:
                in_gbz_file=file_gbz,
                nb_cores=CORES,
                in_extract_mem=HAPLOTYPE_INDEXING_MEM,
                vg_docker=VG_DOCKER
        }
    }
    File file_dist = select_first([createDistanceIndex.output_dist_index, DIST_FILE])

    if (HAPLOTYPE_SAMPLING || !defined(MIN_FILE)) {
        # A minimizer index and its zipcodes are made together and only make
        # sense together, so any provided zipcodes are dropped here.
        call index.createMinimizerIndex {
            input:
                in_gbz_file=file_gbz,
                in_dist_index=file_dist,
                in_minimizer_k=minimizer_k,
                in_minimizer_w=minimizer_w,
                in_minimizer_weighted=INDEX_MINIMIZER_WEIGHTED,
                out_name=sub(basename(GBZ_FILE), "\\.gbz$", ""),
                nb_cores=CORES,
                in_extract_mem=INDEX_MINIMIZER_MEM,
                vg_docker=VG_DOCKER
        }
    }
    File file_min = select_first([createMinimizerIndex.output_minimizer, MIN_FILE])

    # The zipcodes are optional all the way through, so we can't select_first on
    # them; every candidate might be null.
    # The zipcodes might be missing if we got MIN_FILE but no ZIPCODES_FILE.
    # This is an increasingly bad idea, but still allowed for now.
    Array[File] possible_zipcode_files = select_all([createMinimizerIndex.output_zipcodes, ZIPCODES_FILE])
    # WDL 1.0 has no None literal, so we make a File? that is never assigned.
    if (false) {
        Array[File] no_files = []
        #@ except: UnnecessaryFunctionCall
        File NULL_FILE = select_first(no_files)
    }

    output {
        File gbz_file = file_gbz
        File dist_file = file_dist
        File min_file = file_min
        File? zipcodes_file = if length(possible_zipcode_files) > 0 then possible_zipcode_files[0] else NULL_FILE
    }
}
