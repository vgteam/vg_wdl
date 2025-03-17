version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/vg_indexing.wdl" as index

workflow HaplotypeSampling {
    meta {
        author: "Parsa Eskandar"
        email: "seeskand@ucsc.edu"
        description: "Create haplotype sampled graph and the indexes necessary for the vg giraffe. More information at https://github.com/vgteam/vg/wiki/Haplotype-Sampling"
    }

    parameter_meta {
        GBZ_FILE: "Path to .gbz index file"
        INPUT_READ_FILE_FIRST: "Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_SECOND: "Input sample 2st read pair fastq.gz or fastq"
        HAPL_FILE: "Path to .hapl file"
        DIST_FILE: "Path to .dist file"
        R_INDEX_FILE: "Path to .ri file"
        KFF_FILE: "Path to .kff file"
        OUTPUT_NAME_PREFIX: "Name of the output file (Default: haplotype_sampled_graph)"
        KMER_LENGTH: "Size of kmer using for sampling (Up to 31) (Default: 29)"
        CORES: "Number of cores to use with commands. (Default: 16)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320)"

        WINDOW_LENGTH: "Window length used for building the minimizer index for sampling haplotypes. (Default: 11)"
        SUBCHAIN_LENGTH: "Target length (in bp) for subchains. (Default: 10000)"
        HAPLOTYPE_NUMBER: "Number of generated synthetic haplotypes. (Default: 4)"
        PRESENT_DISCOUNT: "Multiplicative factor for discounting scores for present kmers. (Default: 0.9)"
        HET_ADJUST: "Additive term for adjusting scores for heterozygous kmers. (Default: 0.05)"
        ABSENT_SCORE: "Score for absent kmers. (Default: 0.8)"
        INCLUDE_REFERENCE: "Include reference paths and generic paths from the full graph in the sampled graph. (Default: true)"
        SET_REFERENCE: "Name of single reference to include in sampled graph. (Default: all references)"
        DIPLOID: "Activate diploid sampling. (Default: true)"
        INDEX_MINIMIZER_K: "K-mer size of minimizer index to produce for sampled graph. Should be 29 for short read mapping and 31 for long read mapping. (Default: 29)"
        INDEX_MINIMIZER_W: "Window size of minimizer index to produce for sampled graph. Should be 11 for short read mapping and 50 for long read mapping. (Default: 11)"
        INDEX_MINIMIZER_WEIGHTED: "Whether to produce a weighted minimizer index for the sampled graph. (Default: true)"
        VG_DOCKER: "Container image to use when running vg."
    }
    input {
        File GBZ_FILE
        File INPUT_READ_FILE_FIRST
        File? INPUT_READ_FILE_SECOND
        File? HAPL_FILE
        File? DIST_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        String OUTPUT_NAME_PREFIX = "haplotype_sampled_graph"
        Int KMER_LENGTH = 29
        Int CORES = 16
        Int INDEX_MINIMIZER_MEM = 320
        Int WINDOW_LENGTH = 11
        Int SUBCHAIN_LENGTH = 10000
        Int HAPLOTYPE_NUMBER = 4
        Float PRESENT_DISCOUNT = 0.9
        Float HET_ADJUST = 0.05
        Float ABSENT_SCORE = 0.8
        Boolean INCLUDE_REFERENCE = true
        String? SET_REFERENCE
        Boolean DIPLOID = true
        Int INDEX_MINIMIZER_K = 29
        Int INDEX_MINIMIZER_W = 11
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        String? VG_DOCKER


    }

    # Have to create haplotype information
    if (!defined(HAPL_FILE)) {
        # create the dist index file and r-index file to create the haplotype information file .hapl

        if (!defined(DIST_FILE)) {
            call index.createDistanceIndex {
                input:
                    in_gbz_file=GBZ_FILE
            }
        }

        File dist_index_file = select_first([DIST_FILE, createDistanceIndex.output_dist_index])

        if (!defined(R_INDEX_FILE)) {
            call index.createRIndex {
                input:
                    in_gbz_file=GBZ_FILE,
                    nb_cores=CORES
            }
        }

        File r_index_file = select_first([R_INDEX_FILE, createRIndex.output_R_index])

        # create the haplotype information file
        call index.createHaplotypeIndex {
            input:
                in_gbz_file=GBZ_FILE,
                in_dist_index=dist_index_file,
                in_R_index=r_index_file,
                nb_cores=CORES,
                kmer_length=KMER_LENGTH,
                window_length=WINDOW_LENGTH,
                subchain_length=SUBCHAIN_LENGTH
        }
    }

    File haplotype_index = select_first([HAPL_FILE, createHaplotypeIndex.output_hap_index])

    if (!defined(KFF_FILE)) {
        call utils.kmerCountingKMC {
            input:
                input_read_file_1=INPUT_READ_FILE_FIRST,
                input_read_file_2=INPUT_READ_FILE_SECOND,
                output_file_name=OUTPUT_NAME_PREFIX,
                kmer_length=KMER_LENGTH,
                nb_cores=CORES

        }

    }

    File kmer_information = select_first([KFF_FILE, kmerCountingKMC.kff_file])

    call map.samplingHaplotypes {
        input:
            in_gbz_file=GBZ_FILE,
            in_hap_index=haplotype_index,
            in_kmer_info=kmer_information,
            output_file_name=OUTPUT_NAME_PREFIX,
            haplotype_number = HAPLOTYPE_NUMBER,
            present_discount = PRESENT_DISCOUNT,
            het_adjust = HET_ADJUST,
            absent_score = ABSENT_SCORE,
            include_reference = INCLUDE_REFERENCE,
            set_reference = SET_REFERENCE,
            nb_cores=CORES,
            use_diploid_sampling=DIPLOID


    }

    call index.createDistanceIndex as giraffeDist {
                input:
                    in_gbz_file=samplingHaplotypes.output_graph
    }

    call index.createMinimizerIndex {
        input:
            in_gbz_file=samplingHaplotypes.output_graph,
            in_dist_index=giraffeDist.output_dist_index,
            in_minimizer_k = INDEX_MINIMIZER_K,
            in_minimizer_w = INDEX_MINIMIZER_W,
            in_minimizer_weighted = INDEX_MINIMIZER_WEIGHTED,
            out_name=OUTPUT_NAME_PREFIX,
            nb_cores=CORES,
            in_extract_mem=INDEX_MINIMIZER_MEM
    }

    output {
        File sampled_graph = samplingHaplotypes.output_graph
        File sampled_min = createMinimizerIndex.output_minimizer
        File sampled_zipcodes = createMinimizerIndex.output_zipcodes
        File sampled_dist = giraffeDist.output_dist_index
    }

}






