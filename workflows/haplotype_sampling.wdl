version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/vg_map_hts.wdl" as map

workflow HaplotypeSampling {
    meta {
        author: "Parsa Eskandar"
        email: "seeskand@ucsc.edu"
        description: "Create haplotype sampled graph. More information at https://github.com/vgteam/vg/wiki/Haplotype-Sampling"
    }

    parameter_meta {
        GBZ_FILE: "Path to .gbz index file"
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz"
        INPUT_READ_FILE_2: "Input sample 2st read pair fastq.gz"
        HAPL_FILE: "Path to .hapl file"
        DIST_FILE: "Path to .dist file"
        R_INDEX_FILE: "Path to .ri file"
        KFF_FILE: "Path to .kff file"
        IN_OUTPUT_NAME_PREFIX: "Name of the output file (default: haplotype_sampled_graph)"
        IN_KMER_LENGTH: "Size of kmer using for sampling (default: 29)"
        IN_WORKING_DIRECTORY: "Path to a directory that files are written to (default: '.')"
        CORES: "Number of cores to use with commands. (default: 16)"

    }
    input {
        File GBZ_FILE
        File INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? HAPL_FILE
        File? DIST_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        String IN_OUTPUT_NAME_PREFIX = "haplotype_sampled_graph"
        Int IN_KMER_LENGTH = 29
        String IN_WORKING_DIRECTORY = "."
        String GIRAFFE_OPTIONS = ""
        String SAMPLE_NAME
        Int CORES = 16
    }

    String OUTPUT_NAME_PREFIX = select_first([IN_OUTPUT_NAME_PREFIX, "haplotype_sampled_graph"])
    Int KMER_LENGTH = select_first([IN_KMER_LENGTH, 29])
    String WORKING_DIRECTORY = select_first([IN_WORKING_DIRECTORY, "."])




    # Have to create haplotype information
    if (!defined(HAPL_FILE)) {
        # create the dist index file and r-index file to create the haplotype information file .hapl

        if (!defined(DIST_FILE)) {
            call map.createDistanceIndex {
                input:
                    in_gbz_file=GBZ_FILE
            }
        }

        File dist_index_file = select_first([DIST_FILE, createDistanceIndex.output_dist_index])

        if (!defined(R_INDEX_FILE)) {
            call map.createRIndex {
                input:
                    in_gbz_file=GBZ_FILE
            }
        }

        File r_index_file = select_first([R_INDEX_FILE, createRIndex.output_R_index])

        # create the haplotype information file
        call map.createHaplotypeIndex {
            input:
                in_gbz_file=GBZ_FILE,
                in_dist_index=dist_index_file,
                in_R_index=r_index_file
        }
    }

    File haplotype_index = select_first([HAPL_FILE, createHaplotypeIndex.output_hap_index])

    if (!defined(KFF_FILE)) {
        call utils.kmerCountingKMC {
            input:
                input_read_file_1=INPUT_READ_FILE_1,
                input_read_file_2=INPUT_READ_FILE_2,
                output_file_name=OUTPUT_NAME_PREFIX,
                kmer_length=KMER_LENGTH,
                working_directory=WORKING_DIRECTORY

        }

    }

    File kmer_information = select_first([KFF_FILE, kmerCountingKMC.kff_file])

    call map.samplingHaplotypes {
        input:
            in_gbz_file=GBZ_FILE,
            in_hap_index=haplotype_index,
            in_kmer_info=kmer_information,
            output_file_name=OUTPUT_NAME_PREFIX,
            working_directory=WORKING_DIRECTORY

    }

    call map.createDistanceIndex as giraffeDist {
                input:
                    in_gbz_file=samplingHaplotypes.output_graph
            }

    call map.createMinimizerIndex {
        input:
            in_gbz_file=samplingHaplotypes.output_graph,
            out_name=OUTPUT_NAME_PREFIX,
            in_dist_index=giraffeDist.output_dist_index

    }



#    call map.runVGGIRAFFENoIndex {
#        input:
#            fastq_file_1=INPUT_READ_FILE_1,
#            fastq_file_2=INPUT_READ_FILE_2,
#            in_giraffe_options=GIRAFFE_OPTIONS,
#            in_gbz_file=samplingHaplotypes.output_graph,
#            in_sample_name=SAMPLE_NAME,
#            nb_cores=CORES,
#            mem_gb=120,
#            out_prefix=OUTPUT_NAME_PREFIX
#    }

    output {
        File? sampled_graph = samplingHaplotypes.output_graph
        File? sampled_min = giraffeDist.output_dist_index
        File? sampled_dist = createMinimizerIndex.output_minimizer
    }

}






