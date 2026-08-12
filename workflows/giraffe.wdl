version 1.0

import "../tasks/validation.wdl" as validation
import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/vg_indexing.wdl" as index
import "./haplotype_sampling.wdl" as hapl
import "./internal/prepare_reference.wdl" as reference_wf
import "./internal/split_reads.wdl" as reads_wf
import "./internal/surject.wdl" as surject_wf

workflow Giraffe {
    meta {
        description: "## Giraffe workflow \n Core VG Giraffe mapping, usable for DeepVariant. Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). More information at [https://github.com/vgteam/vg_wdl/tree/master#giraffe-workflow](https://github.com/vgteam/vg_wdl/tree/master#giraffe-workflow)."
    }
    parameter_meta {
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz or fastq"
        INPUT_CRAM_FILE: "Input CRAM file to realign"
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        INPUT_BAM_FILE: "Input BAM file to realign"
        READ_CHUNKS_1: "(OPTIONAL) Input reads to map (either all reads or read 1), already split. When used, INPUT_READ_FILE_1 is still used for haplotype sampling."
        READ_CHUNKS_2: "(OPTIONAL) Input reads to map (read 2), in the same order as READ_CHUNKS_1. Only used with READ_CHUNKS_1, when the reads are paired and not interleaved."
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "(OPTIONAL) Path to .dist index file for the graph that will actually be mapped against (the haplotype-sampled graph, if HAPLOTYPE_SAMPLING is used). Generated from that graph if not given."
        MIN_FILE: "(OPTIONAL) Path to .min index file for the graph that will actually be mapped against (the haplotype-sampled graph, if HAPLOTYPE_SAMPLING is used). Generated from that graph if not given."
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? Default is 'true'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs be saved? Default is 'false'."
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'false'."
        OUTPUT_GAF_CHUNKS: "Should the unmerged GAF chunks be saved? Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 million."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. If using REFERENCE_PREFIX, contig names in here should not have the prefix. This is used in BAM processing and not for choosing contigs for the surjection, which uses PATH_LIST_FILE."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'."  
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"  
        GIRAFFE_OPTIONS: "(OPTIONAL) extra command line options for Giraffe mapper"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when splitting the reads into chunks. Default is 50."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'true'"
        DIPLOID:"Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'"
        SET_REFERENCE:"(OPTIONAL) Name of the single reference to keep for haplotype sampling."
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        R_INDEX_FILE: "(OPTIONAL) Path to .ri file used in haplotype sampling"
        KFF_FILE: "(OPTIONAL) Path to .kff file used in haplotype sampling"
        HAPLOTYPE_NUMBER: "Number of generated synthetic haplotypes used in haplotype sampling. (Default: 32)"
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)"
        OUTPUT_HAPL: "Whether or not to output the haplotype index (.hapl) created before haplotype sampling. This is useful if the same pangenome will be used for DeepVariant calling after mapping. Default is 'false'."
        VG_DOCKER: "Container image to use when running vg"
        VG_GIRAFFE_DOCKER: "Alternate container image to use when running vg giraffe mapping"
        VG_SURJECT_DOCKER: "Alternate container image to use when running vg surject"
    }
    input {
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File? INPUT_BAM_FILE
        Array[File]? READ_CHUNKS_1
        Array[File]? READ_CHUNKS_2
        File GBZ_FILE
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean OUTPUT_CALLING_BAMS = false
        Boolean OUTPUT_GAF = false
        Boolean OUTPUT_GAF_CHUNKS = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
        File? PATH_LIST_FILE
        Array[String]+? CONTIGS
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Boolean PRUNE_LOW_COMPLEXITY = true
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = if MAP_MEM < 40 then MAP_MEM else 40
        Boolean HAPLOTYPE_SAMPLING = true
        Boolean DIPLOID = true
        String? SET_REFERENCE
        File? HAPL_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        Int HAPLOTYPE_NUMBER = 32
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int HAPLOTYPE_INDEXING_MEM = 120
        Boolean OUTPUT_HAPL = false

        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? VG_GIRAFFE_DOCKER
        String? VG_SURJECT_DOCKER
    }

    # Set up enough path list stuff to do validation as soon as possible.
    # We don't just check the final path list because that could come from
    # haplotype sampling and we want to check any provided path list before going
    # into haplotype sampling.
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
        if (REFERENCE_PREFIX != "") {
            call validation.checkPathList as checkWrittenPathList {
                input:
                    in_path_list_file=written_path_names_file,
                    in_reference_prefix=REFERENCE_PREFIX
            }
        }
    }
    if (defined(PATH_LIST_FILE) && REFERENCE_PREFIX != "") {
        call validation.checkPathList as checkProvidedPathList {
            input:
                in_path_list_file=select_first([PATH_LIST_FILE]),
                in_reference_prefix=REFERENCE_PREFIX
        }
    }

    # Validate the dict file if fed in.
    if (defined(REFERENCE_DICT_FILE) && REFERENCE_PREFIX != "") {
        call validation.checkDict as checkProvidedDict {
            input:
                in_dict_file=select_first([REFERENCE_DICT_FILE]),
                in_reference_prefix=REFERENCE_PREFIX
        }
    }


    if (!defined(READ_CHUNKS_1)) {
        # Get the reads as FASTQ, whatever container they came in, and split them
        # up for parallel mapping.
        call reads_wf.SplitReads {
            input:
            INPUT_READ_FILE_1=INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=INPUT_READ_FILE_2,
            INPUT_CRAM_FILE=INPUT_CRAM_FILE,
            CRAM_REF=CRAM_REF,
            CRAM_REF_INDEX=CRAM_REF_INDEX,
            INPUT_BAM_FILE=INPUT_BAM_FILE,
            # Haplotype sampling needs to see all of a sample's reads at once,
            # so keep the whole FASTQs when we are going to sample.
            OUTPUT_WHOLE_READS=HAPLOTYPE_SAMPLING,
            PAIRED_READS=PAIRED_READS,
            INTERLEAVED_READS=INTERLEAVED_READS,
            READS_PER_CHUNK=READS_PER_CHUNK,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_MEM=SPLIT_READ_MEM
        }
    }

    Array[File] read_chunks_1 = select_first([READ_CHUNKS_1, SplitReads.read_chunks_1])
    # Null when the reads are single-ended or interleaved, in which case the
    # first set of chunks holds everything.
    Array[File]? read_chunks_2 = if defined(READ_CHUNKS_1) then READ_CHUNKS_2 else SplitReads.read_chunks_2
    # Whole reads for haplotype sampling, if available
    File? whole_read_1_file = if defined(SplitReads.read_1_file) then SplitReads.read_1_file else INPUT_READ_FILE_1
    File? whole_read_2_file = if defined(SplitReads.read_2_file) then SplitReads.read_2_file else INPUT_READ_FILE_2 

    if (HAPLOTYPE_SAMPLING) {
        call hapl.HaplotypeSampling {
        input:
            GBZ_FILE=GBZ_FILE,
            INPUT_READ_FILE_FIRST=select_first([whole_read_1_file]),
            INPUT_READ_FILE_SECOND=whole_read_2_file,
            HAPL_FILE=HAPL_FILE,
            DIST_FILE=DIST_FILE,
            R_INDEX_FILE=R_INDEX_FILE,
            KFF_FILE=KFF_FILE,
            HAPLOTYPE_NUMBER=HAPLOTYPE_NUMBER,
            DIPLOID=DIPLOID,
            SET_REFERENCE=SET_REFERENCE,
            CORES=MAP_CORES,
            KMER_COUNTING_MEM=KMER_COUNTING_MEM,
            HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
            OUTPUT_HAPL=OUTPUT_HAPL,
            VG_DOCKER=VG_DOCKER
        }

    }

    File file_gbz = select_first([HaplotypeSampling.sampled_graph, GBZ_FILE])

    # DIST_FILE and MIN_FILE are indexes of the graph that will actually be
    # mapped against. If haplotype sampling ran, that's a different graph than
    # GBZ_FILE, so any DIST_FILE/MIN_FILE given by the caller (which can only
    # be indexes of GBZ_FILE) can't be reused and we have to rebuild them.
    if (HAPLOTYPE_SAMPLING || !defined(DIST_FILE)) {
        call index.createDistanceIndex as createGiraffeDist {
            input:
                in_gbz_file=file_gbz,
                nb_cores=MAP_CORES,
                in_extract_mem=HAPLOTYPE_INDEXING_MEM,
                vg_docker=VG_DOCKER
        }
    }
    File file_dist = select_first([createGiraffeDist.output_dist_index, DIST_FILE])

    if (HAPLOTYPE_SAMPLING || !defined(MIN_FILE)) {
        call index.createMinimizerIndex as createGiraffeMin {
            input:
                in_gbz_file=file_gbz,
                in_dist_index=file_dist,
                in_minimizer_k = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 29 else 31,
                in_minimizer_w = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 11 else 50,
                in_minimizer_weighted = INDEX_MINIMIZER_WEIGHTED,
                out_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                in_extract_mem=INDEX_MINIMIZER_MEM,
                vg_docker=VG_DOCKER
        }
    }
    File file_min = select_first([createGiraffeMin.output_minimizer, MIN_FILE])
    # The zipcode file is optional but we still have a priority list of places to get it from.
    # But we can't select_first since they all might be null.
    Array[File] possible_zipcode_files = select_all([createGiraffeMin.output_zipcodes, ZIPCODES_FILE])
    # We can't actually use None in WDL 1.0 so we need to use a nonexistent null file.
    if (false) {
        Array[File] no_files = []
        #@ except: UnnecessaryFunctionCall
        File NULL_FILE = select_first(no_files)
    }
    File? file_zipcodes = if length(possible_zipcode_files) > 0 then possible_zipcode_files[0] else NULL_FILE


    # Which path names to work on, and what reference to surject against? These
    # come from the graph we are actually mapping to, which is the sampled one
    # if we sampled.
    call reference_wf.PrepareReference {
        input:
        GBZ_FILE=file_gbz,
        CONTIGS=CONTIGS,
        PATH_LIST_FILE=PATH_LIST_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_DICT_FILE=REFERENCE_DICT_FILE,
        EXTRACT_MEM=MAP_MEM,
        VG_DOCKER=VG_DOCKER
    }
    File pipeline_path_list_file = PrepareReference.path_list_file
    File reference_file = PrepareReference.reference_file
    File reference_index_file = PrepareReference.reference_index_file
    File reference_dict_file = PrepareReference.reference_dict_file

    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    if (PAIRED_READS && !INTERLEAVED_READS) {
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(read_chunks_1, select_first([read_chunks_2]))
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call map.runVGGIRAFFE as runVGGIRAFFE2file {
                input:
                fastq_file_1=read_pair_chunk_files.left,
                fastq_file_2=read_pair_chunk_files.right,
                in_preset=GIRAFFE_PRESET,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=file_gbz,
                in_dist_file=file_dist,
                in_zipcodes_file=file_zipcodes,
                in_min_file=file_min,
                # We always need to pass a full dict file here, with lengths,
                # because if we pass just path lists and the paths are not
                # completely contained in the graph (like if we're working on
                # GRCh38 paths in a CHM13-based graph), giraffe won't be able
                # to get the path lengths and will crash.
                # TODO: Somehow this problem is supposed to go away if we pull
                # any GRCh38. prefix off the path names by setting
                # REFERENCE_PREFIX and making sure the prefix isn't in the
                # truth set.
                # See <https://github.com/adamnovak/giraffe-dv-wdl/pull/2#issuecomment-955096920>
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM,
                vg_docker=select_first([VG_GIRAFFE_DOCKER, VG_DOCKER])
            }
        }
    }
    # TODO: invent else
    if (!(PAIRED_READS && !INTERLEAVED_READS)) {
        scatter (read_pair_chunk_file in read_chunks_1) {
            call map.runVGGIRAFFE as runVGGIRAFFE1file {
                input:
                fastq_file_1=read_pair_chunk_file,
                in_interleaved=INTERLEAVED_READS,
                in_preset=GIRAFFE_PRESET,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=file_gbz,
                in_dist_file=file_dist,
                in_zipcodes_file=file_zipcodes,
                in_min_file=file_min,
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM,
                vg_docker=select_first([VG_GIRAFFE_DOCKER, VG_DOCKER])
            }
        }
    }

    Array[File] gaf_chunks = select_first([runVGGIRAFFE2file.chunk_gaf_file, runVGGIRAFFE1file.chunk_gaf_file])

    if (OUTPUT_GAF) {
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=gaf_chunks
        }
    }

    if (OUTPUT_GAF_CHUNKS) {
        # To avoid always saving all the GAF chunks as part of our output,
        # we make sure this is null unless the caller actually asked for
        # GAF chunks.
        Array[File] wanted_gaf_chunks = gaf_chunks
    }

    if (OUTPUT_SINGLE_BAM || OUTPUT_CALLING_BAMS) {
        # We are outputting BAM so surjection is needed
        call surject_wf.Surject {
            input:
            GAF_CHUNKS=gaf_chunks,
            GBZ_FILE=file_gbz,
            PATH_LIST_FILE=pipeline_path_list_file,
            REFERENCE_DICT_FILE=reference_dict_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            SAMPLE_NAME=SAMPLE_NAME,
            PAIRED_READS=PAIRED_READS,
            PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
            MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
            SURJECT_CORES=MAP_CORES,
            SURJECT_MEM=MAP_MEM,
            BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
            VG_DOCKER=select_first([VG_SURJECT_DOCKER, VG_DOCKER])
        }

        if (OUTPUT_CALLING_BAMS || LEFTALIGN_BAM || REALIGN_INDELS) {
            # We will need to split up the BAM by contig to do processing on it.
            # TODO: Unify BAM processing with deepvariant.wdl

            # Split merged alignment by contigs list
            call utils.splitBAMbyPath {
                input:
                in_sample_name=SAMPLE_NAME,
                in_merged_bam_file=Surject.bam_file,
                in_merged_bam_file_index=Surject.bam_index_file,
                in_path_list_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
                thread_count=SPLIT_READ_CORES,
                mem_gb=SPLIT_READ_MEM
            }

            ##
            ## Prepare each BAM
            ##
            scatter (bam_and_index_for_path in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
                ## Evantually shift and realign reads
                if (LEFTALIGN_BAM){
                    # Just left-shift each read individually
                    call utils.leftShiftBAMFile {
                        input:
                        in_bam_file=bam_and_index_for_path.left,
                        in_reference_file=reference_file,
                        in_reference_index_file=reference_index_file,
                        mem_gb=BAM_PREPROCESS_MEM
                    }
                }
                if (REALIGN_INDELS) {
                    File forrealign_bam = select_first([leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
                    File forrealign_index = select_first([leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
                    # Do indel realignment
                    call utils.prepareRealignTargets {
                        input:
                        in_bam_file=forrealign_bam,
                        in_bam_index_file=forrealign_index,
                        in_reference_file=reference_file,
                        in_reference_index_file=reference_index_file,
                        in_reference_dict_file=reference_dict_file,
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                        thread_count=MAP_CORES,
                        mem_gb=BAM_PREPROCESS_MEM
                    }
                    call utils.runAbraRealigner {
                        input:
                            in_bam_file=forrealign_bam,
                            in_bam_index_file=forrealign_index,
                            in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                            in_reference_file=reference_file,
                            in_reference_index_file=reference_index_file,
                            threadCount=MAP_CORES,
                            memoryGb=REALIGN_MEM
                    }
                }
                File processed_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
                File processed_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
            }
        
            if (OUTPUT_SINGLE_BAM && (LEFTALIGN_BAM || REALIGN_INDELS)){
                # We're outputting one big BAM and we've actually made changes to the chunks. Put them back together.
                call utils.mergeAlignmentBAMChunks as mergeBAM {
                    input:
                    in_sample_name=SAMPLE_NAME,
                    in_alignment_bam_chunk_files=flatten([processed_bam, [splitBAMbyPath.bam_unmapped_file]]),
                    in_cores=SPLIT_READ_CORES,
                    mem_gb=SPLIT_READ_MEM
                }
            }
            
            if (OUTPUT_CALLING_BAMS){
                Array[File] calling_bams = processed_bam
                Array[File] calling_bam_indexes = processed_bam_index
            }
        }

        if (OUTPUT_SINGLE_BAM) {
            # Find the single BAM and index that we want to output.
            # We want the one after postprocessing if we did any, and the plain merged sorted BAM otherwise.
            File single_bam = select_first([mergeBAM.merged_bam_file, Surject.bam_file])
            File single_bam_index = select_first([mergeBAM.merged_bam_file_index, Surject.bam_index_file])
        }
    }

    output {
        File? output_bam = single_bam
        File? output_bam_index = single_bam_index
        File? output_gaf = mergeGAF.output_merged_gaf
        Array[File]? output_gaf_chunks = wanted_gaf_chunks
        Array[File]? output_calling_bams = calling_bams
        Array[File]? output_calling_bam_indexes = calling_bam_indexes
        File? output_hapl_file = HaplotypeSampling.output_hapl_file
    }
}


