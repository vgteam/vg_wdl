version 1.0

import "../tasks/validation.wdl" as validation
import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map
import "./haplotype_sampling.wdl" as hapl

workflow Giraffe {
    meta {
        description: "## Giraffe workflow \n Core VG Giraffe mapping, usable for DeepVariant. Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). More information at [https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-workflow)."
    }
    parameter_meta {
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz or fastq"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz or fastq"
        INPUT_CRAM_FILE: "Input CRAM file to realign"
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        INPUT_BAM_FILE: "Input BAM file to realign"
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "Path to .dist index file. Optional if using haplotype sampling."
        MIN_FILE: "Path to .min index file. Optional if using haplotype sampling."
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? Default is 'true'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs be saved? Default is 'false'."
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'false'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index. If using REFERENCE_PREFIX, contig names in here should have the prefix."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file. If using REFERENCE_PREFIX, contig names in here should not have the prefix."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. If using REFERENCE_PREFIX, contig names in here should not have the prefix. This is used in BAM processing and not for choosing contigs for the surjection, which uses PATH_LIST_FILE."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."  
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
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'true'"
        DIPLOID:"Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'"
        SET_REFERENCE:"(OPTIONAL) Name of the single reference to keep for haplotype sampling."
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        R_INDEX_FILE: "(OPTIONAL) Path to .ri file used in haplotype sampling"
        KFF_FILE: "(OPTIONAL) Path to .kff file used in haplotype sampling"
        HAPLOTYPE_NUMBER: "Number of generated synthetic haplotypes used in haplotype sampling. (Default: 4)"
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)" 

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
        File GBZ_FILE
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean OUTPUT_CALLING_BAMS = false
        Boolean OUTPUT_GAF = false
        Boolean PAIRED_READS = true
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
        Boolean HAPLOTYPE_SAMPLING = true
        Boolean DIPLOID = true
        String? SET_REFERENCE
        File? HAPL_FILE
        File? R_INDEX_FILE
        File? KFF_FILE
        Int HAPLOTYPE_NUMBER = 4
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        
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


    if(defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	    call utils.convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES,
            in_memory=SPLIT_READ_MEM
	    }
    }

    if(defined(INPUT_BAM_FILE)) {
	    call utils.convertBAMtoFASTQ {
            input:
            in_bam_file=INPUT_BAM_FILE,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES,
            in_memory=SPLIT_READ_MEM
	    }
    }

    File read_1_file = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file, convertBAMtoFASTQ.output_fastq_1_file])
    if(PAIRED_READS){
        # We also need the second read in the pair, if paired, for hap sampling.
        File read_2_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file, convertBAMtoFASTQ.output_fastq_2_file])
    }

    if (HAPLOTYPE_SAMPLING) {
        call hapl.HaplotypeSampling {
        input:
            GBZ_FILE=GBZ_FILE,
            INPUT_READ_FILE_FIRST=read_1_file,
            # If we're not doing paired reads the result here is probably null.
            INPUT_READ_FILE_SECOND=if PAIRED_READS then read_2_file else INPUT_READ_FILE_2, 
            HAPL_FILE=HAPL_FILE,
            DIST_FILE=DIST_FILE,
            R_INDEX_FILE=R_INDEX_FILE,
            KFF_FILE=KFF_FILE,
            HAPLOTYPE_NUMBER=HAPLOTYPE_NUMBER,
            DIPLOID=DIPLOID,
            SET_REFERENCE=SET_REFERENCE,
            INDEX_MINIMIZER_K = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 29 else 31,
            INDEX_MINIMIZER_W = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 11 else 50,
            INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
            CORES=MAP_CORES,
            INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
            VG_DOCKER=VG_DOCKER
        }

    }

    File file_gbz = select_first([HaplotypeSampling.sampled_graph, GBZ_FILE])
    File file_min = select_first([HaplotypeSampling.sampled_min, MIN_FILE])
    # The zipcode file is optional but we still have a priority list of places to get it from.
    # But we can't select_first since they all might be null.
    Array[File] possible_zipcode_files = select_all([HaplotypeSampling.sampled_zipcodes, ZIPCODES_FILE])
    # We can't actually use None in WDL 1.0 so we need to use a nonexistent null file.
    if (false) {
        Array[File] no_files = []
        File NULL_FILE = select_first(no_files)
    }
    File? file_zipcodes = if length(possible_zipcode_files) > 0 then possible_zipcode_files[0] else NULL_FILE
    File file_dist = select_first([HaplotypeSampling.sampled_dist, DIST_FILE])


    # Split input reads into chunks for parallelized mapping
    call utils.splitReads as firstReadPair {
        input:
            in_read_file=read_1_file,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
    }
    
    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from GBZ file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call map.extractSubsetPathNames {
                input:
                    in_gbz_file=file_gbz,
                    in_reference_prefix=REFERENCE_PREFIX,
                    in_extract_mem=MAP_MEM,
                    vg_docker=VG_DOCKER
            }
            
            if (REFERENCE_PREFIX != "") {
                call validation.checkPathList as checkExtractedPathList {
                    input:
                        in_path_list_file=extractSubsetPathNames.output_path_list_file,
                        in_reference_prefix=REFERENCE_PREFIX
                }
            }
        }
    } 
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNames.output_path_list_file, written_path_names_file])

    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph (except for removing the name prefix), we
    # generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call map.extractReference {
            input:
            in_gbz_file=file_gbz,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=MAP_MEM,
            vg_docker=VG_DOCKER
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference {
            input:
                in_reference_file=reference_file
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])
    
    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    if(PAIRED_READS){
        call utils.splitReads as secondReadPair {
            input:
            in_read_file=select_first([read_2_file]),
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
        }
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call map.runVGGIRAFFE as runVGGIRAFFEpe {
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
    if (!PAIRED_READS) {
        scatter (read_pair_chunk_file in firstReadPair.output_read_chunks) {
            call map.runVGGIRAFFE as runVGGIRAFFEse {
                input:
                fastq_file_1=read_pair_chunk_file,
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

    Array[File] gaf_chunks = select_first([runVGGIRAFFEpe.chunk_gaf_file, runVGGIRAFFEse.chunk_gaf_file])
    
    if (OUTPUT_SINGLE_BAM || OUTPUT_CALLING_BAMS) {
        # We are outputting BAM so surjection is needed

        scatter (gaf_file in gaf_chunks) {
            call gautils.surjectGAFtoBAM {
                input:
                in_gaf_file=gaf_file,
                in_gbz_file=file_gbz,
                in_path_list_file=pipeline_path_list_file,
                in_sample_name=SAMPLE_NAME,
                in_max_fragment_length=MAX_FRAGMENT_LENGTH,
                in_paired_reads=PAIRED_READS,
                in_prune_low_complexity=PRUNE_LOW_COMPLEXITY,
                mem_gb=MAP_MEM,
                vg_docker=select_first([VG_SURJECT_DOCKER, VG_DOCKER])
            }

            call utils.sortBAM {
                input:
                in_bam_file=surjectGAFtoBAM.output_bam_file,
                in_ref_dict=reference_dict_file,
                in_prefix_to_strip=REFERENCE_PREFIX
            }
        }
        
        # Merge up the unprocessed surjected alignments
        call utils.mergeAlignmentBAMChunks {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=sortBAM.sorted_bam
        }

        if (OUTPUT_CALLING_BAMS || LEFTALIGN_BAM || REALIGN_INDELS) {
            # We will need to split up the BAM by contig to do processing on it.
            # TODO: Unify BAM processing with deepvariant.wdl

            # Split merged alignment by contigs list
            call utils.splitBAMbyPath {
                input:
                in_sample_name=SAMPLE_NAME,
                in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
                in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
                in_path_list_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX
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
                        in_reference_index_file=reference_index_file
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
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES
                    }
                    call utils.runAbraRealigner {
                        input:
                            in_bam_file=forrealign_bam,
                            in_bam_index_file=forrealign_index,
                            in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                            in_reference_file=reference_file,
                            in_reference_index_file=reference_index_file,
                            # If the user has set a very low memory for mapping, don't use more for realignment
                            memoryGb=if MAP_MEM < 40 then MAP_MEM else 40
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
                    in_alignment_bam_chunk_files=select_all(flatten([processed_bam, [splitBAMbyPath.bam_unmapped_file]]))
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
            File single_bam = select_first([mergeBAM.merged_bam_file, mergeAlignmentBAMChunks.merged_bam_file])
            File single_bam_index = select_first([mergeBAM.merged_bam_file_index, mergeAlignmentBAMChunks.merged_bam_file_index])
        }
    }
    
    if (OUTPUT_GAF){
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=gaf_chunks
        }
    }
    
    output {
        File? output_bam = single_bam
        File? output_bam_index = single_bam_index
        File? output_gaf = mergeGAF.output_merged_gaf
        Array[File]? output_calling_bams = calling_bams
        Array[File]? output_calling_bam_indexes = calling_bam_indexes

    }   
}

