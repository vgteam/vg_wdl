version 1.0

### giraffe.wdl ###
## Author: Charles Markello, Jean Monlong, Adam Novak
## Description: Core VG Giraffe mapping, usable for DeepVariant.
## Reference: https://github.com/vgteam/vg/wiki

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map

workflow Giraffe {
    input {
        File? INPUT_READ_FILE_1                         # Input sample 1st read pair fastq.gz
        File? INPUT_READ_FILE_2                         # Input sample 2nd read pair fastq.gz
        File? INPUT_CRAM_FILE                           # Input CRAM file
        File? CRAM_REF                                  # Genome fasta file associated with the CRAM file
        File? CRAM_REF_INDEX                            # Index of the fasta file associated with the CRAM file
        String SAMPLE_NAME                              # The sample name
        Boolean PAIRED_READS = true
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        String GIRAFFE_OPTIONS = ""                     # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths.
        String REFERENCE_PREFIX = ""                    # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File GBZ_FILE                                   # Path to .gbz index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
        Boolean LEFTALIGN_BAM = true                    # Whether or not to left-align reads in the BAM before DV
        Boolean REALIGN_INDELS = true                   # Whether or not to realign reads near indels before DV
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Boolean OUTPUT_SINGLE_BAM = true                # Should a single merged BAM file be saved?
        Boolean OUTPUT_CALLING_BAMS = false             # Should individual contig BAMs be saved?
        Boolean OUTPUT_GAF = false                      # Should a GAF file with the aligned reads be saved?
        Int SPLIT_READ_CORES = 8
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        File? REFERENCE_FILE                            # (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.
        File? REFERENCE_INDEX_FILE                      # (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
        File? REFERENCE_DICT_FILE                       # (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. 
    }

    if(defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	    call utils.convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES
	    }
    }

    File read_1_file = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file])
    
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
                    in_gbz_file=GBZ_FILE,
                    in_extract_mem=MAP_MEM
            }
        }
    } 
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNames.output_path_list_file, written_path_names_file])
    
    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call map.extractReference {
            input:
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=MAP_MEM
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
        File read_2_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file])
        call utils.splitReads as secondReadPair {
            input:
            in_read_file=read_2_file,
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
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_min_file=MIN_FILE,
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
                mem_gb=MAP_MEM
            }
        }
    }
    if (!PAIRED_READS) {
        scatter (read_pair_chunk_file in firstReadPair.output_read_chunks) {
            call map.runVGGIRAFFE as runVGGIRAFFEse {
                input:
                fastq_file_1=read_pair_chunk_file,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_min_file=MIN_FILE,
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
                mem_gb=MAP_MEM
            }
        }
    }

    Array[File] gaf_chunks = select_first([runVGGIRAFFEpe.chunk_gaf_file, runVGGIRAFFEse.chunk_gaf_file])
    scatter (gaf_file in gaf_chunks) {
        call gautils.surjectGAFtoSortedBAM {
            input:
            in_gaf_file=gaf_file,
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_paired_reads=PAIRED_READS
        }
        if (REFERENCE_PREFIX != "") {
            # use samtools to replace the header contigs with those from our dict.
            # this allows the header to contain contigs that are not in the graph,
            # which is more general and lets CHM13-based graphs be used to call on GRCh38
            # also, strip out contig prefixes in the BAM body
            call map.fixBAMContigNaming {
                input:
                in_bam_file=surjectGAFtoSortedBAM.output_bam_file,
                in_ref_dict=reference_dict_file,
                in_prefix_to_strip=REFERENCE_PREFIX
            }
        }
        File properly_named_bam_file = select_first([fixBAMContigNaming.fixed_bam_file, surjectGAFtoSortedBAM.output_bam_file]) 
    }

    call utils.mergeAlignmentBAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=properly_named_bam_file
    }
    
    if (REFERENCE_PREFIX != "") {
        # strip all the GRCh38's off our path list file.  we need them for surject as they are in the path
        # but fixBAMContigNaming above stripped them, so we don't need them downstream
        call map.fixPathNames {
            input:
                in_path_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
        }
    }
    File properly_named_path_list_file = select_first([fixPathNames.fixed_path_list_file, pipeline_path_list_file])
             
    # Split merged alignment by contigs list
    call utils.splitBAMbyPath { 
        input:
    in_sample_name=SAMPLE_NAME,
    in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
    in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
    in_path_list_file=properly_named_path_list_file
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
                    in_reference_index_file=reference_index_file
            }
        }
        File processed_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
        File processed_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
    }
    
    if (OUTPUT_SINGLE_BAM){
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
    
    if (OUTPUT_GAF){
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=gaf_chunks
        }
    }
    
    output {
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        File? output_gaf = mergeGAF.output_merged_gaf
        Array[File]? output_calling_bams = calling_bams
        Array[File]? output_calling_bam_indexes = calling_bam_indexes
    }   
}

