version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map
import "./deepvariant.wdl" as dv_wf
import "./haplotype_sampling.wdl" as hapl

workflow GiraffeDeepVariant {

    meta {
        description: "## Giraffe-DeepVariant workflow \n The full workflow to go from sequencing reads (FASTQs, CRAM) to small variant calls (VCF). Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). DeepVariant then calls small variants. More information at [https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-workflow)."
    }

    parameter_meta {
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz"
        INPUT_CRAM_FILE: "Input CRAM file"
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "Path to .dist index file"
        MIN_FILE: "Path to .min index file"
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignemnt, path to .zipcodes index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'true'."
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling bams' (one per contig) won't be outputed. Default is 'true'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is 1. If null, uses DeepVariant default for the model type."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default or fast)"
        GIRAFFE_OPTIONS: "(OPTIONAL) Extra command line options for Giraffe mapper"
        TRUTH_VCF: "Path to .vcf.gz to compare against"
        TRUTH_VCF_INDEX: "Path to Tabix index for TRUTH_VCF"
        EVALUATION_REGIONS_BED: "BED to restrict comparison against TRUTH_VCF to"
        TARGET_REGION: "contig or region to restrict evaluation to"
        RUN_STANDALONE_VCFEVAL: "whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? Default is 'true'. Should be 'false' for HiFi."
        DV_NORM_READS: "Should DV normalize reads itself? Default is 'false'. Should be 'true' for HiFi."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'false'"
        IN_DIPLOID:"Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'"
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
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
        File GBZ_FILE
        File DIST_FILE
        File MIN_FILE
        File? ZIPCODES_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_GAF = true
        Boolean OUTPUT_SINGLE_BAM = false
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
        Int? MIN_MAPQ = 1
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        File? TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File? EVALUATION_REGIONS_BED
        String? TARGET_REGION
        Boolean RUN_STANDALONE_VCFEVAL = true
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Boolean DV_KEEP_LEGACY_AC = true
        Boolean DV_NORM_READS = false
        String OTHER_MAKEEXAMPLES_ARG = ""
        Int SPLIT_READ_CORES = 8
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Boolean HAPLOTYPE_SAMPLING = false
        Boolean IN_DIPLOID = true
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        String VG_DOCKER = "quay.io/vgteam/vg:v1.51.0"
        String? VG_GIRAFFE_DOCKER
        String? VG_SURJECT_DOCKER
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

    if (HAPLOTYPE_SAMPLING) {
        call hapl.HaplotypeSampling {
        input:
            IN_GBZ_FILE=GBZ_FILE,
            INPUT_READ_FILE_FIRST=read_1_file,
            INPUT_READ_FILE_SECOND=INPUT_READ_FILE_2,
            HAPL_FILE=IN_HAPL_FILE,
            IN_DIST_FILE=DIST_FILE,
            R_INDEX_FILE=IN_R_INDEX_FILE,
            KFF_FILE=IN_KFF_FILE,
            CORES=MAP_CORES,
            HAPLOTYPE_NUMBER=IN_HAPLOTYPE_NUMBER,
            DIPLOID=IN_DIPLOID,
        }

    }



    File file_gbz = select_first([HaplotypeSampling.sampled_graph, GBZ_FILE])
    File file_min = select_first([HaplotypeSampling.sampled_min, MIN_FILE])
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
            in_gbz_file=file_gbz,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=MAP_MEM,
            vg_docker=VG_DOCKER
        }
    }
    if (defined(REFERENCE_FILE)) {
        call utils.uncompressReferenceIfNeeded {
            input:
            # We know REFERENCE_FILE is defined but the WDL type system doesn't.
            in_reference_file=select_first([REFERENCE_FILE]),
        }
    }
    File reference_file = select_first([uncompressReferenceIfNeeded.reference_file, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE) || !defined(REFERENCE_DICT_FILE)) {
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
                in_preset=GIRAFFE_PRESET,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=file_gbz,
                in_dist_file=file_dist,
                in_min_file=file_min,
                in_zipcodes_file=ZIPCODES_FILE,
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
                in_min_file=file_min,
                in_zipcodes_file=ZIPCODES_FILE,
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

    call utils.mergeAlignmentBAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=sortBAM.sorted_bam
    }

    # Run the DeepVariant calling workflow
    call dv_wf.DeepVariant {
        input:
        MERGED_BAM_FILE=mergeAlignmentBAMChunks.merged_bam_file,
        MERGED_BAM_FILE_INDEX=mergeAlignmentBAMChunks.merged_bam_file_index,
        SAMPLE_NAME=SAMPLE_NAME,
        PATH_LIST_FILE=pipeline_path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=reference_file,
        REFERENCE_INDEX_FILE=reference_index_file,
        REFERENCE_DICT_FILE=reference_dict_file,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
        EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
        TARGET_REGION=TARGET_REGION,
        RUN_STANDALONE_VCFEVAL=RUN_STANDALONE_VCFEVAL,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_MODEL_META=DV_MODEL_META,
        DV_MODEL_INDEX=DV_MODEL_INDEX,
        DV_MODEL_DATA=DV_MODEL_DATA,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        SPLIT_READ_CORES=SPLIT_READ_CORES,
        REALIGN_MEM=if MAP_MEM < 40 then MAP_MEM else 40,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM
    }
    
    

    if (OUTPUT_GAF){
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=gaf_chunks,
            vg_docker=VG_DOCKER
        }
    }

    if (OUTPUT_SINGLE_BAM){
        call utils.mergeAlignmentBAMChunks as mergeBAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=select_all(flatten([DeepVariant.output_calling_bams, [DeepVariant.output_unmapped_bam]]))
        }
    }

    if (!OUTPUT_SINGLE_BAM){
        Array[File] output_calling_bam_files = DeepVariant.output_calling_bams
        Array[File] output_calling_bam_index_files = DeepVariant.output_calling_bam_indexes
    }

    output {
        File? output_vcfeval_evaluation_archive = DeepVariant.output_vcfeval_evaluation_archive
        File? output_happy_evaluation_archive = DeepVariant.output_happy_evaluation_archive
        File output_vcf = DeepVariant.output_vcf
        File output_vcf_index = DeepVariant.output_vcf_index
        File output_gvcf = DeepVariant.output_gvcf
        File output_gvcf_index = DeepVariant.output_gvcf_index
        File? output_gaf = mergeGAF.output_merged_gaf
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        Array[File]? output_calling_bams = output_calling_bam_files
        Array[File]? output_calling_bam_indexes = output_calling_bam_index_files
    }   
}
