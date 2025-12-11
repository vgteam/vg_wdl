version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/vg_map_hts.wdl" as map
import "./deepvariant.wdl" as dv_wf
import "./giraffe.wdl" as giraffe_wf


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
        ZIPCODES_FILE: "(OPTIONAL) For chaining-based alignment, path to .zipcodes index file"
        HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        SAMPLE_NAME: "The sample name"
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'true'."
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling bams' (one per contig) won't be outputed by default. Default is 'false'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM."
        OUTPUT_UNMAPPED_BAM: "Should an unmapped reads BAM be saved? Default is false."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        HAPLOID_CONTIGS: "(OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5."
        PAR_REGIONS_BED_FILE: "(OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realingment. Default is 'true'."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"
        GIRAFFE_OPTIONS: "(OPTIONAL) Extra command line options for Giraffe mapper"
        TRUTH_VCF: "Path to .vcf.gz to compare against"
        TRUTH_VCF_INDEX: "Path to Tabix index for TRUTH_VCF"
        EVALUATION_REGIONS_BED: "BED to evaluate against TRUTH_VCF on, where false positives will be counted"
        RESTRICT_REGIONS_BED: "BED to restrict comparison against TRUTH_VCF to"
        TARGET_REGION: "contig or region to restrict evaluation to"
        RUN_STANDALONE_VCFEVAL: "whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_MODEL_FILES: "Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format"
        DV_MODEL_VARIABLES_FILES: "Array of files that need to go in a 'variables' subdirectory for a DV model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+."
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when splitting the reads into chunks. Default is 50."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'true'."
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        EVAL_MEM: "Memory, in GB, to use when evaluating variant calls. Default is 60."
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
        File? DIST_FILE
        File? MIN_FILE
        File? ZIPCODES_FILE
        File? HAPL_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_GAF = true
        Boolean OUTPUT_SINGLE_BAM = false
        Boolean OUTPUT_CALLING_BAMS = !OUTPUT_SINGLE_BAM
        Boolean OUTPUT_UNMAPPED_BAM = false
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = false
        Int READS_PER_CHUNK = 20000000
        Array[String]+? CONTIGS
        File? PATH_LIST_FILE
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Array[String]? HAPLOID_CONTIGS
        File? PAR_REGIONS_BED_FILE
        Boolean PRUNE_LOW_COMPLEXITY = true
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int? MIN_MAPQ
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_PRESET = "default"
        String GIRAFFE_OPTIONS = ""
        File? TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File? EVALUATION_REGIONS_BED
        File? RESTRICT_REGIONS_BED
        String? TARGET_REGION
        Boolean RUN_STANDALONE_VCFEVAL = true
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Array[File]? DV_MODEL_FILES
        Array[File]? DV_MODEL_VARIABLES_FILES
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        String OTHER_MAKEEXAMPLES_ARG = ""
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Boolean HAPLOTYPE_SAMPLING = true
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = if MAP_MEM < 40 then MAP_MEM else 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int EVAL_MEM = 60
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? VG_GIRAFFE_DOCKER
        String? VG_SURJECT_DOCKER

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
            in_gbz_file=GBZ_FILE,
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

    # Run the giraffe mapping workflow.
    # We don't do postprocessing in the Giraffe workflow, just the DV workflow.
    # Otherwise we'd split to contig BAMs, process, re-merge, and re-split.
    call giraffe_wf.Giraffe {
        input:
        INPUT_READ_FILE_1=INPUT_READ_FILE_1,
        INPUT_READ_FILE_2=INPUT_READ_FILE_2,
        INPUT_CRAM_FILE=INPUT_CRAM_FILE,
        CRAM_REF=CRAM_REF,
        CRAM_REF_INDEX=CRAM_REF_INDEX,
        GBZ_FILE=GBZ_FILE,
        DIST_FILE=DIST_FILE,
        MIN_FILE=MIN_FILE,
        ZIPCODES_FILE=ZIPCODES_FILE,
        HAPL_FILE=HAPL_FILE,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=true,
        OUTPUT_CALLING_BAMS=false,
        OUTPUT_GAF=OUTPUT_GAF,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        READS_PER_CHUNK=READS_PER_CHUNK,
        PATH_LIST_FILE=pipeline_path_list_file,
        CONTIGS=CONTIGS,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=reference_file,
        REFERENCE_INDEX_FILE=reference_index_file,
        REFERENCE_DICT_FILE=reference_dict_file,
        PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
        LEFTALIGN_BAM=false,
        REALIGN_INDELS=false,
        MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
        GIRAFFE_PRESET=GIRAFFE_PRESET,
        GIRAFFE_OPTIONS=GIRAFFE_OPTIONS,
        SPLIT_READ_CORES=SPLIT_READ_CORES,
        SPLIT_READ_MEM=SPLIT_READ_MEM,
        MAP_CORES=MAP_CORES,
        MAP_MEM=MAP_MEM,
        HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        VG_DOCKER=VG_DOCKER,
        VG_GIRAFFE_DOCKER=VG_GIRAFFE_DOCKER,
        VG_SURJECT_DOCKER=VG_SURJECT_DOCKER
    }

    # Run the DeepVariant calling workflow
    call dv_wf.DeepVariant {
        input:
        MERGED_BAM_FILE=select_first([Giraffe.output_bam]),
        MERGED_BAM_FILE_INDEX=select_first([Giraffe.output_bam_index]),
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=OUTPUT_SINGLE_BAM,
        OUTPUT_CALLING_BAMS=OUTPUT_CALLING_BAMS,
        OUTPUT_UNMAPPED_BAM=OUTPUT_UNMAPPED_BAM,
        PATH_LIST_FILE=pipeline_path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=reference_file,
        REFERENCE_INDEX_FILE=reference_index_file,
        REFERENCE_DICT_FILE=reference_dict_file,
        HAPLOID_CONTIGS=HAPLOID_CONTIGS,
        PAR_REGIONS_BED_FILE=PAR_REGIONS_BED_FILE,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
        EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
        RESTRICT_REGIONS_BED=RESTRICT_REGIONS_BED,
        TARGET_REGION=TARGET_REGION,
        RUN_STANDALONE_VCFEVAL=RUN_STANDALONE_VCFEVAL,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_MODEL_META=DV_MODEL_META,
        DV_MODEL_INDEX=DV_MODEL_INDEX,
        DV_MODEL_DATA=DV_MODEL_DATA,
        DV_MODEL_FILES=DV_MODEL_FILES,
        DV_MODEL_VARIABLES_FILES=DV_MODEL_VARIABLES_FILES,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        REALIGN_MEM=REALIGN_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        MAKE_EXAMPLES_CORES=MAKE_EXAMPLES_CORES,
        MAKE_EXAMPLES_MEM=MAKE_EXAMPLES_MEM,
        EVAL_MEM=EVAL_MEM
    }
    
    output {
        File? output_vcfeval_evaluation_archive = DeepVariant.output_vcfeval_evaluation_archive
        File? output_happy_evaluation_archive = DeepVariant.output_happy_evaluation_archive
        File output_vcf = DeepVariant.output_vcf
        File output_vcf_index = DeepVariant.output_vcf_index
        File output_gvcf = DeepVariant.output_gvcf
        File output_gvcf_index = DeepVariant.output_gvcf_index
        File? output_gaf = Giraffe.output_gaf
        File? output_bam = DeepVariant.output_bam
        File? output_bam_index = DeepVariant.output_bam_index
        Array[File]? output_calling_bams = DeepVariant.output_calling_bams
        Array[File]? output_calling_bam_indexes = DeepVariant.output_calling_bam_indexes
        File? output_unmapped_bam = DeepVariant.output_unmapped_bam
    }   
}
