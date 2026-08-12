version 1.0

import "./giraffe.wdl" as giraffe_wf
import "./giraffe_and_deepvariant_fromGAF.wdl" as gaf_wf
import "./internal/prepare_reference.wdl" as reference_wf


workflow GiraffeDeepVariant {

    meta {
        description: "## Giraffe-DeepVariant workflow \n The full workflow to go from sequencing reads (FASTQs, CRAM) to small variant calls (VCF). Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). DeepVariant then calls small variants. More information at [https://github.com/vgteam/vg_wdl/tree/master#giraffe-deepvariant-workflow](https://github.com/vgteam/vg_wdl/tree/master#giraffe-deepvariant-workflow)."
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
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted by default. Default is 'false'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM."
        OUTPUT_UNMAPPED_BAM: "Should an unmapped reads BAM be saved? Default is false."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 million."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        HAPLOID_CONTIGS: "(OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5."
        PAR_REGIONS_BED_FILE: "(OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type."
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)"
        GIRAFFE_OPTIONS: "(OPTIONAL) Extra command line options for Giraffe mapper"
        TRUTH_VCF: "Path to .vcf.gz to compare against"
        TRUTH_VCF_INDEX: "Path to Tabix index for TRUTH_VCF"
        EVALUATION_REGIONS_BED: "BED to evaluate against TRUTH_VCF on, where false positives will be counted. Required when EVALUATE_WITH_AARDVARK is set."
        EVALUATE_WITH_AARDVARK: "Should the calls be compared to TRUTH_VCF with Aardvark instead of hap.py? Default is 'false'."
        STRATIFICATION_ARCHIVE: "(OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the results down by. Only used when EVALUATE_WITH_AARDVARK is set."
        RESTRICT_REGIONS_BED: "BED to restrict comparison against TRUTH_VCF to"
        TARGET_REGION: "contig or region to restrict evaluation to"
        RUN_STANDALONE_VCFEVAL: "whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_MODEL_FILES: "Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format"
        DV_MODEL_VARIABLES_FILES: "Array of files that need to go in a 'variables' subdirectory for a DV model"
        DV_PANGENOME_GBZ: "(OPTIONAL) Path to a pangenome graph in GBZ format for pangenome-aware DV."
        DV_PANGENOME_IMAGE_HEIGHT: "(OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware DV. It will be used only if DV_PANGENOME_GBZ is set. If DV_PANGENOME_HAPLOTYPE_SAMPLING is done by this workflow and this is not set, it defaults to DV_PANGENOME_HAPLOTYPE_NUMBER + 5, DeepVariant's convention for a graph with that many haplotypes. If passing in an already-sampled DV_PANGENOME_GBZ instead, set this explicitly to (haplotype count + 5); leaving it unset then gets DeepVariant's own default, which is tuned for the un-sampled reference pangenome."
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB: "(OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set."
        DV_PANGENOME_REFERENCE_PREFIX: "(OPTIONAL) Prefix on chromosome names in the pangenome GBZ (like 'GRCh38.') that isn't on the corresponding names in the BAM, analogous to REFERENCE_PREFIX but for the pangenome reference instead of the calling reference. Empty by default."
        DV_PANGENOME_REF_NAME: "(OPTIONAL) The name of the reference to keep in the pangenome gbz file for pangenome-aware DV; all other reference-sense paths are removed before calling. Required if DV_PANGENOME_GBZ is set."
        DV_PANGENOME_HAPLOTYPE_SAMPLING: "Should haplotype sampling of DV_PANGENOME_GBZ be done before pangenome-aware DV calling? This is a separate round of sampling from HAPLOTYPE_SAMPLING, which (if used) samples GBZ_FILE before mapping. Default is 'false'."
        DV_PANGENOME_DIPLOID_SAMPLING: "Should the DV_PANGENOME_HAPLOTYPE_SAMPLING round of haplotype sampling be done in diploid mode? Default is 'false'."
        DV_PANGENOME_HAPLOTYPE_NUMBER: "Number of haplotypes to sample for DV_PANGENOME_HAPLOTYPE_SAMPLING. Also used, if DV_PANGENOME_IMAGE_HEIGHT is not set, to size the pangenome-aware DV pileup images, so set it to the actual haplotype count even when passing in an already-sampled DV_PANGENOME_GBZ. Default is 32."
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 200)"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        DV_USE_GPUS: "Should DeepVariant use GPUs for calling variants? Default is 'true'."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+."
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        SPLIT_READ_MEM: "Memory, in GB, to use when splitting the reads into chunks. Default is 50."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling before running giraffe. Default is 'true'."
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        EVAL_CORES: "Number of cores to use when evaluating variant calls. Default is 8."
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
        Boolean EVALUATE_WITH_AARDVARK = false
        File? STRATIFICATION_ARCHIVE
        File? RESTRICT_REGIONS_BED
        String? TARGET_REGION
        Boolean RUN_STANDALONE_VCFEVAL = true
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Array[File]? DV_MODEL_FILES
        Array[File]? DV_MODEL_VARIABLES_FILES
        File? DV_PANGENOME_GBZ
        Int? DV_PANGENOME_IMAGE_HEIGHT
        Int? DV_PANGENOME_SHARED_MEMORY_SIZE_GB
        String? DV_PANGENOME_REFERENCE_PREFIX
        String? DV_PANGENOME_REF_NAME
        Boolean DV_PANGENOME_HAPLOTYPE_SAMPLING = false
        Boolean DV_PANGENOME_DIPLOID_SAMPLING = false
        Int DV_PANGENOME_HAPLOTYPE_NUMBER = 32
        Int HAPLOTYPE_INDEXING_MEM = 200
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        String OTHER_MAKEEXAMPLES_ARG = ""
        Boolean DV_USE_GPUS = true
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_MEM = 50
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Boolean HAPLOTYPE_SAMPLING = true
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = if MAP_MEM < 40 then MAP_MEM else 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int EVAL_CORES = 8
        Int EVAL_MEM = 60
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? VG_GIRAFFE_DOCKER
        String? VG_SURJECT_DOCKER

    }

    # Get the path names to operate on, and the FASTA reference.
    call reference_wf.PrepareReference {
        input:
        GBZ_FILE=GBZ_FILE,
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

    # If the same pangenome is used for mapping and for DeepVariant, the
    # .hapl file (and haplotype-sampled graph) made for mapping can be reused
    # for DeepVariant instead of being recomputed from scratch.
    Boolean pangenomes_are_same = defined(DV_PANGENOME_GBZ) && GBZ_FILE == select_first([DV_PANGENOME_GBZ])

    # Map the reads to GAF chunks.
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
        OUTPUT_SINGLE_BAM=false,
        OUTPUT_CALLING_BAMS=false,
        OUTPUT_GAF=OUTPUT_GAF,
        OUTPUT_GAF_CHUNKS=true,
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
        OUTPUT_HAPL=pangenomes_are_same,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
        INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
        KMER_COUNTING_MEM=KMER_COUNTING_MEM,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        VG_DOCKER=VG_DOCKER,
        VG_GIRAFFE_DOCKER=VG_GIRAFFE_DOCKER
    }

    # Surject the alignments, call variants, and compare to a truth set if one
    # was given.
    call gaf_wf.GiraffeDeepVariantFromGAF {
        input:
        GAF_CHUNKS=select_first([Giraffe.output_gaf_chunks]),
        GBZ_FILE=GBZ_FILE,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=OUTPUT_SINGLE_BAM,
        PAIRED_READS=PAIRED_READS,
        PATH_LIST_FILE=pipeline_path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=reference_file,
        REFERENCE_INDEX_FILE=reference_index_file,
        REFERENCE_DICT_FILE=reference_dict_file,
        HAPLOID_CONTIGS=HAPLOID_CONTIGS,
        PAR_REGIONS_BED_FILE=PAR_REGIONS_BED_FILE,
        PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
        MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
        EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
        EVALUATE_WITH_AARDVARK=EVALUATE_WITH_AARDVARK,
        STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
        RESTRICT_REGIONS_BED=RESTRICT_REGIONS_BED,
        TARGET_REGION=TARGET_REGION,
        RUN_STANDALONE_VCFEVAL=RUN_STANDALONE_VCFEVAL,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_MODEL_META=DV_MODEL_META,
        DV_MODEL_INDEX=DV_MODEL_INDEX,
        DV_MODEL_DATA=DV_MODEL_DATA,
        DV_MODEL_FILES=DV_MODEL_FILES,
        DV_MODEL_VARIABLES_FILES=DV_MODEL_VARIABLES_FILES,
        PANGENOME_GBZ=DV_PANGENOME_GBZ,
        DV_PANGENOME_IMAGE_HEIGHT=DV_PANGENOME_IMAGE_HEIGHT,
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB=DV_PANGENOME_SHARED_MEMORY_SIZE_GB,
        DV_PANGENOME_REFERENCE_PREFIX=DV_PANGENOME_REFERENCE_PREFIX,
        DV_PANGENOME_REF_NAME=DV_PANGENOME_REF_NAME,
        DV_PANGENOME_HAPLOTYPE_SAMPLING=DV_PANGENOME_HAPLOTYPE_SAMPLING,
        DV_PANGENOME_DIPLOID_SAMPLING=DV_PANGENOME_DIPLOID_SAMPLING,
        DV_PANGENOME_HAPLOTYPE_NUMBER=DV_PANGENOME_HAPLOTYPE_NUMBER,
        DV_PANGENOME_READS_FOR_SAMPLING_1=INPUT_READ_FILE_1,
        DV_PANGENOME_READS_FOR_SAMPLING_2=INPUT_READ_FILE_2,
        DV_PANGENOME_HAPL_FILE=select_first([Giraffe.output_hapl_file, HAPL_FILE]),
        DV_PANGENOME_DIST_FILE=DIST_FILE,
        DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES=MAP_CORES,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        VG_CORES=MAP_CORES,
        VG_MEM=MAP_MEM,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        REALIGN_MEM=REALIGN_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        MAKE_EXAMPLES_CORES=MAKE_EXAMPLES_CORES,
        MAKE_EXAMPLES_MEM=MAKE_EXAMPLES_MEM,
        EVAL_CORES=EVAL_CORES,
        EVAL_MEM=EVAL_MEM,
        VG_DOCKER=VG_DOCKER,
        VG_SURJECT_DOCKER=VG_SURJECT_DOCKER
    }

    output {
        File? output_vcfeval_evaluation_archive = GiraffeDeepVariantFromGAF.output_vcfeval_evaluation_archive
        File? output_happy_evaluation_archive = GiraffeDeepVariantFromGAF.output_happy_evaluation_archive
        File? output_aardvark_summary = GiraffeDeepVariantFromGAF.output_aardvark_summary
        Array[File]? output_aardvark_all_files = GiraffeDeepVariantFromGAF.output_aardvark_all_files
        File output_vcf = GiraffeDeepVariantFromGAF.output_vcf
        File output_vcf_index = GiraffeDeepVariantFromGAF.output_vcf_index
        File output_gvcf = GiraffeDeepVariantFromGAF.output_gvcf
        File output_gvcf_index = GiraffeDeepVariantFromGAF.output_gvcf_index
        File? output_gaf = Giraffe.output_gaf
        File? output_bam = GiraffeDeepVariantFromGAF.output_bam
        File? output_bam_index = GiraffeDeepVariantFromGAF.output_bam_index
        Array[File]? output_calling_bams = GiraffeDeepVariantFromGAF.output_calling_bams
        Array[File]? output_calling_bam_indexes = GiraffeDeepVariantFromGAF.output_calling_bam_indexes
        File? output_unmapped_bam = GiraffeDeepVariantFromGAF.output_unmapped_bam
    }
}

