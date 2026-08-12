version 1.0

import "../tasks/gam_gaf_utils.wdl" as gautils
import "./deepvariant.wdl" as dv_wf
import "./internal/prepare_reference.wdl" as reference_wf
import "./internal/surject.wdl" as surject_wf

workflow GiraffeDeepVariantFromGAF {
    meta {
        description: "## Giraffe-DeepVariant-fromGAF workflow \n Surject a GAF and call small variants with DeepVariant, optionally comparing the calls to a truth set. More information at [https://github.com/vgteam/vg_wdl/tree/master#giraffe-deepvariant-from-gaf-workflow](https://github.com/vgteam/vg_wdl/tree/master#giraffe-deepvariant-from-gaf-workflow)."
    }

    parameter_meta {
        INPUT_GAF: "(OPTIONAL) Input gzipped GAF file, which is split up to surject in parallel. Give this or GAF_CHUNKS."
        GAF_CHUNKS: "(OPTIONAL) Input gzipped GAF, already split into chunks that can be surjected in parallel. Give this or INPUT_GAF."
        READS_PER_CHUNK: "Number of reads to put in each chunk when splitting INPUT_GAF. Unused if GAF_CHUNKS is given. Default 20 million."
        GBZ_FILE: "Path to .gbz index file. Has to be the graph the reads were mapped to, since the alignments name its nodes."
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted. Default is 'true'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM."
        OUTPUT_UNMAPPED_BAM: "Should an unmapped reads BAM be saved? Default is false."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set"
        HAPLOID_CONTIGS: "(OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5."
        PAR_REGIONS_BED_FILE: "(OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5."
        PRUNE_LOW_COMPLEXITY: "Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is 1"
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        SURJECT_OPTIONS: "Extra command line options for vg surject."
        TRUTH_VCF: "(OPTIONAL) Path to .vcf.gz to compare the calls against. Evaluation only runs if this and TRUTH_VCF_INDEX are given."
        TRUTH_VCF_INDEX: "(OPTIONAL) Path to Tabix index for TRUTH_VCF"
        EVALUATION_REGIONS_BED: "(OPTIONAL) BED to evaluate against TRUTH_VCF on, where false positives will be counted. Required when EVALUATE_WITH_AARDVARK is set."
        EVALUATE_WITH_AARDVARK: "Should the calls be compared to TRUTH_VCF with Aardvark instead of hap.py? Default is 'false'."
        STRATIFICATION_ARCHIVE: "(OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the results down by. Only used when EVALUATE_WITH_AARDVARK is set."
        RESTRICT_REGIONS_BED: "(OPTIONAL) BED to restrict comparison against TRUTH_VCF to"
        TARGET_REGION: "(OPTIONAL) Contig or region to restrict evaluation to"
        RUN_STANDALONE_VCFEVAL: "Whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: "(OPTIONAL) .meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: "(OPTIONAL) .index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: "(OPTIONAL) .data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        DV_MODEL_FILES: "Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format"
        DV_MODEL_VARIABLES_FILES: "Array of files that need to go in a 'variables' subdirectory for a DV model"
        PANGENOME_GBZ: "(OPTIONAL) Path to a pangenome graph in GBZ format for pangenome-aware DV. All reference-sense paths in it other than DV_PANGENOME_REF_NAME are removed before calling, since they are not part of this sample and would otherwise show up as uninformative extra tracks in the pileup images."
        DV_PANGENOME_IMAGE_HEIGHT: "(OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware models. It will be used only if PANGENOME_GBZ is set. If DV_PANGENOME_HAPLOTYPE_SAMPLING is done by this workflow and this is not set, it defaults to DV_PANGENOME_HAPLOTYPE_NUMBER + 5, DeepVariant's convention for a graph with that many haplotypes. If passing in an already-sampled PANGENOME_GBZ instead, set this explicitly to (haplotype count + 5); leaving it unset then gets DeepVariant's own default, which is tuned for the un-sampled reference pangenome."
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB: "(OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set."
        DV_PANGENOME_REFERENCE_PREFIX: "(OPTIONAL) Prefix on chromosome names in the pangenome GBZ (like 'GRCh38.') that isn't on the corresponding names in the BAM, analogous to REFERENCE_PREFIX but for the pangenome reference instead of the calling reference. Empty by default."
        DV_PANGENOME_REF_NAME: "(OPTIONAL) The name of the reference to keep in the pangenome gbz file for pangenome-aware DV; all other reference-sense paths are removed before calling. Required if PANGENOME_GBZ is set."
        DV_PANGENOME_HAPLOTYPE_SAMPLING: "Should haplotype sampling of PANGENOME_GBZ be done before pangenome-aware DV calling? Default is 'false'."
        DV_PANGENOME_READS_FOR_SAMPLING_1: "(OPTIONAL) First input read file for haplotype sampling"
        DV_PANGENOME_READS_FOR_SAMPLING_2: "(OPTIONAL) Second input read file for haplotype sampling (if paired)"
        DV_PANGENOME_DIPLOID_SAMPLING: "Should haplotype sampling be done in diploid mode? Default is 'false'."
        DV_PANGENOME_HAPLOTYPE_NUMBER: "Number of haplotypes to sample for haplotype sampling. Also used, if DV_PANGENOME_IMAGE_HEIGHT is not set, to size the pangenome-aware DV pileup images, so set it to the actual haplotype count even when passing in an already-sampled PANGENOME_GBZ. Default is 32."
        DV_PANGENOME_HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        DV_PANGENOME_DIST_FILE: "(OPTIONAL) Path to .dist file used in haplotype sampling"
        DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES: "Number of cores to use for haplotype sampling. Default is 16."
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling). (Default: 120)"
        DV_USE_GPUS: "Should DeepVariant use GPUs for calling variants? Default is 'true'."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        VG_CORES: "Number of cores to use when projecting the reads. Default is 16."
        VG_MEM: "Memory, in GB, to use when projecting the reads. Default is 120."
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        EVAL_CORES: "Number of cores to use when evaluating variant calls. Default is 8."
        EVAL_MEM: "Memory, in GB, to use when evaluating variant calls. Default is 60."
        VG_DOCKER: "Container image to use when running vg"
        VG_SURJECT_DOCKER: "(OPTIONAL) Alternate container image to use when running vg surject"
    }

    input {
        File? INPUT_GAF
        Array[File]? GAF_CHUNKS
        Int READS_PER_CHUNK = 20000000
        File GBZ_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean OUTPUT_CALLING_BAMS = !OUTPUT_SINGLE_BAM
        Boolean OUTPUT_UNMAPPED_BAM = false
        Boolean PAIRED_READS = true
        File? PATH_LIST_FILE
        Array[String]+? CONTIGS
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
        Int MIN_MAPQ = 1
        Int MAX_FRAGMENT_LENGTH = 3000
        String SURJECT_OPTIONS = ""
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
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        Array[File]? DV_MODEL_FILES
        Array[File]? DV_MODEL_VARIABLES_FILES
        File? PANGENOME_GBZ
        Int? DV_PANGENOME_IMAGE_HEIGHT
        Int? DV_PANGENOME_SHARED_MEMORY_SIZE_GB
        String? DV_PANGENOME_REFERENCE_PREFIX
        String? DV_PANGENOME_REF_NAME
        Boolean DV_PANGENOME_HAPLOTYPE_SAMPLING = false
        File? DV_PANGENOME_READS_FOR_SAMPLING_1
        File? DV_PANGENOME_READS_FOR_SAMPLING_2
        Boolean DV_PANGENOME_DIPLOID_SAMPLING = false
        Int DV_PANGENOME_HAPLOTYPE_NUMBER = 32
        File? DV_PANGENOME_HAPL_FILE
        File? DV_PANGENOME_DIST_FILE
        Int DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES = 16
        Int HAPLOTYPE_INDEXING_MEM = 120
        Boolean DV_USE_GPUS = true
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        String OTHER_MAKEEXAMPLES_ARG = ""
        Int VG_CORES = 16
        Int VG_MEM = 120
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int EVAL_CORES = 8
        Int EVAL_MEM = 60
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String? VG_SURJECT_DOCKER
    }

    # Get the contigs to work on, and the reference.
    call reference_wf.PrepareReference {
        input:
        GBZ_FILE=GBZ_FILE,
        CONTIGS=CONTIGS,
        PATH_LIST_FILE=PATH_LIST_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_DICT_FILE=REFERENCE_DICT_FILE,
        EXTRACT_MEM=VG_MEM,
        VG_DOCKER=VG_DOCKER
    }

    if (!defined(GAF_CHUNKS)) {
        # We were handed one GAF, so we have to chunk it ourselves to surject in
        # parallel.
        call gautils.splitGAF {
            input:
            in_gaf_file=select_first([INPUT_GAF]),
            in_read_per_chunk=READS_PER_CHUNK
        }
    }

    call surject_wf.Surject {
        input:
        GAF_CHUNKS=select_first([GAF_CHUNKS, splitGAF.gaf_chunk_files]),
        GBZ_FILE=GBZ_FILE,
        PATH_LIST_FILE=PrepareReference.path_list_file,
        REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        SAMPLE_NAME=SAMPLE_NAME,
        PAIRED_READS=PAIRED_READS,
        PRUNE_LOW_COMPLEXITY=PRUNE_LOW_COMPLEXITY,
        MAX_FRAGMENT_LENGTH=MAX_FRAGMENT_LENGTH,
        SURJECT_OPTIONS=SURJECT_OPTIONS,
        SURJECT_CORES=VG_CORES,
        SURJECT_MEM=VG_MEM,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        VG_DOCKER=select_first([VG_SURJECT_DOCKER, VG_DOCKER])
    }

    call dv_wf.DeepVariant {
        input:
        MERGED_BAM_FILE=Surject.bam_file,
        MERGED_BAM_FILE_INDEX=Surject.bam_index_file,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_SINGLE_BAM=OUTPUT_SINGLE_BAM,
        OUTPUT_CALLING_BAMS=OUTPUT_CALLING_BAMS,
        OUTPUT_UNMAPPED_BAM=OUTPUT_UNMAPPED_BAM,
        PATH_LIST_FILE=PrepareReference.path_list_file,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        # We surjected to a BAM without the prefix.
        REFERENCE_PREFIX_ON_BAM=false,
        REFERENCE_FILE=PrepareReference.reference_file,
        REFERENCE_INDEX_FILE=PrepareReference.reference_index_file,
        REFERENCE_DICT_FILE=PrepareReference.reference_dict_file,
        HAPLOID_CONTIGS=HAPLOID_CONTIGS,
        PAR_REGIONS_BED_FILE=PAR_REGIONS_BED_FILE,
        LEFTALIGN_BAM=LEFTALIGN_BAM,
        REALIGN_INDELS=REALIGN_INDELS,
        REALIGNMENT_EXPANSION_BASES=REALIGNMENT_EXPANSION_BASES,
        MIN_MAPQ=MIN_MAPQ,
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
        PANGENOME_GBZ=PANGENOME_GBZ,
        DV_PANGENOME_IMAGE_HEIGHT=DV_PANGENOME_IMAGE_HEIGHT,
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB=DV_PANGENOME_SHARED_MEMORY_SIZE_GB,
        DV_PANGENOME_REFERENCE_PREFIX=DV_PANGENOME_REFERENCE_PREFIX,
        DV_PANGENOME_REF_NAME=DV_PANGENOME_REF_NAME,
        DV_PANGENOME_HAPLOTYPE_SAMPLING=DV_PANGENOME_HAPLOTYPE_SAMPLING,
        DV_PANGENOME_READS_FOR_SAMPLING_1=DV_PANGENOME_READS_FOR_SAMPLING_1,
        DV_PANGENOME_READS_FOR_SAMPLING_2=DV_PANGENOME_READS_FOR_SAMPLING_2,
        DV_PANGENOME_DIPLOID_SAMPLING=DV_PANGENOME_DIPLOID_SAMPLING,
        DV_PANGENOME_HAPLOTYPE_NUMBER=DV_PANGENOME_HAPLOTYPE_NUMBER,
        DV_PANGENOME_HAPL_FILE=DV_PANGENOME_HAPL_FILE,
        DV_PANGENOME_DIST_FILE=DV_PANGENOME_DIST_FILE,
        DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES=DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        DV_KEEP_LEGACY_AC=DV_KEEP_LEGACY_AC,
        DV_NORM_READS=DV_NORM_READS,
        OTHER_MAKEEXAMPLES_ARG=OTHER_MAKEEXAMPLES_ARG,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        VG_DOCKER=VG_DOCKER,
        BAM_PREPROCESS_MEM=BAM_PREPROCESS_MEM,
        REALIGN_MEM=REALIGN_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        MAKE_EXAMPLES_CORES=MAKE_EXAMPLES_CORES,
        MAKE_EXAMPLES_MEM=MAKE_EXAMPLES_MEM,
        EVAL_CORES=EVAL_CORES,
        EVAL_MEM=EVAL_MEM
    }

    output {
        File output_vcf = DeepVariant.output_vcf
        File output_vcf_index = DeepVariant.output_vcf_index
        File output_gvcf = DeepVariant.output_gvcf
        File output_gvcf_index = DeepVariant.output_gvcf_index
        File? output_vcfeval_evaluation_archive = DeepVariant.output_vcfeval_evaluation_archive
        File? output_happy_evaluation_archive = DeepVariant.output_happy_evaluation_archive
        File? output_aardvark_summary = DeepVariant.output_aardvark_summary
        Array[File]? output_aardvark_all_files = DeepVariant.output_aardvark_all_files
        File? output_bam = DeepVariant.output_bam
        File? output_bam_index = DeepVariant.output_bam_index
        Array[File]? output_calling_bams = DeepVariant.output_calling_bams
        Array[File]? output_calling_bam_indexes = DeepVariant.output_calling_bam_indexes
        File? output_unmapped_bam = DeepVariant.output_unmapped_bam
    }
}
