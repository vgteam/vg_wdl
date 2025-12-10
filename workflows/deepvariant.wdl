version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/deepvariant.wdl" as dv
import "happy_evaluation.wdl" as happy

workflow DeepVariant {

    meta {
        description: "## DeepVariant workflow \n Partial workflow to go from mapped reads (BAM) to small variant calls (VCF). Reads are pre-processed (e.g. indel realignment). DeepVariant then calls small variants. Includes optional comparison to a truth set."
    }

    parameter_meta {
        MERGED_BAM_FILE: "The all-contigs sorted BAM to call with."
        MERGED_BAM_FILE_INDEX: "The .bai index for the input BAM file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling bams' (one per contig) won't be outputed by default. Default is 'false'."
        OUTPUT_CALLING_BAMS: "Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM."
        OUTPUT_UNMAPPED_BAM: "Should an unmapped reads BAM be saved? Default is false."
        CONTIGS: "Contig path names to use as PATH_LIST_FILE. Must be set if PATH_LIST_FILE is not."
        PATH_LIST_FILE: "Text file where each line is a contig name to evaluate on. Must be set if CONTIGS is not."
        REFERENCE_PREFIX: "Remove this off the beginning of path names to get contig names in the BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_PREFIX_ON_BAM: "If true, the REFERENCE_PREFIX is also on the sequence names in the BAM header and needs to be removed."
        REFERENCE_FILE: "FASTA reference to call against."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        HAPLOID_CONTIGS: "(OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5."
        PAR_REGIONS_BED_FILE: "(OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'. If true, all input reads, including secondaries, must have the read sequence given."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'. If true, all input reads must be in a read group."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type."
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
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        EVAL_MEM: "Memory, in GB, to use when evaluating variant calls. Default is 60."
    }

    input {
        File MERGED_BAM_FILE
        File MERGED_BAM_FILE_INDEX
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = false
        Boolean OUTPUT_CALLING_BAMS = !OUTPUT_SINGLE_BAM
        Boolean OUTPUT_UNMAPPED_BAM = false
        Array[String]+? CONTIGS
        File PATH_LIST_FILE = write_lines(select_first([CONTIGS]))
        String REFERENCE_PREFIX = ""
        Boolean REFERENCE_PREFIX_ON_BAM = false
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Array[String]? HAPLOID_CONTIGS
        File? PAR_REGIONS_BED_FILE
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int? MIN_MAPQ
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
        Array[File] DV_MODEL_FILES = select_all([DV_MODEL_META, DV_MODEL_INDEX, DV_MODEL_DATA])
        Array[File] DV_MODEL_VARIABLES_FILES = []
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        String OTHER_MAKEEXAMPLES_ARG = ""
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int EVAL_MEM = 60
    }

    call utils.uncompressReferenceIfNeeded {
        input:
        # We know REFERENCE_FILE is defined but the WDL type system doesn't.
        in_reference_file=REFERENCE_FILE
    }
    File reference_file = uncompressReferenceIfNeeded.reference_file

    if (!defined(REFERENCE_INDEX_FILE) || !defined(REFERENCE_DICT_FILE)) {
        call utils.indexReference {
            input:
                in_reference_file=reference_file
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])

    # Split merged alignment by contigs list
    call utils.splitBAMbyPath {
        input:
        in_sample_name=SAMPLE_NAME,
        in_merged_bam_file=MERGED_BAM_FILE,
        in_merged_bam_file_index=MERGED_BAM_FILE_INDEX,
        in_path_list_file=PATH_LIST_FILE,
        in_prefix_to_strip=REFERENCE_PREFIX,
        strip_from_bam=REFERENCE_PREFIX_ON_BAM
    }

    ##
    ## Call variants with DeepVariant in each contig
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
                mem_gb=BAM_PREPROCESS_MEM
            }
            call utils.runAbraRealigner {
                input:
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    memoryGb=REALIGN_MEM
            }
        }
        File calling_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
        File calling_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
        ## DeepVariant calling
        call dv.runDeepVariantMakeExamples {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=calling_bam,
                in_bam_file_index=calling_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_model_type=DV_MODEL_TYPE,
                in_model_files=DV_MODEL_FILES,
                in_model_variables_files=DV_MODEL_VARIABLES_FILES,
                in_min_mapq=MIN_MAPQ,
                in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                in_norm_reads=DV_NORM_READS,
                in_other_makeexamples_arg=OTHER_MAKEEXAMPLES_ARG,
                in_dv_container=DV_NO_GPU_DOCKER,
                in_call_cores=MAKE_EXAMPLES_CORES,
                in_call_mem=MAKE_EXAMPLES_MEM
        }
        call dv.runDeepVariantCallVariants {
            input:
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=runDeepVariantMakeExamples.examples_file,
                in_nonvariant_site_tf_file=runDeepVariantMakeExamples.nonvariant_site_tf_file,
                in_model_type=DV_MODEL_TYPE,
                in_model_files=DV_MODEL_FILES,
                in_model_variables_files=DV_MODEL_VARIABLES_FILES,
                in_haploid_contigs=HAPLOID_CONTIGS,
                in_par_regions_bed_file=PAR_REGIONS_BED_FILE,
                in_dv_gpu_container=DV_GPU_DOCKER,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }
    }

    # Merge distributed variant called VCFs
    call utils.concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file
    }
    call utils.concatClippedVCFChunks as concatClippedGVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_gvcf_file
    }

    if (defined(TRUTH_VCF) && defined(TRUTH_VCF_INDEX)) {
        call happy.HappyEvaluation {
            input:
                VCF=concatClippedVCFChunks.output_merged_vcf,
                VCF_INDEX=concatClippedVCFChunks.output_merged_vcf_index,
                TRUTH_VCF=select_first([TRUTH_VCF]),
                TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
                REFERENCE_FILE=reference_file,
                REFERENCE_INDEX_FILE=reference_index_file,
                EVALUATION_REGIONS_BED=EVALUATION_REGIONS_BED,
                RESTRICT_REGIONS_BED=RESTRICT_REGIONS_BED,
                TARGET_REGION=TARGET_REGION,
                # Don't forward the reference prefix; we did it already on the BAMs.
                REMOVE_HOM_REFS=false,
                RUN_STANDALONE_VCFEVAL=RUN_STANDALONE_VCFEVAL,
                EVAL_MEM=EVAL_MEM
        }
    }

    if (OUTPUT_SINGLE_BAM) {
        call utils.mergeAlignmentBAMChunks as mergeBAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=flatten([calling_bam, [splitBAMbyPath.bam_unmapped_file]])
        }
    }

    if (OUTPUT_CALLING_BAMS) {
        Array[File] output_calling_bam_files = calling_bam
        Array[File] output_calling_bam_index_files = calling_bam_index
    }

    if (OUTPUT_UNMAPPED_BAM) {
        File output_unmapped_bam_file = splitBAMbyPath.bam_unmapped_file
    }

    output {
        File? output_vcfeval_evaluation_archive = HappyEvaluation.output_vcfeval_evaluation_archive
        File? output_happy_evaluation_archive = HappyEvaluation.output_happy_evaluation_archive
        File output_vcf = concatClippedVCFChunks.output_merged_vcf
        File output_vcf_index = concatClippedVCFChunks.output_merged_vcf_index
        File output_gvcf = concatClippedGVCFChunks.output_merged_vcf
        File output_gvcf_index = concatClippedGVCFChunks.output_merged_vcf_index
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        Array[File]? output_calling_bams = output_calling_bam_files
        Array[File]? output_calling_bam_indexes = output_calling_bam_index_files
        File? output_unmapped_bam = output_unmapped_bam_file
    }

}
