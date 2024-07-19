version 1.0

import "../tasks/variant_evaluation.wdl" as eval
import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/deepvariant.wdl" as dv

workflow DeepVariant {

    meta {
        description: "## DeepVariant workflow \n Partial workflow to go from mapped reads (BAM) to small variant calls (VCF). Reads are pre-processed (e.g. indel realignment). DeepVariant then calls small variants. Includes optional comparison to a truth set."
    }

    parameter_meta {
        MERGED_BAM_FILE: "The all-contigs sorted BAM to call with."
        MERGED_BAM_FILE_INDEX: "The .bai index for the input BAM file"
        SAMPLE_NAME: "The sample name"
        PATH_LIST_FILE: "Text file where each line is a path name."
        REFERENCE_PREFIX: "Remove this off the beginning of path names to get contig names in the BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "FASTA reference to call against."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths."
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is 1. If null, uses DeepVariant default for the model type."
        TRUTH_VCF: "Path to .vcf.gz to compare against"
        TRUTH_VCF_INDEX: "Path to Tabix index for TRUTH_VCF"
        EVALUATION_REGIONS_BED: "BED to restrict comparison against TRUTH_VCF to"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? Default is 'true'."
        DV_NORM_READS: "Should DV normalize reads itself? Default is 'false'."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        REALIGN_MEM: "Memory, in GB, to use when realigning the reads. Default is 40."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
    }

    input {
        File MERGED_BAM_FILE
        File MERGED_BAM_FILE_INDEX
        String SAMPLE_NAME
        File PATH_LIST_FILE
        String REFERENCE_PREFIX = ""
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int? MIN_MAPQ = 1
        File? TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File? EVALUATION_REGIONS_BED
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Boolean DV_KEEP_LEGACY_AC = true
        Boolean DV_NORM_READS = false
        String OTHER_MAKEEXAMPLES_ARG = ""
        Int SPLIT_READ_CORES = 8
        Int REALIGN_MEM = 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50
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
        in_prefix_to_strip=REFERENCE_PREFIX
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
                in_min_mapq=MIN_MAPQ,
                in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                in_norm_reads=DV_NORM_READS,
                in_other_makeexamples_arg=OTHER_MAKEEXAMPLES_ARG,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }
        call dv.runDeepVariantCallVariants {
            input:
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=runDeepVariantMakeExamples.examples_file,
                in_nonvariant_site_tf_file=runDeepVariantMakeExamples.nonvariant_site_tf_file,
                in_model_type=DV_MODEL_TYPE,
                in_model_meta_file=DV_MODEL_META,
                in_model_index_file=DV_MODEL_INDEX,
                in_model_data_file=DV_MODEL_DATA,
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
    
        # To evaluate the VCF we need a template of the reference
        call eval.buildReferenceTemplate {
            input:
                in_reference_file=reference_file
        }
        
        # Direct vcfeval comparison makes an archive with FP and FN VCFs
        call eval.compareCalls {
            input:
                in_sample_vcf_file=concatClippedVCFChunks.output_merged_vcf,
                in_sample_vcf_index_file=concatClippedVCFChunks.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_template_archive=buildReferenceTemplate.output_template_archive,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_mem=CALL_MEM
        }
        
        # Hap.py comparison makes accuracy results stratified by SNPs and indels
        call eval.compareCallsHappy {
            input:
                in_sample_vcf_file=concatClippedVCFChunks.output_merged_vcf,
                in_sample_vcf_index_file=concatClippedVCFChunks.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_template_archive=buildReferenceTemplate.output_template_archive,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_mem=CALL_MEM
        }
    }

    output {
        File? output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File? output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_vcf = concatClippedVCFChunks.output_merged_vcf
        File output_vcf_index = concatClippedVCFChunks.output_merged_vcf_index
        File output_gvcf = concatClippedGVCFChunks.output_merged_vcf
        File output_gvcf_index = concatClippedGVCFChunks.output_merged_vcf_index
        Array[File] output_calling_bams = calling_bam
        Array[File] output_calling_bam_indexes = calling_bam_index
        File output_unmapped_bam = splitBAMbyPath.bam_unmapped_file
    }

}
