version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/deepvariant.wdl" as dv

workflow GiraffeDeepVariantFromGAF {
    meta {
        description: "## Giraffe-DeepVariant-fromGAF workflow \n Surject a GAF and call small variants with DeepVariant. More information at [https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-from-gaf-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-from-gaf-workflow)."
    }

    parameter_meta {
        INPUT_GAF: "Input gzipped GAF file"
        GBZ_FILE: "Path to .gbz index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling bams' (one per contig) won't be outputed. Default is 'true'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set"
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is 1"
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        DV_MODEL_FILES: "Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format"
        DV_MODEL_VARIABLES_FILES: "Array of files that need to go in a 'variables' subdirectory for a DV model"
        DV_IS_1_7_OR_NEWER: "Flag to use DeepVariant 1.7+ command line syntax and recommended flags. Must be true if providing a DV 1.7+ Docker image, and false if providing an older one."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs"
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? Default is 'true'."
        DV_NORM_READS: "Should DV normalize reads itself? Default is 'fasle'."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        VG_CORES: "Number of cores to use when projecting the reads. Default is 16."
        VG_MEM: "Memory, in GB, to use when projecting the reads. Default is 120."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        VG_DOCKER: "Container image to use when running vg"
    }

    input {
        File INPUT_GAF
        File GBZ_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean PAIRED_READS = true
        File? PATH_LIST_FILE
        Array[String]+? CONTIGS
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int MIN_MAPQ = 1
        Int MAX_FRAGMENT_LENGTH = 3000
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        Array[File]? DV_MODEL_FILES
        Array[File]? DV_MODEL_VARIABLES_FILES
        Boolean? DV_IS_1_7_OR_NEWER
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Boolean DV_KEEP_LEGACY_AC = true
        Boolean DV_NORM_READS = false
        String OTHER_MAKEEXAMPLES_ARG = ""
        Int VG_CORES = 16
        Int VG_MEM = 120
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
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
                    in_extract_mem=VG_MEM,
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
            in_extract_mem=VG_MEM,
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

    call gautils.surjectGAFtoBAM {
        input:
        in_gaf_file=INPUT_GAF,
        in_gbz_file=GBZ_FILE,
        in_path_list_file=pipeline_path_list_file,
        in_sample_name=SAMPLE_NAME,
        in_max_fragment_length=MAX_FRAGMENT_LENGTH,
        in_paired_reads=PAIRED_READS,
        nb_cores=VG_CORES,
        mem_gb=VG_MEM,
        vg_docker=VG_DOCKER
    }

    call utils.sortBAM {
        input:
        in_bam_file=surjectGAFtoBAM.output_bam_file,
        in_ref_dict=reference_dict_file,
        in_prefix_to_strip=REFERENCE_PREFIX
    }
    
    # Split merged alignment by contigs list
    call utils.splitBAMbyPath {
        input:
        in_sample_name=SAMPLE_NAME,
        in_merged_bam_file=sortBAM.sorted_bam,
        in_merged_bam_file_index=sortBAM.sorted_bam_index,
        in_path_list_file=pipeline_path_list_file,
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
                    in_reference_index_file=reference_index_file
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
                in_min_mapq=MIN_MAPQ,
                in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                in_norm_reads=DV_NORM_READS,
                in_other_makeexamples_arg=OTHER_MAKEEXAMPLES_ARG,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM,
                in_dv_is_1_7_or_newer=DV_IS_1_7_OR_NEWER,
                in_dv_container=DV_NO_GPU_DOCKER,
                in_model_type=DV_MODEL_TYPE,
                in_model_files=DV_MODEL_FILES,
                in_model_variables_files=DV_MODEL_VARIABLES_FILES

        }
        call dv.runDeepVariantCallVariants {
            input:
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=runDeepVariantMakeExamples.examples_file,
                in_nonvariant_site_tf_file=runDeepVariantMakeExamples.nonvariant_site_tf_file,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM,
                in_dv_gpu_container=DV_GPU_DOCKER,
                in_model_type=DV_MODEL_TYPE,
                in_model_files=DV_MODEL_FILES,
                in_model_variables_files=DV_MODEL_VARIABLES_FILES
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

    if (OUTPUT_SINGLE_BAM){
        call utils.mergeAlignmentBAMChunks as mergeBAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=select_all(flatten([calling_bam, [splitBAMbyPath.bam_unmapped_file]]))
        }
    }

    if (!OUTPUT_SINGLE_BAM){
        Array[File] output_calling_bam_files = calling_bam
        Array[File] output_calling_bam_index_files = calling_bam_index
    }
    
    output {
        File output_vcf = concatClippedVCFChunks.output_merged_vcf
        File output_vcf_index = concatClippedVCFChunks.output_merged_vcf_index
        File output_gvcf = concatClippedGVCFChunks.output_merged_vcf
        File output_gvcf_index = concatClippedGVCFChunks.output_merged_vcf_index
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        Array[File]? output_calling_bams = output_calling_bam_files
        Array[File]? output_calling_bam_indexes = output_calling_bam_index_files
    }   
}
