version 1.0

import "../tasks/bioinfo_utils.wdl" as utils
import "../tasks/deepvariant.wdl" as dv
import "aardvark_evaluation.wdl" as aardvark
import "happy_evaluation.wdl" as happy
import "../tasks/vg_map_hts.wdl" as map
import "./haplotype_sampling.wdl" as hapl

workflow DeepVariant {

    meta {
        description: "## DeepVariant workflow \n Partial workflow to go from mapped reads (BAM) to small variant calls (VCF). Reads are pre-processed (e.g. indel realignment). DeepVariant then calls small variants. Includes optional comparison to a truth set. More information at [https://github.com/vgteam/vg_wdl/tree/master#deepvariant-workflow](https://github.com/vgteam/vg_wdl/tree/master#deepvariant-workflow)."
    }

    parameter_meta {
        MERGED_BAM_FILE: "The all-contigs sorted BAM to call with."
        MERGED_BAM_FILE_INDEX: "The .bai index for the input BAM file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted by default. Default is 'false'."
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
        PANGENOME_GBZ: "(OPTIONAL) Path to a pangenome graph in GBZ format"
        DV_PANGENOME_IMAGE_HEIGHT: "(OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware models. It will be used only if PANGENOME_GBZ is set."
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB: "(OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set."
        DV_PANGENOME_REF_CHROM_PREFIX: "(OPTIONAL) The prefix to add to the chromosome name in the pangenome gbz file. It is empty by default. However sometimes we need to add a prefix (like 'GRCh38.') to the chromosome name in the pangenome gbz file to match the chromosome name in the BAM file. It is empty by default."
        DV_PANGENOME_REF_NAME: "The name of the reference in the pangenome gbz file. It is 'GRCh38' by default."
        REF_NAME_TO_REMOVE_FROM_PANGENOME: "(OPTIONAL) Name of the reference to remove from the graph before haplotype sampling"
        DV_PANGENOME_HAPLOTYPE_SAMPLING: "Should haplotype sampling be done on PANGENOME_GBZ? Default is 'false'."
        READS_FOR_SAMPLING_1: "(OPTIONAL) First input read file for haplotype sampling"
        READS_FOR_SAMPLING_2: "(OPTIONAL) Second input read file for haplotype sampling (if paired)"    
        PAIRED_READS_FOR_SAMPLING: "Are the reads for sampling paired-end?, Default is 'false'."
        DV_PANGENOME_DIPLOID_SAMPLING: "Should haplotype sampling be done in diploid mode? Default is 'false'."
        DV_PANGENOME_HAPLOTYPE_NUMBER: "Number of haplotypes for haplotype sampling. Default is 32."
        DV_PANGENOME_HAPL_FILE: "(OPTIONAL) Path to .hapl file used in haplotype sampling"
        DV_PANGENOME_DIST_FILE: "(OPTIONAL) Path to .dist file used in haplotype sampling"
        DV_PANGENOME_R_INDEX_FILE: "(OPTIONAL) Path to .ri file used in haplotype sampling"
        DV_PANGENOME_KFF_FILE: "(OPTIONAL) Path to .kff file used in haplotype sampling"
        INDEX_MINIMIZER_WEIGHTED: "Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)"
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)"
        KMER_COUNTING_MEM: "Memory, in GB, to use when counting kmers. (Default: 64)"
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)"
        HAPLOTYPE_SAMPLE_CORES: "Number of cores to use for haplotype sampling. Default is 16."
        GIRAFFE_PRESET: "(OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10). It is used for haplotype sampling the input pangenome. Default is 'default'."  
        CREATE_INDEX_OPTIONS_BEFORE_SAMPLING: "Additional options to pass to vg index before haplotype sampling. It is recommended to set --snarl-limit 1 to speed up index creation."
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads."
        DV_NORM_READS: "Should DV normalize reads itself? If unspecified this is not done, unless set in the model."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
        DV_USE_GPUS: "Should DeepVariant use GPUs for calling variants? Default is 'true'."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+."
        VG_DOCKER: "Container image to use when running vg (used for pangenome haplotype sampling)"
        BAM_PREPROCESS_MEM: "Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20."
        REALIGN_MEM: "Memory, in GB, to use for Abra indel realignment. Default is 40."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        MAKE_EXAMPLES_CORES: "Number of cores to use when making DeepVariant examples. Default is CALL_CORES."
        MAKE_EXAMPLES_MEM: "Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM."
        EVAL_CORES: "Number of cores to use when evaluating variant calls. Default is 8."
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
        Boolean EVALUATE_WITH_AARDVARK = false
        File? STRATIFICATION_ARCHIVE
        File? RESTRICT_REGIONS_BED
        String? TARGET_REGION
        Boolean RUN_STANDALONE_VCFEVAL = true
        String DV_MODEL_TYPE = "WGS"
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Array[File] DV_MODEL_FILES = select_all([DV_MODEL_META, DV_MODEL_INDEX, DV_MODEL_DATA])
        Array[File] DV_MODEL_VARIABLES_FILES = []
        File? PANGENOME_GBZ
        Int? DV_PANGENOME_IMAGE_HEIGHT
        Int? DV_PANGENOME_SHARED_MEMORY_SIZE_GB
        String? REF_NAME_TO_REMOVE_FROM_PANGENOME
        String? DV_PANGENOME_REF_CHROM_PREFIX
        String DV_PANGENOME_REF_NAME = "GRCh38"
        Boolean DV_PANGENOME_HAPLOTYPE_SAMPLING = false
        File? READS_FOR_SAMPLING_1
        File? READS_FOR_SAMPLING_2
        Boolean PAIRED_READS_FOR_SAMPLING = false
        Boolean DV_PANGENOME_DIPLOID_SAMPLING = false
        File? DV_PANGENOME_HAPL_FILE
        File? DV_PANGENOME_DIST_FILE
        File? DV_PANGENOME_R_INDEX_FILE
        File? DV_PANGENOME_KFF_FILE
        Int DV_PANGENOME_HAPLOTYPE_NUMBER = 32
        Boolean INDEX_MINIMIZER_WEIGHTED = true
        Int INDEX_MINIMIZER_MEM = if INDEX_MINIMIZER_WEIGHTED then 320 else 120
        Int KMER_COUNTING_MEM = 64
        Int HAPLOTYPE_INDEXING_MEM = 120
        Int HAPLOTYPE_SAMPLE_CORES = 16
        String GIRAFFE_PRESET = "default"
        String CREATE_INDEX_OPTIONS_BEFORE_SAMPLING = ""
        Boolean? DV_KEEP_LEGACY_AC
        Boolean? DV_NORM_READS
        String OTHER_MAKEEXAMPLES_ARG = ""
        Boolean DV_USE_GPUS = true
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        String VG_DOCKER = "quay.io/vgteam/vg:v1.68.0"
        Int BAM_PREPROCESS_MEM = 20
        Int REALIGN_MEM = 40
        Int CALL_CORES = 8
        Int CALL_MEM = 50 + (if defined(PANGENOME_GBZ) then CALL_CORES * 2 else 0)
        Int MAKE_EXAMPLES_CORES = CALL_CORES
        Int MAKE_EXAMPLES_MEM = CALL_MEM
        Int EVAL_CORES = 8
        Int EVAL_MEM = 60
    }

    call utils.uncompressReferenceIfNeeded {
        input:
        # We know REFERENCE_FILE is defined but the WDL type system doesn't.
        in_reference_file=REFERENCE_FILE,
        in_uncompress_cores=CALL_CORES
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
        strip_from_bam=REFERENCE_PREFIX_ON_BAM,
        thread_count=CALL_CORES,
        mem_gb=BAM_PREPROCESS_MEM
    }

    ##
    ## Haplotype sampling and then remove unnecessary sample (like CHM13 if the analysis is in GRCh38 coordinate space) from GBZ
    ## 

    if (DV_PANGENOME_HAPLOTYPE_SAMPLING && defined(PANGENOME_GBZ)) {
        call hapl.HaplotypeSampling {
        input:
            GBZ_FILE=select_first([PANGENOME_GBZ]),
            INPUT_READ_FILE_FIRST=select_first([READS_FOR_SAMPLING_1]),
            # If we're not doing paired reads the result here is probably null.
            INPUT_READ_FILE_SECOND=READS_FOR_SAMPLING_2,
            HAPLOTYPE_NUMBER=DV_PANGENOME_HAPLOTYPE_NUMBER,
            HAPL_FILE=DV_PANGENOME_HAPL_FILE,
            DIST_FILE=DV_PANGENOME_DIST_FILE,
            R_INDEX_FILE=DV_PANGENOME_R_INDEX_FILE,
            KFF_FILE=DV_PANGENOME_KFF_FILE,
            DIPLOID=DV_PANGENOME_DIPLOID_SAMPLING,
            INDEX_MINIMIZER_K = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 29 else 31,
            INDEX_MINIMIZER_W = if GIRAFFE_PRESET == "default" || GIRAFFE_PRESET == "fast" then 11 else 50,
            INDEX_MINIMIZER_WEIGHTED=INDEX_MINIMIZER_WEIGHTED,
            CREATE_INDEX_OPTIONS_BEFORE_SAMPLING=CREATE_INDEX_OPTIONS_BEFORE_SAMPLING,
            CORES=HAPLOTYPE_SAMPLE_CORES,
            KMER_COUNTING_MEM=KMER_COUNTING_MEM,
            HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
            INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
            SKIP_DIST_GENERATION=true,
            SKIP_MIN_GENERATION=true,
            VG_DOCKER=VG_DOCKER
        }

    }

    # The reason to remove unnecessary reference like CHM13 from the graph is that it's not haplotype sampled, and when we 
    # run pangenome-aware DV it would be included in pangenome section of the image though it is not relevant to the sample.
    if (defined(REF_NAME_TO_REMOVE_FROM_PANGENOME) && defined(PANGENOME_GBZ)) {
        call map.removeSampleFromGraph {
            input:
                in_gbz_file=select_first([HaplotypeSampling.sampled_graph, PANGENOME_GBZ]),
                sample_name=select_first([REF_NAME_TO_REMOVE_FROM_PANGENOME]),
                docker_image=VG_DOCKER
        }
    }

    if (defined(PANGENOME_GBZ)) {
            File DV_PANGENOME_GBZ = select_first([removeSampleFromGraph.output_graph_gbz,HaplotypeSampling.sampled_graph, PANGENOME_GBZ])
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
                thread_count=CALL_CORES,
                mem_gb=BAM_PREPROCESS_MEM
            }
            call utils.runAbraRealigner {
                input:
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    threadCount=CALL_CORES,
                    memoryGb=REALIGN_MEM
            }
        }
        File calling_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
        File calling_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
        # Here the name of the shared memory segment is based on the name of the BAM file. This is just 
        # to make sure that the shared memory segment name is unique for each contig to avoid memory 
        # interference between the DV jobs being run in the same file system.
        String DV_PANGENOME_SHARED_MEMORY_NAME = "SHARED_MEM_GBZ_" + sub(basename(calling_bam), "(left_shifted\\.)?(indel_realigned)?\\.bam$", "")  

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
                in_pangenome_gbz_file=DV_PANGENOME_GBZ,
                in_pangenome_height=DV_PANGENOME_IMAGE_HEIGHT,
                in_pangenome_shared_memory_name=DV_PANGENOME_SHARED_MEMORY_NAME,
                in_pangenome_shared_memory_size_gb=DV_PANGENOME_SHARED_MEMORY_SIZE_GB,
                in_pangenome_ref_chrom_prefix=DV_PANGENOME_REF_CHROM_PREFIX,
                in_pangenome_ref_name=DV_PANGENOME_REF_NAME,
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
                in_use_gpus=DV_USE_GPUS,
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
        if (!EVALUATE_WITH_AARDVARK) {
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
                    EVAL_CORES=EVAL_CORES,
                    EVAL_MEM=EVAL_MEM
            }
        }
        if (EVALUATE_WITH_AARDVARK) {
            
            # Aardvark takes one mandatory BED to actually run on.
            # We've decided this means EVALUATION_REGIONS_BED is now mandatory.
            # But we still need to intersect in RESTRICT_REGIONS_BED if provided.
            if (defined(RESTRICT_REGIONS_BED)) {
                call utils.intersectBeds as makeAardvarkBed {
                    input:
                        in_bed_1 = select_first([EVALUATION_REGIONS_BED]),
                        in_bed_2 = select_first([RESTRICT_REGIONS_BED])
                }
            }

            call aardvark.AardvarkEvaluation {
                input:
                    QUERY_VCF=concatClippedVCFChunks.output_merged_vcf,
                    QUERY_VCF_INDEX=concatClippedVCFChunks.output_merged_vcf_index,
                    TRUTH_VCF=select_first([TRUTH_VCF]),
                    TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
                    REFERENCE_FILE=reference_file,
                    REFERENCE_INDEX_FILE=reference_index_file,
                    REGIONS_BED=select_first([makeAardvarkBed.output_bed_file, EVALUATION_REGIONS_BED]),
                    STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
                    SAMPLE_NAME=SAMPLE_NAME,
                    THREADS=EVAL_CORES,
                    EVAL_MEM=EVAL_MEM
            }
        }
    }

    if (OUTPUT_SINGLE_BAM) {
        call utils.mergeAlignmentBAMChunks as mergeBAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=flatten([calling_bam, [splitBAMbyPath.bam_unmapped_file]]),
            in_cores=CALL_CORES,
            mem_gb=BAM_PREPROCESS_MEM
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
        File? output_aardvark_summary = AardvarkEvaluation.aardvark_summary
        Array[File]? output_aardvark_all_files = AardvarkEvaluation.aardvark_all_files
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

