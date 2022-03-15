version 1.0

import "giraffe_and_deepvariant.wdl" as main

### giraffe_and_deepvariant.wdl ###
## Author: Charles Markello
## Description: Core VG Giraffe mapping and DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMap {
    input {
        File INPUT_CRAM_FILE                            # Input CRAM file
	File CRAM_REF                                   # Genome fasta file associated with the CRAM file
        File CRAM_REF_INDEX                             # Index of the fasta file associated with the CRAM file
        Int NB_CRAM_CHUNKS = 15                         # Number of chunks to split the reads in
        Int MAX_CRAM_CHUNKS = 15                        # Number of chunks to actually analyze
        String SAMPLE_NAME                              # The sample name
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.37.0" # VG Container used in the pipeline
        String GIRAFFE_OPTIONS = ""                     # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        String REFERENCE_PREFIX = ""                    # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File XG_FILE                                    # Path to .xg index file
        File GBWT_FILE                                  # Path to .gbwt index file
        File GGBWT_FILE                                 # Path to .gg index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
        File? TRUTH_VCF                                 # Path to .vcf.gz to compare against
        File? TRUTH_VCF_INDEX                           # Path to Tabix index for TRUTH_VCF
        File? EVALUATION_REGIONS_BED                    # BED to restrict comparison against TRUTH_VCF to
        File? DV_MODEL_META                             # .meta file for a custom DeepVariant calling model
        File? DV_MODEL_INDEX                            # .index file for a custom DeepVariant calling model
        File? DV_MODEL_DATA                             # .data-00000-of-00001 file for a custom DeepVariant calling model
        Boolean LEFTALIGN_BAM = true                    # Whether or not to left-align reads in the BAM before DV
        Boolean REALIGN_INDELS = true                   # Whether or not to realign reads near indels before DV
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        # DeepVariant tontainer to use for CPU steps
        String DV_CONTAINER = "google/deepvariant:1.3.0"
        # DeepVariant container to use for GPU steps
        String DV_GPU_CONTAINER = "google/deepvariant:1.3.0-gpu"
        Boolean DV_KEEP_LEGACY_AC = true                # Should DV use the legacy allele counter behavior?
        Boolean DV_NORM_READS = false                   # Should DV normalize reads itself?
        Boolean OUTPUT_GAF = true                       # Should a GAF file with the aligned reads be saved?
        Boolean OUTPUT_GAM = true                       # Should a GAM file with the aligned reads be saved?
        Int CRAM_CONVERT_CORES = 8
        Int CRAM_CONVERT_DISK = 60
        Int MAP_CORES = 16
        Int MAP_DISK = 200
        Int MAP_MEM = 120
        Int CALL_CORES = 8
        Int CALL_DISK = 40
        Int CALL_MEM = 50
        File? REFERENCE_FILE                            # (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.
        File? REFERENCE_INDEX_FILE                      # (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
        File? REFERENCE_DICT_FILE                       # (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. 
    }

    # Split the reads (CRAM file) in chunks (FASTQ)
    call splitCramFastq {
        input:
        in_cram_file=INPUT_CRAM_FILE,
        in_ref_file=CRAM_REF,
        in_ref_index_file=CRAM_REF_INDEX,
        in_cram_convert_cores=CRAM_CONVERT_CORES,
        in_cram_convert_disk=CRAM_CONVERT_DISK,
        in_nb_chunks=NB_CRAM_CHUNKS,
        in_max_chunks=MAX_CRAM_CHUNKS
    }

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call main.extractSubsetPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_extract_disk=MAP_DISK,
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
        call main.extractReference {
            input:
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_vg_container=VG_CONTAINER,
                in_extract_disk=MAP_DISK,
                in_extract_mem=MAP_MEM
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call main.indexReference {
            input:
                in_reference_file=reference_file,
                in_index_disk=MAP_DISK,
                in_index_mem=MAP_MEM
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])

    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    scatter (chunk_id in range(length(splitCramFastq.output_read_chunks_1))) {
        call main.runVGGIRAFFE {
            input:
                in_left_read_pair_chunk_file=splitCramFastq.output_read_chunks_1[chunk_id],
                in_right_read_pair_chunk_file=splitCramFastq.output_read_chunks_2[chunk_id],
                in_vg_container=VG_CONTAINER,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_xg_file=XG_FILE,
                in_gbwt_file=GBWT_FILE,
                in_ggbwt_file=GGBWT_FILE,
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
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        call main.surjectGAMtoBAM {
                input:
                in_gam_file=runVGGIRAFFE.chunk_gam_file,
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_sample_name=SAMPLE_NAME,
                in_max_fragment_length=MAX_FRAGMENT_LENGTH,
                in_vg_container=VG_CONTAINER,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
            }
        if (REFERENCE_PREFIX != "") {
            # use samtools to replace the header contigs with those from our dict.
            # this allows the header to contain contigs that are not in the graph,
            # which is more general and lets CHM13-based graphs be used to call on GRCh38
            # also, strip out contig prefixes in the BAM body
            call main.fixBAMContigNaming {
                input:
                    in_bam_file=surjectGAMtoBAM.chunk_bam_file,
                    in_ref_dict=reference_dict_file,
                    in_prefix_to_strip=REFERENCE_PREFIX,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        File properly_named_bam_file = select_first([fixBAMContigNaming.fixed_bam_file, surjectGAMtoBAM.chunk_bam_file]) 
        call main.sortBAMFile {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_chunk_file=properly_named_bam_file,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
    }
    Array[File] alignment_chunk_bam_files = select_all(sortBAMFile.sorted_chunk_bam)

    call main.mergeAlignmentBAMChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=alignment_chunk_bam_files,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }
    
    if (REFERENCE_PREFIX != "") {
        # strip all the GRCh38's off our path list file.  we need them for surject as they are in the path
        # but fixBAMContigNaming above stripped them, so we don't need them downstream
        call main.fixPathNames {
            input:
                in_path_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
        }
    }
    File properly_named_path_list_file = select_first([fixPathNames.fixed_path_list_file, pipeline_path_list_file])
             
    # Split merged alignment by contigs list
    call main.splitBAMbyPath { 
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
            in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
            in_path_list_file=properly_named_path_list_file,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }

    ##
    ## Call variants with DeepVariant in each contig
    ##
    scatter (deepvariant_caller_input_files in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        ## Evantually shift and realign reads
        if (LEFTALIGN_BAM){
            # Just left-shift each read individually
            call main.leftShiftBAMFile {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=deepvariant_caller_input_files.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_call_disk=CALL_DISK
            }
            # This tool can't make an index itself so we need to re-index the BAM
            call main.indexBAMFile {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=leftShiftBAMFile.left_shifted_bam,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
            }
        }
        if (REALIGN_INDELS) {
            File forrealign_bam = select_first([leftShiftBAMFile.left_shifted_bam, deepvariant_caller_input_files.left])
            File forrealign_index = select_first([indexBAMFile.bam_index, deepvariant_caller_input_files.right])
            # Do indel realignment
            call main.runGATKRealignerTargetCreator {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_reference_dict_file=reference_dict_file,
                    in_call_disk=CALL_DISK
            }
            if (REALIGNMENT_EXPANSION_BASES != 0) {
                # We want the realignment targets to be wider
                call main.widenRealignmentTargets {
                    input:
                        in_target_bed_file=runGATKRealignerTargetCreator.realigner_target_bed,
                        in_reference_index_file=reference_index_file,
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                        in_call_disk=CALL_DISK
                }
            }
            File target_bed_file = select_first([widenRealignmentTargets.output_target_bed_file, runGATKRealignerTargetCreator.realigner_target_bed])
            call main.runAbraRealigner {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_call_disk=CALL_DISK
            }
        }
        File calling_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.left_shifted_bam, deepvariant_caller_input_files.left])
        File calling_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, indexBAMFile.bam_index, deepvariant_caller_input_files.right])
        ## DeepVariant calling
        call main.runDeepVariantMakeExamples {
            input:
                in_dv_container=DV_CONTAINER,
                in_sample_name=SAMPLE_NAME,
                in_bam_file=calling_bam,
                in_bam_file_index=calling_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_min_mapq=MIN_MAPQ,
                in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                in_norm_reads=DV_NORM_READS,
                in_call_cores=CALL_CORES,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
        call main.runDeepVariantCallVariants {
            input:
                in_dv_gpu_container=DV_GPU_CONTAINER,
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=runDeepVariantMakeExamples.examples_file,
                in_nonvariant_site_tf_file=runDeepVariantMakeExamples.nonvariant_site_tf_file,
                in_model_meta_file=DV_MODEL_META,
                in_model_index_file=DV_MODEL_INDEX,
                in_model_data_file=DV_MODEL_DATA,
                in_call_cores=CALL_CORES,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
    }
    # Merge distributed variant called VCFs
    call main.concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call main.bgzipMergedVCF {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    
    if (defined(TRUTH_VCF) && defined(TRUTH_VCF_INDEX)) {
    
        # To evaluate the VCF we need a template of the reference
        call main.buildReferenceTemplate {
            input:
                in_reference_file=reference_file
        }
        
        # Direct vcfeval comparison makes an archive with FP and FN VCFs
        call main.compareCalls {
            input:
                in_sample_vcf_file=bgzipMergedVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipMergedVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_template_archive=buildReferenceTemplate.output_template_archive,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
        
        # Hap.py comparison makes accuracy results stratified by SNPs and indels
        call main.compareCallsHappy {
            input:
                in_sample_vcf_file=bgzipMergedVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipMergedVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
    }

    if (OUTPUT_GAF){
        scatter (gam_chunk_file in runVGGIRAFFE.chunk_gam_file) {
            call main.convertGAMtoGAF {
                input:
                in_xg_file=XG_FILE,
                in_gam_file= gam_chunk_file,
                in_vg_container=VG_CONTAINER,
                in_cores=MAP_CORES,
                in_disk=MAP_DISK,
                in_mem=MAP_MEM
            }
        }
        call main.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=convertGAMtoGAF.output_gaf,
            in_vg_container=VG_CONTAINER,
            in_disk=2*MAP_DISK
        }
    }

    if (OUTPUT_GAM){
        call main.mergeGAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gam_chunk_files=runVGGIRAFFE.chunk_gam_file,
            in_vg_container=VG_CONTAINER,
            in_disk=2*MAP_DISK
        }
    }
    
    output {
        File? output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File? output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_vcf = bgzipMergedVCF.output_merged_vcf
        File output_vcf_index = bgzipMergedVCF.output_merged_vcf_index
        File? output_gaf = mergeGAF.output_merged_gaf
        File? output_gam = mergeGAM.output_merged_gam
        Array[File] output_calling_bams = calling_bam
        Array[File] output_calling_bam_indexes = calling_bam_index
    }   
}

########################
### TASK DEFINITIONS ###
########################

# Retrieve chunks from a CRAM file and convert them to FASTQ
# CRAM + (k,N) -> FASTQ_1_of_N, FASTQ_2_of_N, ..., FASTQ_k_of_N
task splitCramFastq {
    input {
        File in_cram_file
        File in_ref_file
        File in_ref_index_file
        Int in_cram_convert_cores
        Int in_cram_convert_disk
        Int in_nb_chunks
        Int in_max_chunks
    }
    command <<< 
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        seq 0 ~{in_nb_chunks} | head -n ~{in_max_chunks} | parallel -j ~{in_cram_convert_cores} "samtools collate -k {} -K ~{in_nb_chunks} --reference ~{in_ref_file} -Ouf ~{in_cram_file} {} | samtools fastq -1 reads.{}.R1.fastq.gz -2 reads.{}.R2.fastq.gz -0 reads.{}.o.fq.gz -s reads.{}.s.fq.gz -c 1 -N -"
    >>>
    output {
        Array[File] output_read_chunks_1 = glob("reads.*.R1.fastq.gz")
        Array[File] output_read_chunks_2 = glob("reads.*.R2.fastq.gz")
    }
    runtime {
        cpu: in_cram_convert_cores
        memory: "50 GB"
        disks: "local-disk " + in_cram_convert_disk + " SSD"
        docker: "jmonlong/samtools-jm:release-1.19jm0.2.2"
        preemptible: 3
    }
}
