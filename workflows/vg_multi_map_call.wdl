version 1.0

### vg_multi_map_call.wdl ###
workflow vgMultiMapCall {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        description: "Core VG mapping and variant calling workflow for single sample datasets."
    }
    
    input {
        String MAPPER = "GIRAFFE"               # Set to 'MAP' to use the "VG MAP" algorithm, set to 'MPMAP' to use "VG MPMAP" algorithm, set to 'GIRAFFE' to use "VG GIRAFFE".
        Boolean SURJECT_MODE = true             # Set to 'true' to run pipeline using alignmed BAM files surjected from GAM. Set to 'false' to output graph aligned GAM files.
        Boolean DEEPVARIANT_MODE = false        # Set to 'true' to use the DeepVariant variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
        Boolean GVCF_MODE = false               # Set to 'true' to process and output gVCFs instead of VCFs.
        Boolean SNPEFF_ANNOTATION = true        # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        Boolean SV_CALLER_MODE = false          # Set to 'true' to run structural variant calling from graph aligned GAMs (SURJECT_MODE must be 'false' for this feature to be used)
        Boolean CLEANUP_FILES = true            # Set to 'false' to turn off intermediate file cleanup.
        Boolean GOOGLE_CLEANUP_MODE = false     # Set to 'true' to use google cloud compatible script for intermediate file cleanup. Set to 'false' to use local unix filesystem compatible script for intermediate file cleanup.
        File INPUT_READ_FILE_1                  # Input sample 1st read pair fastq.gz
        File INPUT_READ_FILE_2                  # Input sample 2nd read pair fastq.gz
        String SAMPLE_NAME                      # The sample name
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.28.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        String PCR_INDEL_MODEL = "CONSERVATIVE" # PCR indel model used in GATK Haplotypecaller (NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE)
        Int READS_PER_CHUNK = 20000000          # Number of reads contained in each mapping chunk (20000000 for wgs).
        Int CHUNK_BASES = 50000000              # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                      # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                    # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                            # Path to .xg index file
        File? GCSA_FILE                         # Path to .gcsa index file
        File? GCSA_LCP_FILE                     # Path to .gcsa.lcp index file
        File? GBWT_FILE                         # (OPTIONAL) Path to .gbwt index file
        File? GGBWT_FILE                        # (OPTIONAL) Path to .gg index file
        File? DIST_FILE                         # (OPTIONAL) Path to .dist index file
        File? MIN_FILE                          # (OPTIONAL) Path to .min index file
        File? SNARLS_FILE                       # (OPTIONAL) Path to .snarls index file
        File REF_FILE                           # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                     # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                      # Path to .dict file of the REF_FILE fasta reference
        File? SNPEFF_DATABASE                   # Path to snpeff database .zip file for snpEff annotation functionality.
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 32
        Int MAP_DISK = 10
        Int MAP_MEM = 60
        Int MERGE_GAM_CORES = 56
        Int MERGE_GAM_DISK = 100
        Int MERGE_GAM_MEM = 40
        Int MERGE_GAM_TIME = 2400
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 100
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 80
    }
    
    # Split input reads into chunks for parallelized mapping
    call splitReads as firstReadPair {
        input:
            in_read_file=INPUT_READ_FILE_1,
            in_pair_id="1",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }
    call splitReads as secondReadPair {
        input:
            in_read_file=INPUT_READ_FILE_2,
            in_pair_id="2",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }
    
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    if (!defined(PATH_LIST_FILE)) {
        call extractPathNames {
            input:
                in_xg_file=XG_FILE,
                in_vg_container=VG_CONTAINER
        }
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractPathNames.output_path_list])
    
    ################################################################
    # Distribute vg mapping opperation over each chunked read pair #
    ################################################################
    Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
    scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
        if (MAPPER == "MPMAP") {
            call runVGMPMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_snarls_file=SNARLS_FILE,
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        } 
        if (MAPPER == "MAP") {
            call runVGMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        if (MAPPER == "GIRAFFE") {
            call runVGGIRAFFE {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_ggbwt_file=GGBWT_FILE,
                    in_dist_file=DIST_FILE,
                    in_min_file=MIN_FILE,
                    in_ref_dict=REF_DICT_FILE,
                    in_sample_name=SAMPLE_NAME,
                    surject_output=SURJECT_MODE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        
        # Surject GAM alignment files to BAM if SURJECT_MODE set to true
        File vg_map_algorithm_chunk_gam_output = select_first([runVGMAP.chunk_gam_file, runVGMPMAP.chunk_gam_file, runVGGIRAFFE.chunk_gam_file])
        # Cleanup input reads after use
        if (CLEANUP_FILES) {
            if (GOOGLE_CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpVGMapperInputsGoogle {
                    input:
                        previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                        current_task_output = vg_map_algorithm_chunk_gam_output
                }
            }
            if (!GOOGLE_CLEANUP_MODE) {
                call cleanUpUnixFilesystem as cleanUpVGMapperInputsUnix {
                    input:
                        previous_task_outputs = [read_pair_chunk_files.left, read_pair_chunk_files.right],
                        current_task_output = vg_map_algorithm_chunk_gam_output
                }
            }
        }
        if (SURJECT_MODE) {
            call sortMDTagBAMFile {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_chunk_file=vg_map_algorithm_chunk_gam_output,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM,
            }
            # Cleanup intermediate surject files after use
            if (CLEANUP_FILES) {
                if (GOOGLE_CLEANUP_MODE) {
                    call cleanUpGoogleFilestore as cleanUpVGSurjectInputsGoogle {
                        input:
                            previous_task_outputs = [vg_map_algorithm_chunk_gam_output],
                            current_task_output = sortMDTagBAMFile.mark_dupped_reordered_bam
                    }
                }
                if (!GOOGLE_CLEANUP_MODE) {
                    call cleanUpUnixFilesystem as cleanUpVGSurjectInputsUnix {
                        input:
                            previous_task_outputs = [vg_map_algorithm_chunk_gam_output],
                            current_task_output = sortMDTagBAMFile.mark_dupped_reordered_bam
                    }
                }
            }
        }
    }
    
    ##############################################
    # Run the linear alignment calling procedure #
    ##############################################
    if (SURJECT_MODE) {
        # Merge chunked alignments from surjected GAM files
        Array[File?] alignment_chunk_bam_files_maybes = sortMDTagBAMFile.mark_dupped_reordered_bam
        Array[File] alignment_chunk_bam_files_valid = select_all(alignment_chunk_bam_files_maybes)
        call mergeAlignmentBAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=alignment_chunk_bam_files_valid,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        File merged_bam_file_output = mergeAlignmentBAMChunks.merged_bam_file
        File merged_bam_file_index_output = mergeAlignmentBAMChunks.merged_bam_file_index
        if (CLEANUP_FILES) {
            if (GOOGLE_CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpAlignmentBAMChunksGoogle {
                    input:
                        previous_task_outputs = alignment_chunk_bam_files_valid,
                        current_task_output = merged_bam_file_output
                }
            }
            if (!GOOGLE_CLEANUP_MODE) {
                call cleanUpUnixFilesystem as cleanUpAlignmentBAMChunksUnix {
                    input:
                        previous_task_outputs = alignment_chunk_bam_files_valid,
                        current_task_output = merged_bam_file_output
                }
            }
        }
         
        # Run VCF variant calling procedure
        # Split merged alignment by contigs list
        call splitBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_bam_file=merged_bam_file_output,
                in_merged_bam_file_index=merged_bam_file_index_output,
                in_path_list_file=pipeline_path_list_file,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        # Run distributed GATK linear variant calling
        scatter (gatk_caller_input_files in splitBAMbyPath.bams_and_indexes_by_contig) { 
            if (!GVCF_MODE) {
                # Run regular VCF genotypers
                if (!DEEPVARIANT_MODE) {
                    call runGATKHaplotypeCaller {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_pcr_indel_model=PCR_INDEL_MODEL,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call cleanUpGoogleFilestore as cleanUpGATKCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call cleanUpUnixFilesystem as cleanUpGATKCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCaller.genotyped_vcf
                            }
                        }
                    }
                }
                if (DEEPVARIANT_MODE) {
                    call runDeepVariantCaller {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call cleanUpGoogleFilestore as cleanUpDeepVariantCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCaller.genotyped_vcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call cleanUpUnixFilesystem as cleanUpDeepVariantCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCaller.genotyped_vcf
                            }
                        }
                    }
                }
            }    
            # Run GVCF variant calling procedure
            if (GVCF_MODE) {
                # Run GVCF genotypers
                if (!DEEPVARIANT_MODE) {
                    call runGATKHaplotypeCallerGVCF {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_pcr_indel_model=PCR_INDEL_MODEL,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call cleanUpGoogleFilestore as cleanUpGATKGVCFCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCallerGVCF.genotyped_gvcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call cleanUpUnixFilesystem as cleanUpGATKGVCFCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runGATKHaplotypeCallerGVCF.genotyped_gvcf
                            }
                        }
                    }
                }
                if (DEEPVARIANT_MODE) {
                    call runDeepVariantCallerGVCF {
                        input:
                            in_sample_name=SAMPLE_NAME,
                            in_bam_file=gatk_caller_input_files.left,
                            in_bam_file_index=gatk_caller_input_files.right,
                            in_reference_file=REF_FILE,
                            in_reference_index_file=REF_INDEX_FILE,
                            in_reference_dict_file=REF_DICT_FILE,
                            in_vgcall_cores=VGCALL_CORES,
                            in_vgcall_disk=VGCALL_DISK,
                            in_vgcall_mem=VGCALL_MEM
                    }
                    # Cleanup intermediate variant calling files after use
                    if (CLEANUP_FILES) {
                        if (GOOGLE_CLEANUP_MODE) {
                            call cleanUpGoogleFilestore as cleanUpDeepVariantGVCFCallerInputsGoogle {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCallerGVCF.genotyped_gvcf
                            }
                        }
                        if (!GOOGLE_CLEANUP_MODE) {
                            call cleanUpUnixFilesystem as cleanUpDeepVariantGVCFCallerInputsUnix {
                                input:
                                    previous_task_outputs = [gatk_caller_input_files.left, gatk_caller_input_files.right],
                                    current_task_output = runDeepVariantCallerGVCF.genotyped_gvcf
                            }
                        }
                    }
                }
            }
            File final_contig_vcf = select_first([runGATKHaplotypeCaller.genotyped_vcf, runDeepVariantCaller.genotyped_vcf, runGATKHaplotypeCallerGVCF.genotyped_gvcf, runDeepVariantCallerGVCF.genotyped_gvcf])
        }
        # Merge linear-based called VCFs
        Array[File?] final_contig_vcf_output = final_contig_vcf
        Array[File] final_contig_vcf_output_list = select_all(final_contig_vcf_output)
        call concatClippedVCFChunks as concatLinearVCFChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_clipped_vcf_chunk_files=final_contig_vcf_output_list,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call bgzipMergedVCF as bgzipLinearCalledVCF { 
            input: 
                in_sample_name=SAMPLE_NAME, 
                in_merged_vcf_file=concatLinearVCFChunks.output_merged_vcf,
                in_vg_container=VG_CONTAINER, 
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        
        File final_vcf_output = bgzipLinearCalledVCF.output_merged_vcf
    }
    
    ################################
    # Run the VG calling procedure #
    ################################
    if (!SURJECT_MODE) {
        # Merge chunked graph alignments
        call mergeAlignmentGAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_vg_container=VG_CONTAINER,
                in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
                in_merge_gam_cores=MERGE_GAM_CORES,
                in_merge_gam_disk=MERGE_GAM_DISK,
                in_merge_gam_mem=MERGE_GAM_MEM,
                in_merge_gam_time=MERGE_GAM_TIME
        }
        # Cleanup gam chunk files after use
        if (CLEANUP_FILES) {
            if (GOOGLE_CLEANUP_MODE) {
                call cleanUpGoogleFilestore as cleanUpGAMChunksGoogle {
                    input:
                        previous_task_outputs = vg_map_algorithm_chunk_gam_output,
                        current_task_output = mergeAlignmentGAMChunks.merged_sorted_gam_file
                }
            }
            if (!GOOGLE_CLEANUP_MODE) {
                call cleanUpUnixFilesystem as cleanUpGAMChunksUnix {
                    input:
                        previous_task_outputs = vg_map_algorithm_chunk_gam_output,
                        current_task_output = mergeAlignmentGAMChunks.merged_sorted_gam_file
                }
            }
        }
        # Set chunk context amount depending on structural variant or default variant calling mode
        Int default_or_sv_chunk_context = if SV_CALLER_MODE then 2500 else 50
        # Chunk GAM alignments by contigs list
        call chunkAlignmentsByPathNames {
            input:
            in_sample_name=SAMPLE_NAME,
            in_merged_sorted_gam=mergeAlignmentGAMChunks.merged_sorted_gam_file,
            in_merged_sorted_gam_gai=mergeAlignmentGAMChunks.merged_sorted_gam_gai_file,
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_chunk_context=default_or_sv_chunk_context,
            in_chunk_bases=CHUNK_BASES,
            in_overlap=OVERLAP,
            in_vg_container=VG_CONTAINER,
            in_chunk_gam_cores=CHUNK_GAM_CORES,
            in_chunk_gam_disk=CHUNK_GAM_DISK,
            in_chunk_gam_mem=CHUNK_GAM_MEM
        }
        Array[Pair[File,File]] vg_caller_input_files_list = zip(chunkAlignmentsByPathNames.output_vg_chunks, chunkAlignmentsByPathNames.output_gam_chunks) 
        # Run distributed graph-based variant calling
        scatter (vg_caller_input_files in vg_caller_input_files_list) { 
            call runVGCaller { 
                input: 
                    in_sample_name=SAMPLE_NAME, 
                    in_chunk_bed_file=chunkAlignmentsByPathNames.output_bed_chunk_file, 
                    in_vg_file=vg_caller_input_files.left, 
                    in_gam_file=vg_caller_input_files.right, 
                    in_chunk_bases=CHUNK_BASES, 
                    in_overlap=OVERLAP, 
                    in_vg_container=VG_CONTAINER,
                    in_vgcall_cores=VGCALL_CORES,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM,
                    in_sv_mode=SV_CALLER_MODE
            } 
            call runVCFClipper { 
                input: 
                    in_chunk_vcf=runVGCaller.output_vcf, 
                    in_chunk_vcf_index=runVGCaller.output_vcf_index, 
                    in_chunk_clip_string=runVGCaller.clip_string,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM
            } 
            # Cleanup vg call input files after use
            if (CLEANUP_FILES) {
                if (GOOGLE_CLEANUP_MODE) {
                    call cleanUpGoogleFilestore as cleanUpVGCallInputsGoogle {
                        input:
                            previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                            current_task_output = runVCFClipper.output_clipped_vcf
                    }
                }
                if (!GOOGLE_CLEANUP_MODE) {
                    call cleanUpUnixFilesystem as cleanUpVGCallInputsUnix {
                        input:
                            previous_task_outputs = [vg_caller_input_files.left, vg_caller_input_files.right, runVGCaller.output_vcf],
                            current_task_output = runVCFClipper.output_clipped_vcf
                    }
                }
            }
        } 
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksVGCall { 
            input: 
                in_sample_name=SAMPLE_NAME, 
                in_clipped_vcf_chunk_files=runVCFClipper.output_clipped_vcf,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        } 
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_vcf_file=concatVCFChunksVGCall.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    
    # Extract either the linear-based or graph-based VCF
    File variantcaller_vcf_output = select_first([bgzipVGCalledVCF.output_merged_vcf, final_vcf_output])
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION && defined(SNPEFF_DATABASE)) {
        call normalizeVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bgzip_vcf_file=variantcaller_vcf_output,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call snpEffAnnotateVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_normalized_vcf_file=normalizeVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    output {
        File output_vcf = select_first([snpEffAnnotateVCF.output_snpeff_annotated_vcf, variantcaller_vcf_output])
        File? output_bam = merged_bam_file_output
        File? output_bam_index = merged_bam_file_index_output
        File? output_gam = mergeAlignmentGAMChunks.merged_sorted_gam_file
        File? output_gam_index = mergeAlignmentGAMChunks.merged_sorted_gam_gai_file
    }
}

########################
### TASK DEFINITIONS ###
########################

# Tasks for intermediate file cleanup
task cleanUpUnixFilesystem {
    input {
        Array[String] previous_task_outputs
        String current_task_output
    }
    command <<<
        set -eux -o pipefail
        cat ~{write_lines(previous_task_outputs)} | sed 's/.*\(\/cromwell-executions\)/\1/g' | xargs -I {} ls -li {} | cut -f 1 -d ' ' | xargs -I {} find ../../../ -xdev -inum {} | xargs -I {} rm -v {}
    >>>
    runtime {
        docker: "ubuntu@sha256:2695d3e10e69cc500a16eae6d6629c803c43ab075fa5ce60813a0fc49c47e859"
        continueOnReturnCode: true
    }
}

task cleanUpGoogleFilestore {
    input {
        Array[String] previous_task_outputs
        String current_task_output
    }
    command {
        set -eux -o pipefail
        gsutil rm -I < ${write_lines(previous_task_outputs)}
    }
    runtime {
        docker: "google/cloud-sdk@sha256:4ef6b0e969fa96f10acfd893644d100469e979f4384e5e70f58be5cb80593a8a"
        continueOnReturnCode: true
    }
}

task splitReads {
    input {
        File in_read_file
        String in_pair_id
        String in_vg_container
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk
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

        CHUNK_LINES=$(( ~{in_reads_per_chunk} * 4 ))
        gzip -cd ~{in_read_file} | split -l $CHUNK_LINES --filter='pigz -p ~{in_split_read_cores} > ${FILE}.fq.gz' - "fq_chunk_~{in_pair_id}.part."
    >>>
    output {
        Array[File] output_read_chunks = glob("fq_chunk_~{in_pair_id}.part.*")
    }
    runtime {
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: in_vg_container
    }
}

task extractPathNames {
    input {
        File in_xg_file
        String in_vg_container
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt
    }
    output {
        File output_path_list = "path_list.txt"
    }
    runtime {
        memory: "20 GB"
        disks: "local-disk 50 SSD"
        docker: in_vg_container
    }
}

task runVGMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File? in_gcsa_file
        File? in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }
    
    Boolean gbwt_options = defined(in_gbwt_file)
    
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

        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        GBWT_OPTION_STRING=""
        if [ ~{gbwt_options} == true ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file}"
        fi
        ln -s ~{in_gcsa_file} input_gcsa_file.gcsa
        ln -s ~{in_gcsa_lcp_file} input_gcsa_file.gcsa.lcp
        if [ ~{surject_output} == false ]; then
            vg map \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else    
           vg map \
              --ref-paths ~{in_ref_dict} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --surject-to bam \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task runVGMPMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File? in_gcsa_file
        File? in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_snarls_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }
    
    Boolean gbwt_options = defined(in_gbwt_file)
    Boolean snarl_options = defined(in_snarls_file)
    
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
        
        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        GBWT_OPTION_STRING=""
        if [ ~{gbwt_options} == true ] && [ ~{snarl_options} == false ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file} --recombination-penalty 5.0"
        elif [ ~{gbwt_options} == true ] && [ ~{snarl_options} == true ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file} -s ~{in_snarls_file} --recombination-penalty 5.0"
        fi
        ln -s ~{in_gcsa_file} input_gcsa_file.gcsa
        ln -s ~{in_gcsa_lcp_file} input_gcsa_file.gcsa.lcp
        if [ ~{surject_output} == false ]; then
            vg mpmap \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -S \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else    
            vg mpmap \
              --ref-paths ~{in_ref_dict} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --output-fmt BAM \
              -S \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -g input_gcsa_file.gcsa \
              ${GBWT_OPTION_STRING} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}
task runVGGIRAFFE {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File? in_gbwt_file
        File? in_ggbwt_file
        File? in_dist_file
        File? in_min_file
        File? in_ref_dict
        String in_vg_container
        String in_sample_name
        Boolean surject_output
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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
        
        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        if [ ~{surject_output} == false ]; then
            vg giraffe \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -H ~{in_gbwt_file} \
              -g ~{in_ggbwt_file} \
              -d ~{in_dist_file} \
              -m ~{in_min_file} \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
        else    
            vg giraffe \
              --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
              --sample "~{in_sample_name}" \
              --output-format BAM \
              --ref-paths ~{in_ref_dict} \
              -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
              -x ~{in_xg_file} \
              -H ~{in_gbwt_file} \
              -g ~{in_ggbwt_file} \
              -d ~{in_dist_file} \
              -m ~{in_min_file} \
              -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
        fi
    >>>
    output {
        File chunk_gam_file = glob("*am")[0]
    }
    runtime {
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task sortMDTagBAMFile {
    input {
        String in_sample_name
        File in_bam_chunk_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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
        samtools sort \
          --threads ~{in_map_cores} \
          ~{in_bam_chunk_file} \
          -O BAM \
        | samtools calmd \
          -b \
          - \
          ~{in_reference_file} \
          > ~{in_sample_name}_positionsorted.mdtag.bam \
        && samtools index \
          ~{in_sample_name}_positionsorted.mdtag.bam
        java -Xmx~{in_map_mem}g -XX:ParallelGCThreads=~{in_map_cores} -jar /usr/picard/picard.jar MarkDuplicates \
          PROGRAM_RECORD_ID=null \
          VALIDATION_STRINGENCY=LENIENT \
          I=~{in_sample_name}_positionsorted.mdtag.bam \
          O=~{in_sample_name}.mdtag.dupmarked.bam \
          M=marked_dup_metrics.txt 2> mark_dup_stderr.txt \
        && rm -f ~{in_sample_name}_positionsorted.mdtag.bam ~{in_sample_name}_positionsorted.mdtag.bam.bai \
    >>>
    output {
        File mark_dupped_reordered_bam = "~{in_sample_name}.mdtag.dupmarked.bam"
    }
    runtime {
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
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
        samtools merge \
          -f -p -c --threads ~{in_map_cores} \
          ~{in_sample_name}_merged.positionsorted.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.positionsorted.bam
    >>>
    output {
        File merged_bam_file = "~{in_sample_name}_merged.positionsorted.bam"
        File merged_bam_file_index = "~{in_sample_name}_merged.positionsorted.bam.bai"
    }
    runtime {
        time: 240
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }
    
    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai
        
        while IFS=$'\t' read -ra path_list_line; do
            path_name="${path_list_line[0]}"
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${path_name} \
              -o ~{in_sample_name}.${path_name}.bam \
            && samtools index \
              ~{in_sample_name}.${path_name}.bam
        done < ~{in_path_list_file}
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai") 
        Array[Pair[File, File]] bams_and_indexes_by_contig = zip(bam_contig_files, bam_contig_files_index)
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKHaplotypeCaller {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        String in_pcr_indel_model
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        
        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        gatk HaplotypeCaller \
          --native-pair-hmm-threads ~{in_vgcall_cores} \
          --pcr-indel-model ~{in_pcr_indel_model} \
          -L ${CONTIG_ID} \
          --reference ~{in_reference_file} \
          --input input_bam_file.bam \
          --output ~{in_sample_name}.vcf \
        && bgzip ~{in_sample_name}.vcf
    >>>
    output {
        File genotyped_vcf = "~{in_sample_name}.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "broadinstitute/gatk@sha256:cc8981d0527e716775645b04a7f59e96a52ad59a7ae9788ddc47902384bf35aa"
    }
}

task runGATKHaplotypeCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        String in_pcr_indel_model
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        
        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        gatk HaplotypeCaller \
          --native-pair-hmm-threads ~{in_vgcall_cores} \
          -ERC GVCF \
          -L ${CONTIG_ID} \
          --pcr-indel-model ~{in_pcr_indel_model} \
          --reference ~{in_reference_file} \
          --input input_bam_file.bam \
          --output ~{in_sample_name}.rawLikelihoods.gvcf \
        && bgzip ~{in_sample_name}.rawLikelihoods.gvcf
    >>>
    output {
        File genotyped_gvcf = "~{in_sample_name}.rawLikelihoods.gvcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "broadinstitute/gatk@sha256:cc8981d0527e716775645b04a7f59e96a52ad59a7ae9788ddc47902384bf35aa"
    }
}

task runDeepVariantCaller {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        
        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        
        /opt/deepvariant/bin/run_deepvariant \
        --num_shards ~{in_vgcall_cores} \
        --model_type WGS \
        --regions ${CONTIG_ID} \
        --ref ~{in_reference_file} \
        --reads input_bam_file.bam \
        --output_vcf ~{in_sample_name}.vcf \
        && bgzip ~{in_sample_name}.vcf
    >>>
    output {
        File genotyped_vcf = "~{in_sample_name}.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:1.0.0"
    }
}

task runDeepVariantCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        
        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | awk -F'.' '{print $(NF-1)}'))
        
        /opt/deepvariant/bin/run_deepvariant \
        --num_shards ~{in_vgcall_cores} \
        --model_type WGS \
        --regions ${CONTIG_ID} \
        --ref ~{in_reference_file} \
        --reads input_bam_file.bam \
        --output_gvcf ~{in_sample_name}.rawLikelihoods.gvcf \
        && bgzip ~{in_sample_name}.rawLikelihoods.gvcf
    >>>
    output {
        File genotyped_gvcf = "~{in_sample_name}.rawLikelihoods.gvcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:1.0.0"
    }
}

task mergeAlignmentGAMChunks {
    input {
        String in_sample_name
        String in_vg_container
        Array[File] in_alignment_gam_chunk_files
        Int in_merge_gam_cores
        Int in_merge_gam_disk
        Int in_merge_gam_mem
        Int in_merge_gam_time
    }
    
    command {
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
        
        cat ${sep=" " in_alignment_gam_chunk_files} > ${in_sample_name}_merged.gam \
        && vg gamsort \
            ${in_sample_name}_merged.gam \
            -i ${in_sample_name}_merged.sorted.gam.gai \
            -t ${in_merge_gam_cores} > ${in_sample_name}_merged.sorted.gam
    }
    output {
        File merged_sorted_gam_file = "${in_sample_name}_merged.sorted.gam"
        File merged_sorted_gam_gai_file = "${in_sample_name}_merged.sorted.gam.gai"
    }
    runtime {
        memory: in_merge_gam_mem + " GB"
        cpu: in_merge_gam_cores
        disks: "local-disk " + in_merge_gam_disk  + " SSD"
        time: in_merge_gam_time
        docker: in_vg_container
    }
}

task chunkAlignmentsByPathNames {
    input {
        String in_sample_name
        File in_merged_sorted_gam
        File in_merged_sorted_gam_gai
        File in_xg_file
        Int in_chunk_context
        File in_path_list_file
        Int in_chunk_bases
        Int in_overlap
        String in_vg_container
        Int in_chunk_gam_cores
        Int in_chunk_gam_disk
        Int in_chunk_gam_mem
    }
    
    command { 
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
    
        vg chunk \
            -x ${in_xg_file} \
            -a ${in_merged_sorted_gam} \
            -c ${in_chunk_context} \
            -P ${in_path_list_file} \
            -g \
            -s ${in_chunk_bases} -o ${in_overlap} \
            -b call_chunk \
            -t ${in_chunk_gam_cores} \
            -E output_bed_chunks.bed -f
    }
    output {
        Array[File] output_vg_chunks = glob("call_chunk*.vg")
        Array[File] output_gam_chunks = glob("call_chunk*.gam")
        File output_bed_chunk_file = "output_bed_chunks.bed"
    }   
    runtime {
        time: 900
        memory: in_chunk_gam_mem + " GB"
        cpu: in_chunk_gam_cores
        disks: "local-disk " + in_chunk_gam_disk + " SSD"
        docker: in_vg_container
    }   
}

task runVGCaller {
    input {
        String in_sample_name
        File in_chunk_bed_file
        File in_vg_file
        File in_gam_file
        Int in_chunk_bases
        Int in_overlap
        String in_vg_container
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
        Boolean in_sv_mode
    }

    String chunk_tag = basename(in_vg_file, ".vg")
    command <<<
        set -eux -o pipefail

        VG_INDEX_XG_COMMAND=""
        VG_FILTER_COMMAND=""
        VG_AUGMENT_SV_OPTIONS=""
        VG_CALL_SV_OPTIONS=""
        if [ ~{in_sv_mode} == false ]; then
            VG_INDEX_XG_COMMAND="vg index ~{in_vg_file} -x ~{chunk_tag}.xg -t ~{in_vgcall_cores} && \\"
            VG_FILTER_COMMAND="vg filter ~{in_gam_file} -t 1 -r 0.9 -fu -s 1000 -m 1 -q 15 -D 999 -x ~{chunk_tag}.xg > ~{chunk_tag}.filtered.gam && \\"
            VG_AUGMENT_SV_OPTIONS="~{chunk_tag}.filtered_gam"
        else
            VG_AUGMENT_SV_OPTIONS="~{in_gam_file} --recall"
            VG_CALL_SV_OPTIONS="-u -n 0 -e 1000 -G 3"
        fi

        PATH_NAME="$(echo ~{chunk_tag} | cut -f 4 -d '_')"
        OFFSET=""
        BED_CHUNK_LINE_RECORD=""
        FIRST_Q=""
        LAST_Q=""
        CHUNK_INDEX_COUNTER="0"
        CHUNK_TAG_INDEX="0"
        CHR_LENGTH=""
        while IFS=$'\t' read -ra bed_chunk_line; do
            bed_line_chunk_tag="$(echo ${bed_chunk_line[3]} | cut -f 1 -d '.')"
            if [ "~{chunk_tag}" = "${bed_line_chunk_tag}" ]; then
                OFFSET="${bed_chunk_line[1]}"
                BED_CHUNK_LINE_RECORD="${bed_chunk_line[@]//$'\t'/ }"
                CHUNK_TAG_INDEX="${CHUNK_INDEX_COUNTER}"
            fi
            bed_line_path_name="$(echo ${bed_line_chunk_tag} | cut -f 4 -d '_')"
            if [ "${PATH_NAME}" = "${bed_line_path_name}" ]; then
                if [ "${FIRST_Q}" = "" ]; then
                    FIRST_Q="${bed_chunk_line[@]//$'\t'/ }"
                fi
                LAST_Q="${bed_chunk_line[@]//$'\t'/ }"
                CHUNK_INDEX_COUNTER="$(( ${CHUNK_INDEX_COUNTER} + 1 ))"
                CHR_LENGTH="${bed_chunk_line[2]}"
            fi
        done < ~{in_chunk_bed_file}
        
        ${VG_INDEX_XG_COMMAND}
        ${VG_FILTER_COMMAND}
        vg augment \
            ~{in_vg_file} \
            ${VG_AUGMENT_SV_OPTIONS} \
            -t ~{in_vgcall_cores} \
            -a pileup \
            -Z ~{chunk_tag}.trans \
            -S ~{chunk_tag}.support > ~{chunk_tag}.aug.vg && \
        vg call \
            ~{chunk_tag}.aug.vg \
            -t ~{in_vgcall_cores} \
            -S ~{in_sample_name} \
            -z ~{chunk_tag}.trans \
            -s ~{chunk_tag}.support \
            -b ~{in_vg_file} \
            -r ${PATH_NAME} \
            -c ${PATH_NAME} \
            -l ${CHR_LENGTH} \
            ${VG_CALL_SV_OPTIONS} \
            -o ${OFFSET} > ~{chunk_tag}.vcf && \
        head -10000 ~{chunk_tag}.vcf | grep "^#" >> ~{chunk_tag}.sorted.vcf && \
        if [ "$(cat ~{chunk_tag}.vcf | grep -v '^#')" ]; then
            cat ~{chunk_tag}.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ~{chunk_tag}.sorted.vcf
        fi
        bgzip ~{chunk_tag}.sorted.vcf && \
        tabix -f -p vcf ~{chunk_tag}.sorted.vcf.gz

        # Compile clipping string
        clipped_chunk_offset="$(( (${CHUNK_TAG_INDEX} * ~{in_chunk_bases}) - (${CHUNK_TAG_INDEX} * ~{in_overlap}) ))"
        if [ "${BED_CHUNK_LINE_RECORD}" = "${FIRST_Q}" ]; then
          START="$(( ${clipped_chunk_offset} + 1 ))"
        else
          START="$(( ${clipped_chunk_offset} + 1 + ~{in_overlap}/2 ))"
        fi
        if [ "${BED_CHUNK_LINE_RECORD}" = "${LAST_Q}" ]; then
          END="$(( ${clipped_chunk_offset} + ~{in_chunk_bases}))"
        else
          END="$(( ${clipped_chunk_offset} + ~{in_chunk_bases} - (~{in_overlap}/2) ))"
        fi
        echo "${PATH_NAME}:${START}-${END}" > ~{chunk_tag}_clip_string.txt
    >>>
    output {
        File output_vcf = "~{chunk_tag}.sorted.vcf.gz"
        File output_vcf_index = "~{chunk_tag}.sorted.vcf.gz.tbi"
        String clip_string = read_lines("~{chunk_tag}_clip_string.txt")[0]
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: in_vg_container
    }
}

task runVCFClipper {
    input {
        File in_chunk_vcf
        File in_chunk_vcf_index
        String in_chunk_clip_string
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    String chunk_tag = basename(in_chunk_vcf, ".sorted.vcf.gz")
    command {
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

        bcftools view -O b ${in_chunk_vcf} -t ${in_chunk_clip_string} > ${chunk_tag}.clipped.vcf.gz
    }
    output {
        File output_clipped_vcf = "${chunk_tag}.clipped.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    command {
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

        for vcf_file in ${sep=" " in_clipped_vcf_chunk_files} ; do
            bcftools index "$vcf_file"
        done
        bcftools concat -a ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort - > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        time: 60
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    # TODO:
    #   If GVCF in in_merged_vcf_file then output_vcf_extension="gvcf" else output_vcf_extension="vcf"
    command {
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

        bgzip -c ${in_merged_vcf_file} > ${in_sample_name}_merged.vcf.gz && \
        tabix -f -p vcf ${in_sample_name}_merged.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}_merged.vcf.gz.tbi"
    }
    runtime {
        time: 30
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: in_vg_container
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        bcftools norm -m-both --threads ~{in_vgcall_cores} -o ~{in_sample_name}.unrolled.vcf ~{in_bgzip_vcf_file}
    >>>
    output {
        File output_normalized_vcf = "~{in_sample_name}.unrolled.vcf"
    }
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File? in_snpeff_database
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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
        
        unzip ~{in_snpeff_database}
        snpEff -Xmx40g -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data GRCh37.75 ~{in_normalized_vcf_file} > ~{in_sample_name}.snpeff.unrolled.vcf
    >>>
    output {
        File output_snpeff_annotated_vcf = "~{in_sample_name}.snpeff.unrolled.vcf"
    }
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:4.3.1t--2"
    }
}

