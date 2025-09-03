version 1.0

### giraffe_and_deeptrio.wdl ###
## Author: Charles Markello
## Description: Core VG Giraffe mapping and DeepTrio calling workflow for maternal-paternal-child sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

import "./giraffe.wdl" as GiraffeWorkflow
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/bioinfo_utils.wdl" as utils

workflow vgGiraffeDeeptrio {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        Updated: "Parsa Eskandar"
        description: "Core VG Giraffe mapping and DeepTrio calling workflow for maternal-paternal-child sample datasets. It takes as inputs reads in FASTQ and graphs containing the population-based haplotypes to genotype. The graphs files required include the XG, GCSA, GBWT, graph GBWT, Distance and Minimizer indexes. It outputs a VCF file and BAM file for the child along with optional RTG and hap.py vcf evaluation if the user provides benchmark truth-set VCFs."
    }
    input {
        File CHILD_INPUT_READ_FILE_1                    # Input child sample 1st read pair fastq.gz
        File CHILD_INPUT_READ_FILE_2                    # Input child sample 2nd read pair fastq.gz
        File MATERNAL_INPUT_READ_FILE_1                 # Input maternal sample 1st read pair fastq.gz
        File MATERNAL_INPUT_READ_FILE_2                 # Input maternal sample 2nd read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_1                 # Input paternal sample 1st read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_2                 # Input paternal sample 2nd read pair fastq.gz
        String SAMPLE_NAME                              # The child sample name
        String MATERNAL_NAME                            # The maternal sample name
        String PATERNAL_NAME                            # The paternal sample name
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        String? GIRAFFE_OPTIONS                         # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS.
        String REFERENCE_PREFIX = ""                    # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File GBZ_FILE                                   # Path to .gbz graph file
        File? DIST_FILE                                 # Path to .dist index file
        File? MIN_FILE                                  # Path to .min index file
        File? ZIPCODES_FILE                             # (OPTIONAL) Zipcodes index file
        File? TRUTH_VCF                                 # Path to .vcf.gz to compare against
        File? TRUTH_VCF_INDEX                           # Path to Tabix index for TRUTH_VCF
        File? EVALUATION_REGIONS_BED                    # BED to restrict comparison against TRUTH_VCF to
        File? DV_MODEL_META                             # .meta file for a custom DeepVariant calling model
        File? DV_MODEL_INDEX                            # .index file for a custom DeepVariant calling model
        File? DV_MODEL_DATA                             # .data-00000-of-00001 file for a custom DeepVariant calling model
        File? CHILD_DT_MODEL_META                       # .meta file for a custom DeepVariant calling model
        File? CHILD_DT_MODEL_INDEX                      # .index file for a custom DeepVariant calling model
        File? CHILD_DT_MODEL_DATA                       # .data-00000-of-00001 file for a custom DeepVariant calling model
        File? PARENT_DT_MODEL_META                      # .meta file for a custom DeepVariant calling model
        File? PARENT_DT_MODEL_INDEX                     # .index file for a custom DeepVariant calling model
        File? PARENT_DT_MODEL_DATA                      # .data-00000-of-00001 file for a custom DeepVariant calling model
        Boolean LEFTALIGN_BAM = false                   # Whether or not to left-align reads in the BAM before DV
        Boolean REALIGN_INDELS = false                  # Whether or not to realign reads near indels
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        Boolean? DV_KEEP_LEGACY_AC                      # Should DV use the legacy allele counter behavior?
        Boolean? DV_NORM_READS                          # Should DV normalize reads?
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 16
        Int MAP_DISK = 10
        Int MAP_MEM = 50
        Int CALL_CORES = 8
        Int CALL_DISK = 40
        Int CALL_MEM = 50
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE
        String VG_DOCKER = "quay.io/vgteam/vg:v1.64.0"
        String DEEPTRIO_DOCKER = "google/deepvariant:deeptrio-1.9.0"
        Boolean MERGE_TRIO_GVCFS = false                 # Optionally merge trio gVCFs using GLnexus
    }

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from GBZ file if PATH_LIST_FILE input not provided
            call map.extractSubsetPathNames as extractSubsetPathNamesGBZ {
                input:
                    in_gbz_file=GBZ_FILE,
                    in_reference_prefix=REFERENCE_PREFIX,
                    in_extract_mem=MAP_MEM,
                    vg_docker=VG_DOCKER
            }
        }
    }
    if (defined(CONTIGS)) {
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNamesGBZ.output_path_list_file, written_path_names_file])
    
    # Generate reference FASTA from GBZ to match chosen paths
    if (!defined(REFERENCE_FILE)) {
        call map.extractReference as extractReferenceFromGBZ {
            input:
                in_gbz_file=GBZ_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
                in_extract_mem=MAP_MEM,
                vg_docker=VG_DOCKER
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReferenceFromGBZ.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference as buildRefIndex {
            input:
                in_reference_file=reference_file,
                in_index_disk=MAP_DISK,
                in_index_mem=MAP_MEM
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, buildRefIndex.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, buildRefIndex.reference_dict_file])

    #######################################################
    ############ Run mapping workflows on Trio ############
    #######################################################
    call GiraffeWorkflow.Giraffe as maternalMapWorkflow {
        input:
            INPUT_READ_FILE_1=MATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=MATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=MATERNAL_NAME,
            READS_PER_CHUNK=READS_PER_CHUNK,
            GBZ_FILE=GBZ_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            ZIPCODES_FILE=ZIPCODES_FILE,
            PATH_LIST_FILE=pipeline_path_list_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            VG_DOCKER=VG_DOCKER,
            OUTPUT_SINGLE_BAM=true,
            OUTPUT_CALLING_BAMS=true,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            MAP_CORES=MAP_CORES,
            MAP_MEM=MAP_MEM
    }
    call GiraffeWorkflow.Giraffe as paternalMapWorkflow {
        input:
            INPUT_READ_FILE_1=PATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=PATERNAL_NAME,
            READS_PER_CHUNK=READS_PER_CHUNK,
            GBZ_FILE=GBZ_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            ZIPCODES_FILE=ZIPCODES_FILE,
            PATH_LIST_FILE=pipeline_path_list_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            VG_DOCKER=VG_DOCKER,
            OUTPUT_SINGLE_BAM=true,
            OUTPUT_CALLING_BAMS=true,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            MAP_CORES=MAP_CORES,
            MAP_MEM=MAP_MEM
    }
    call GiraffeWorkflow.Giraffe as childMapWorkflow {
        input:
            INPUT_READ_FILE_1=CHILD_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=CHILD_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME,
            READS_PER_CHUNK=READS_PER_CHUNK,
            GBZ_FILE=GBZ_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            ZIPCODES_FILE=ZIPCODES_FILE,
            PATH_LIST_FILE=pipeline_path_list_file,
            REFERENCE_PREFIX=REFERENCE_PREFIX,
            VG_DOCKER=VG_DOCKER,
            OUTPUT_SINGLE_BAM=true,
            OUTPUT_CALLING_BAMS=true,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            MAP_CORES=MAP_CORES,
            MAP_MEM=MAP_MEM
    }


    ##################################################
    # Distribute deeptrio operation over each contig #
    ##################################################

    ##
    ## Call variants with DeepTrio in each contig
    ##
    
    # Run distributed DeepTrio linear variant calling for each chromosomal contig
    Array[File] maternal_calling_bams = select_first([maternalMapWorkflow.output_calling_bams, []])
    Array[File] maternal_calling_bam_indexes = select_first([maternalMapWorkflow.output_calling_bam_indexes, []])
    Array[File] paternal_calling_bams = select_first([paternalMapWorkflow.output_calling_bams, []])
    Array[File] paternal_calling_bam_indexes = select_first([paternalMapWorkflow.output_calling_bam_indexes, []])
    Array[File] child_calling_bams = select_first([childMapWorkflow.output_calling_bams, []])
    Array[File] child_calling_bam_indexes = select_first([childMapWorkflow.output_calling_bam_indexes, []])

    Array[Pair[File, File]] maternal_bams_and_indexes_by_contig = zip(maternal_calling_bams, maternal_calling_bam_indexes)
    Array[Pair[File, File]] paternal_bams_and_indexes_by_contig = zip(paternal_calling_bams, paternal_calling_bam_indexes)
    Array[Pair[File, File]] child_bams_and_indexes_by_contig = zip(child_calling_bams, child_calling_bam_indexes)
    Array[Pair[Pair[File,File],Pair[Pair[File,File],Pair[File,File]]]] trio_bam_index_by_contigs_pair = zip(child_bams_and_indexes_by_contig, zip(maternal_bams_and_indexes_by_contig, paternal_bams_and_indexes_by_contig))
    #              trio_bam_index_by_contigs_pair
    #           _________________|_________________
    #       ___/__                     ____________\_____________
    #      /      \             _____ /____                ______\____
    # child.bam child.bai      /            \             /           \       
    #                     maternal.bam maternal.bai paternal.bam paternal.bai        

    scatter (deeptrio_caller_input_files in trio_bam_index_by_contigs_pair) {
        File child_bam_file = deeptrio_caller_input_files.left.left
        File child_bam_file_index = deeptrio_caller_input_files.left.right
        File maternal_bam_file = deeptrio_caller_input_files.right.left.left
        File maternal_bam_file_index = deeptrio_caller_input_files.right.left.right
        File paternal_bam_file = deeptrio_caller_input_files.right.right.left
        File paternal_bam_file_index = deeptrio_caller_input_files.right.right.right

        ## DeepTrio calling
        String contig_name = sub(
            sub(
                sub(
                    sub(
                        sub(basename(child_bam_file), "\\.bam$", ""),
                        "\\.indel_realigned$", ""
                    ),
                    "\\.left_shifted$", ""
                ),
                "^" + SAMPLE_NAME + "\\.", ""
            ),
            "^\\.", ""
        )
        if ((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT")) {
            call runDeepVariantMakeExamples as callDeepVariantMakeExamplesChild {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=child_bam_file,
                    in_bam_file_index=child_bam_file_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_min_mapq=MIN_MAPQ,
                    in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                    in_norm_reads=DV_NORM_READS,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
            call runDeepVariantCallVariants as callDeepVariantCallVariantsChild {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_examples_file=callDeepVariantMakeExamplesChild.examples_file,
                    in_nonvariant_site_tf_file=callDeepVariantMakeExamplesChild.nonvariant_site_tf_file,
                    in_model_meta_file=DV_MODEL_META,
                    in_model_index_file=DV_MODEL_INDEX,
                    in_model_data_file=DV_MODEL_DATA,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
            call runDeepVariantMakeExamples as callDeepVariantMakeExamplesMaternal {
                input:
                    in_sample_name=MATERNAL_NAME,
                    in_bam_file=maternal_bam_file,
                    in_bam_file_index=maternal_bam_file_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_min_mapq=MIN_MAPQ,
                    in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                    in_norm_reads=DV_NORM_READS,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
            call runDeepVariantCallVariants as callDeepVariantCallVariantsMaternal {
                input:
                    in_sample_name=MATERNAL_NAME,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_examples_file=callDeepVariantMakeExamplesMaternal.examples_file,
                    in_nonvariant_site_tf_file=callDeepVariantMakeExamplesMaternal.nonvariant_site_tf_file,
                    in_model_meta_file=DV_MODEL_META,
                    in_model_index_file=DV_MODEL_INDEX,
                    in_model_data_file=DV_MODEL_DATA,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
            call runDeepVariantMakeExamples as callDeepVariantMakeExamplesPaternal {
                input:
                    in_sample_name=PATERNAL_NAME,
                    in_bam_file=paternal_bam_file,
                    in_bam_file_index=paternal_bam_file_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_min_mapq=MIN_MAPQ,
                    in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                    in_norm_reads=DV_NORM_READS,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
            call runDeepVariantCallVariants as callDeepVariantCallVariantsPaternal {
                input:
                    in_sample_name=PATERNAL_NAME,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_examples_file=callDeepVariantMakeExamplesPaternal.examples_file,
                    in_nonvariant_site_tf_file=callDeepVariantMakeExamplesPaternal.nonvariant_site_tf_file,
                    in_model_meta_file=DV_MODEL_META,
                    in_model_index_file=DV_MODEL_INDEX,
                    in_model_data_file=DV_MODEL_DATA,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM
            }
        }
        if (!((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT"))) {
            call runDeepTrioOneCommand as trioCall {
                input:
                    in_child_name=SAMPLE_NAME,
                    in_maternal_name=MATERNAL_NAME,
                    in_paternal_name=PATERNAL_NAME,
                    in_child_bam_file=child_bam_file,
                    in_child_bam_file_index=child_bam_file_index,
                    in_maternal_bam_file=maternal_bam_file,
                    in_maternal_bam_file_index=maternal_bam_file_index,
                    in_paternal_bam_file=paternal_bam_file,
                    in_paternal_bam_file_index=paternal_bam_file_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_region=contig_name,
                    in_call_cores=CALL_CORES,
                    in_call_disk=CALL_DISK,
                    in_call_mem=CALL_MEM,
                    in_deeptrio_docker=DEEPTRIO_DOCKER
            }
        }
    }
    Array[File] childDeepVarGVCF = select_all(callDeepVariantCallVariantsChild.output_gvcf_file)
    Array[File] childDeepTrioGVCF = select_all(trioCall.output_child_gvcf_file)
    Array[File] child_contig_gvcf_output_list = select_all(flatten([childDeepTrioGVCF, childDeepVarGVCF]))
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks as concatVCFChunksChild {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=child_contig_gvcf_output_list,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF as bgzipVGCalledChildVCF {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_vcf_file=concatVCFChunksChild.output_merged_vcf,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    Array[File] maDeepVarGVCF = select_all(callDeepVariantCallVariantsMaternal.output_gvcf_file)
    Array[File] maDeepTrioGVCF = select_all(trioCall.output_parent2_gvcf_file)
    Array[File] maternal_contig_gvcf_output_list = select_all(flatten([maDeepTrioGVCF, maDeepVarGVCF]))
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks as concatVCFChunksMaternal {
        input:
            in_sample_name=MATERNAL_NAME,
            in_clipped_vcf_chunk_files=maternal_contig_gvcf_output_list,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF as bgzipVGCalledMaternalVCF {
        input:
            in_sample_name=MATERNAL_NAME,
            in_merged_vcf_file=concatVCFChunksMaternal.output_merged_vcf,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    Array[File] paDeepVarGVCF = select_all(callDeepVariantCallVariantsPaternal.output_gvcf_file)
    Array[File] paDeepTrioGVCF = select_all(trioCall.output_parent1_gvcf_file)
    Array[File] paternal_contig_gvcf_output_list = select_all(flatten([paDeepTrioGVCF, paDeepVarGVCF]))
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks as concatVCFChunksPaternal {
        input:
            in_sample_name=PATERNAL_NAME,
            in_clipped_vcf_chunk_files=paternal_contig_gvcf_output_list,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF as bgzipVGCalledPaternalVCF {
        input:
            in_sample_name=PATERNAL_NAME,
            in_merged_vcf_file=concatVCFChunksPaternal.output_merged_vcf,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    if (MERGE_TRIO_GVCFS) {
        call glnexusMergeTrioBCF {
            input:
                in_child_gvcf_file=bgzipVGCalledChildVCF.output_merged_vcf,
                in_parent1_gvcf_file=bgzipVGCalledPaternalVCF.output_merged_vcf,
                in_parent2_gvcf_file=bgzipVGCalledMaternalVCF.output_merged_vcf,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
        call convertBCFtoVCF as convertMergedTrioBCFtoVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bcf_file=glnexusMergeTrioBCF.output_bcf_file,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM,
                in_call_cores=CALL_CORES
        }
    }
    if (defined(TRUTH_VCF) && defined(TRUTH_VCF_INDEX)) {
    
        # To evaluate the VCF we need a template of the reference
        call buildReferenceTemplate {
            input:
                in_reference_file=reference_file
        }
        
        # Direct vcfeval comparison makes an archive with FP and FN VCFs
        call compareCalls {
            input:
                in_sample_vcf_file=bgzipVGCalledChildVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipVGCalledChildVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_template_archive=buildReferenceTemplate.output_template_archive,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
        
        # Hap.py comparison makes accuracy results stratified by SNPs and indels
        call compareCallsHappy {
            input:
                in_sample_vcf_file=bgzipVGCalledChildVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipVGCalledChildVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
    }
    
    output {
        File? output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File? output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_vcf = bgzipVGCalledChildVCF.output_merged_vcf
        File output_vcf_index = bgzipVGCalledChildVCF.output_merged_vcf_index
        Array[File] output_calling_bams = child_calling_bams
        Array[File] output_calling_bam_indexes = child_calling_bam_indexes
        File? output_trio_merged_vcf = convertMergedTrioBCFtoVCF.output_vcf_file
        File? output_trio_merged_vcf_index = convertMergedTrioBCFtoVCF.output_vcf_index
    }   
}

########################
### TASK DEFINITIONS ###
########################

task splitReads {
    input {
        File in_read_file
        String in_pair_id
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
        preemptible: 2
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: "quay.io/glennhickey/pigz:2.3.1"
    }
}

task extractSubsetPathNames {
    input {
        File in_xg_file
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt

        grep -v _decoy path_list.txt | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM > path_list.sub.txt
    }
    output {
        File output_path_list_file = "path_list.sub.txt"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.64.0"
    }
}

task extractReference {
    input {
        File in_xg_file
        File in_path_list_file
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        # Subset to just the paths we care about (may be the whole file) so we
        # get a good dict with just those paths later
        vg paths \
           --extract-fasta \
           -p ${in_path_list_file} \
           --xg ${in_xg_file} > ref.fa
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.64.0"
    }
}

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem
        Int in_index_disk
    }

    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_reference_file} ref.fa
                
        # Index the subset reference
        samtools faidx ref.fa 
        
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict
    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task runVGGIRAFFE {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gbwt_file
        File in_ggbwt_file
        File in_dist_file
        File in_min_file
        File in_ref_dict
        String? in_giraffe_options
        String in_sample_name
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
        vg giraffe \
          --progress \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --sample "~{in_sample_name}" \
          --output-format BAM \
          ~{in_giraffe_options} \
          --ref-paths ~{in_ref_dict} \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -H ~{in_gbwt_file} \
          -g ~{in_ggbwt_file} \
          -d ~{in_dist_file} \
          -m ~{in_min_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
    >>>
    output {
        File chunk_bam_file = glob("*bam")[0]
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.64.0"
    }
}

task sortBAMFile {
    input {
        String in_sample_name
        File in_bam_chunk_file
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
          -O BAM > ~{in_sample_name}.positionsorted.bam
    >>>
    output {
        File sorted_chunk_bam = "~{in_sample_name}.positionsorted.bam"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task indexBAMFile {
    input {
        String in_sample_name
        File in_bam_file
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
        
        # Never use Samtools 1.4.1 here! See https://github.com/samtools/samtools/issues/687
        samtools index \
          ~{in_bam_file} \
          index.bai
    >>>
    output {
        File bam_index = "index.bai"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task fixBAMContigNaming {
    input {
        File in_bam_file
        File in_ref_dict
        String in_prefix_to_strip
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

        # patch the SQ fields from the dict into a new header
        samtools view -H ~{in_bam_file} | grep ^@HD > new_header.sam
        grep ^@SQ ~{in_ref_dict} | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
        samtools view -H ~{in_bam_file}  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam

        # insert the new header, and strip all instances of the prefix
        samtools reheader -P new_header.sam  ~{in_bam_file} | \
          samtools view -h | \
          sed -e "s/~{in_prefix_to_strip}//g" | \
          samtools view --threads ~{in_map_cores} -O BAM > fixed.bam   
    >>>
    output {
        File fixed_bam_file = "fixed.bam"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task fixPathNames {
    input {
        File in_path_file
        String in_prefix_to_strip
     }

     command <<<
        sed -e "s/~{in_prefix_to_strip}//g" ~{in_path_file}  > "fixed_names_file"
     >>>
    output {
        File fixed_path_list_file = "fixed_names_file"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: 2 + " GB"
        cpu: 1
        disks: "local-disk " + 32 + " SSD"
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
        preemptible: 2
        time: 240
        memory: 5 + " GB"
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

        while read -r contig; do
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < "~{in_path_list_file}"
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
    }
    runtime {
        preemptible: 2
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKRealignerTargetCreator {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_call_disk
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

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        # And the dict must be adjacent to both
        ln -f -s "~{in_reference_dict_file}" reference.dict

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "32" \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed
    >>>
    output {
        File realigner_target_bed = glob("*.bed")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}

task widenRealignmentTargets {
    input {
        File in_target_bed_file
        File in_reference_index_file
        Int in_expansion_bases
        Int in_call_disk
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
        
        BASE_NAME=($(ls ~{in_target_bed_file} | rev | cut -f1 -d'/' | rev | sed s/.bed$//g))

        # Widen the BED regions, but don't escape the chromosomes
        bedtools slop -i "~{in_target_bed_file}" -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > "${BASE_NAME}.widened.bed"
    >>>
    output {
        File output_target_bed_file = glob("*.widened.bed")[0]
    }
    runtime {
        preemptible: 2
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }
}

task runAbraRealigner {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk
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

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java -Xmx20G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads 32
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned*bai")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_call_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task leftShiftBAMFile {
    input {
        String in_sample_name
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk
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

        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        bamleftalign \
            <~{in_bam_file} \
            >~{in_sample_name}.${CONTIG_ID}.left_shifted.bam \
            --fasta-reference reference.fa \
            --compressed
    >>>
    output {
        File left_shifted_bam = glob("~{in_sample_name}.*.left_shifted.bam")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 1
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "biocontainers/freebayes:v1.2.0-2-deb_cv1"
    }
}


task runDeepVariantMakeExamples {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        Int in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
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
        # Files may or may not be indel realigned or left shifted in the names.
        # TODO: move tracking of contig ID to WDL variables!
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
                
        NORM_READS_ARG=""
        if [ ~{in_norm_reads} == true ]; then
          NORM_READS_ARG="--normalize_reads"
        fi

        KEEP_LEGACY_AC_ARG=""
        if [ ~{in_keep_legacy_ac} == true ]; then
          KEEP_LEGACY_AC_ARG="--keep_legacy_allele_counter_behavior"
        fi

        seq 0 $((~{in_call_cores}-1)) | \
        parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref reference.fa \
        --reads input_bam_file.bam \
        --examples ./make_examples.tfrecord@~{in_call_cores}.gz \
        --sample_name ~{in_sample_name} \
        --gvcf ./gvcf.tfrecord@~{in_call_cores}.gz \
        --min_mapping_quality ~{in_min_mapq} \
        ${KEEP_LEGACY_AC_ARG} ${NORM_READS_ARG} \
        --regions ${CONTIG_ID} \
        --task {}
        ls | grep 'make_examples.tfrecord-' | tar -czf 'make_examples.tfrecord.tar.gz' -T -
        ls | grep 'gvcf.tfrecord-' | tar -czf 'gvcf.tfrecord.tar.gz' -T -
    >>>
    output {
        File examples_file = "make_examples.tfrecord.tar.gz"
        File nonvariant_site_tf_file = "gvcf.tfrecord.tar.gz"
    }
    runtime {
        preemptible: 5
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "google/deepvariant:1.3.0"
    }
}

task runDeepVariantCallVariants {
    input {
        String in_sample_name
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
        File? in_model_meta_file
        File? in_model_index_file
        File? in_model_data_file
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
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
        
        tar -xzf ~{in_examples_file}
        tar -xzf ~{in_nonvariant_site_tf_file}
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        # We should use an array here, but that doesn't seem to work the way I
        # usually do them (because of a set -u maybe?)
        if [[ ! -z "~{in_model_meta_file}" ]] ; then
            # Model files must be adjacent and not at arbitrary paths
            ln -f -s "~{in_model_meta_file}" model.meta
            ln -f -s "~{in_model_index_file}" model.index
            ln -f -s "~{in_model_data_file}" model.data-00000-of-00001
        else
            # use default WGS models
            ln -f -s "/opt/models/wgs/model.ckpt.meta" model.meta
            ln -f -s "/opt/models/wgs/model.ckpt.index" model.index
            ln -f -s "/opt/models/wgs/model.ckpt.data-00000-of-00001" model.data-00000-of-00001
        fi
        
        /opt/deepvariant/bin/call_variants \
        --outfile call_variants_output.tfrecord.gz \
        --examples "make_examples.tfrecord@~{in_call_cores}.gz" \
        --checkpoint model && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref reference.fa \
        --infile call_variants_output.tfrecord.gz \
        --nonvariant_site_tfrecord_path "gvcf.tfrecord@~{in_call_cores}.gz" \
        --outfile "~{in_sample_name}_deepvariant.vcf.gz" \
        --gvcf_outfile "~{in_sample_name}_deepvariant.g.vcf.gz"
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deepvariant.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deepvariant.g.vcf.gz"
    }
    runtime {
        preemptible: 5
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "google/deepvariant:1.3.0-gpu"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_call_disk
        Int in_call_mem
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
        preemptible: 2
        time: 60
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task glnexusMergeTrioBCF {
    input {
        File in_child_gvcf_file
        File in_parent1_gvcf_file
        File in_parent2_gvcf_file
        Int in_call_disk
        Int in_call_mem
    }
    command <<<'
        set -eux -o pipefail

        glnexus_cli \
            --config DeepVariant_unfiltered \
            "~{in_child_gvcf_file}" \
            "~{in_parent1_gvcf_file}" \
            "~{in_parent2_gvcf_file}" > merged_trio.bcf
    '>>>
    output {
        File output_bcf_file = "merged_trio.bcf"
    }
    runtime {
        preemptible: 2
        time: 60
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/mlin/glnexus:v1.2.7"
    }
}

task convertBCFtoVCF {
    input {
        String in_sample_name
        File in_bcf_file
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
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

        bcftools view \
          --threads ~{in_call_cores} \
          -O z \
          -o "~{in_sample_name}_trio_merged.vcf.gz" \
          "~{in_bcf_file}"

        bcftools index -t "~{in_sample_name}_trio_merged.vcf.gz"
    }
    output {
        File output_vcf_file = "~{in_sample_name}_trio_merged.vcf.gz"
        File output_vcf_index = "~{in_sample_name}_trio_merged.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        time: 60
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        Int in_call_disk
        Int in_call_mem
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
        preemptible: 2
        time: 30
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.64.0"
    }
}

task buildReferenceTemplate {
    input {
        File in_reference_file
    }
    command <<<
        set -eux -o pipefail
    
        rtg format -o template.sdf "~{in_reference_file}"
        tar -czf template.sdf.tar.gz template.sdf/
    >>>
    output {
        File output_template_archive = "template.sdf.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "realtimegenomics/rtg-tools:3.12.1"
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + 10 + " SSD" 
    }
}

task compareCalls {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_template_archive
        File? in_evaluation_regions_file
        Int in_call_disk
        Int in_call_mem
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Set up template; we assume it drops a "template.sdf"
        tar -xf "~{in_template_archive}"
    
        rtg vcfeval \
            --baseline truth.vcf.gz \
            --calls sample.vcf.gz \
            ~{"--evaluation-regions=" + in_evaluation_regions_file} \
            --template template.sdf \
            --threads 32 \
            --output vcfeval_results
            
        tar -czf vcfeval_results.tar.gz vcfeval_results/
    >>>
    output {
        File output_evaluation_archive = "vcfeval_results.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "realtimegenomics/rtg-tools:3.12.1"
        cpu: 32
        disks: "local-disk " + in_call_disk + " SSD"
        memory: in_call_mem + " GB"
    }
}

task compareCallsHappy {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_reference_file
        File in_reference_index_file
        File? in_evaluation_regions_file
        Int in_call_disk
        Int in_call_mem
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        
        mkdir happy_results
   
        /opt/hap.py/bin/hap.py \
            truth.vcf.gz \
            sample.vcf.gz \
            ~{"-f " + in_evaluation_regions_file} \
            --reference reference.fa \
            --threads 32 \
            --engine=vcfeval \
            -o happy_results/eval
    
        tar -czf happy_results.tar.gz happy_results/
    >>>
    output {
        File output_evaluation_archive = "happy_results.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "jmcdani20/hap.py:v0.3.12"
        cpu: 32
        disks: "local-disk " + in_call_disk + " SSD"
        memory: in_call_mem + " GB"
    }
}

task runDeepTrioMakeExamples {
    input {
        String in_child_name
        String in_maternal_name
        String in_paternal_name
        File in_child_bam_file
        File in_child_bam_file_index
        File in_maternal_bam_file
        File in_maternal_bam_file_index
        File in_paternal_bam_file
        File in_paternal_bam_file_index
        File in_reference_file
        File in_reference_index_file
        Int in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
        String in_deeptrio_docker
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
        
        ln -s ~{in_child_bam_file} input_bam_file.child.bam
        ln -s ~{in_child_bam_file_index} input_bam_file.child.bam.bai
        ln -s ~{in_maternal_bam_file} input_bam_file.maternal.bam
        ln -s ~{in_maternal_bam_file_index} input_bam_file.maternal.bam.bai
        ln -s ~{in_paternal_bam_file} input_bam_file.paternal.bam
        ln -s ~{in_paternal_bam_file_index} input_bam_file.paternal.bam.bai
        ln -f -s ~{in_reference_file} ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        # Files may or may not be indel realigned or left shifted in the names.
        # TODO: move tracking of contig ID to WDL variables!
        CONTIG_ID=($(ls ~{in_child_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_child_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))
        echo ${CONTIG_ID}
        NORM_READS_ARG=""
        if [ ~{in_norm_reads} == true ]; then
          NORM_READS_ARG="--normalize_reads"
        fi

        KEEP_LEGACY_AC_ARG=""
        if [ ~{in_keep_legacy_ac} == true ]; then
          KEEP_LEGACY_AC_ARG="--keep_legacy_allele_counter_behavior"
        fi

        /opt/deepvariant/bin/deeptrio/make_examples \
        --mode calling \
        --ref ref.fna \
        --reads input_bam_file.child.bam \
        --reads_parent1 input_bam_file.paternal.bam \
        --reads_parent2 input_bam_file.maternal.bam \
        --examples ./make_examples.tfrecord@~{in_call_cores}.gz \
        --sample_name ~{in_child_name} \
        --sample_name_parent1 ~{in_paternal_name} \
        --sample_name_parent2 ~{in_maternal_name} \
        --gvcf ./gvcf.tfrecord@~{in_call_cores}.gz \
        --min_mapping_quality ~{in_min_mapq} \
        ${KEEP_LEGACY_AC_ARG} ${NORM_READS_ARG} \
        --regions ${CONTIG_ID} \
        --num_shards ~{in_call_cores}
        #--pileup_image_height_child 60 \
        #--pileup_image_height_parent 40 \
        ls | grep 'make_examples_child.tfrecord-' | tar -czf 'make_examples_child.tfrecord.tar.gz' -T -
        ls | grep 'make_examples_parent1.tfrecord-' | tar -czf 'make_examples_parent1.tfrecord.tar.gz' -T -
        ls | grep 'make_examples_parent2.tfrecord-' | tar -czf 'make_examples_parent2.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_child.tfrecord-' | tar -czf 'gvcf_child.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_parent1.tfrecord-' | tar -czf 'gvcf_parent1.tfrecord.tar.gz' -T -
        ls | grep 'gvcf_parent2.tfrecord-' | tar -czf 'gvcf_parent2.tfrecord.tar.gz' -T -
    >>>
    output { 
        File child_examples_file = "make_examples_child.tfrecord.tar.gz"
        File paternal_examples_file = "make_examples_parent1.tfrecord.tar.gz"
        File maternal_examples_file = "make_examples_parent2.tfrecord.tar.gz"
        File child_nonvariant_site_tf_file = "gvcf_child.tfrecord.tar.gz"
        File paternal_nonvariant_site_tf_file = "gvcf_parent1.tfrecord.tar.gz"
        File maternal_nonvariant_site_tf_file = "gvcf_parent2.tfrecord.tar.gz"
    }
    runtime {
        preemptible: 5
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores 
        disks: "local-disk " + in_call_disk + " SSD"
        docker: in_deeptrio_docker
    }
}

task runDeepTrioCallVariants {
    input {
        String in_sample_name
        String in_sample_type
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
        File? in_model_meta_file
        File? in_model_index_file
        File? in_model_data_file
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
        String DEEPTRIO_DOCKER
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
        tar -xzf ~{in_examples_file}
        tar -xzf ~{in_nonvariant_site_tf_file}
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        # We should use an array here, but that doesn't seem to work the way I
        # usually do them (because of a set -u maybe?)
        if [[ ! -z "~{in_model_meta_file}" ]] ; then
            # Model files must be adjacent and not at arbitrary paths
            ln -f -s "~{in_model_meta_file}" model.meta
            ln -f -s "~{in_model_index_file}" model.index
            ln -f -s "~{in_model_data_file}" model.data-00000-of-00001
        elif [ ~{in_sample_type} == "child" ]; then
            # use default child WGS models
            ln -f -s "/opt/models/deeptrio/wgs/child/model.ckpt.meta" model.meta
            ln -f -s "/opt/models/deeptrio/wgs/child/model.ckpt.index" model.index
            ln -f -s "/opt/models/deeptrio/wgs/child/model.ckpt.data-00000-of-00001" model.data-00000-of-00001
        else
            # use default parent WGS models
            ln -f -s "/opt/models/deeptrio/wgs/parent/model.ckpt.meta" model.meta
            ln -f -s "/opt/models/deeptrio/wgs/parent/model.ckpt.index" model.index
            ln -f -s "/opt/models/deeptrio/wgs/parent/model.ckpt.data-00000-of-00001" model.data-00000-of-00001
        fi
        
        # Define expected examples and nonvariant site tensor flow files 
        if [ ~{in_sample_type} == "child" ]; then
            EXAMPLES_FILE="make_examples_child.tfrecord@~{in_call_cores}.gz"
            NONVARIANT_SITE_FILE="gvcf_child.tfrecord@~{in_call_cores}.gz"
        else
            EXAMPLES_FILE="make_examples_~{in_sample_type}.tfrecord@~{in_call_cores}.gz"
            NONVARIANT_SITE_FILE="gvcf_~{in_sample_type}.tfrecord@~{in_call_cores}.gz"
        fi
        
        /opt/deepvariant/bin/call_variants \
        --outfile call_variants_output.tfrecord.gz \
        --examples ${EXAMPLES_FILE} \
        --checkpoint model && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref reference.fa \
        --infile call_variants_output.tfrecord.gz \
        --nonvariant_site_tfrecord_path ${NONVARIANT_SITE_FILE} \
        --outfile "~{in_sample_name}_deeptrio.vcf.gz" \
        --gvcf_outfile "~{in_sample_name}_deeptrio.g.vcf.gz"
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deeptrio.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deeptrio.g.vcf.gz"
    }
    runtime {
        preemptible: 5
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: DEEPTRIO_DOCKER
    }
}

task runDeepTrioOneCommand {
    input {
        String in_child_name
        String in_maternal_name
        String in_paternal_name
        File in_child_bam_file
        File in_child_bam_file_index
        File in_maternal_bam_file
        File in_maternal_bam_file_index
        File in_paternal_bam_file
        File in_paternal_bam_file_index
        File in_reference_file
        File in_reference_index_file
        String in_region
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
        String in_deeptrio_docker
    }
    command <<<
        set -eux -o pipefail

        ln -f -s ~{in_child_bam_file} input_bam_file.child.bam
        ln -f -s ~{in_child_bam_file_index} input_bam_file.child.bam.bai
        ln -f -s ~{in_maternal_bam_file} input_bam_file.maternal.bam
        ln -f -s ~{in_maternal_bam_file_index} input_bam_file.maternal.bam.bai
        ln -f -s ~{in_paternal_bam_file} input_bam_file.paternal.bam
        ln -f -s ~{in_paternal_bam_file_index} input_bam_file.paternal.bam.bai
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        mkdir -p intermediate_results_dir

        # Validate required inputs are present and non-empty
        for f in \
          input_bam_file.child.bam \
          input_bam_file.child.bam.bai \
          input_bam_file.paternal.bam \
          input_bam_file.paternal.bam.bai \
          input_bam_file.maternal.bam \
          input_bam_file.maternal.bam.bai \
          reference.fa reference.fa.fai; do
          if [ ! -s "$f" ]; then
            echo "Missing or empty required file: $f" >&2
            ls -lah || true
            exit 1
          fi
        done

        # Find DeepTrio runner (script may or may not be executable)
        RUNNER=""
        if [ -x "/opt/deepvariant/bin/deeptrio/run_deeptrio" ]; then
          RUNNER="/opt/deepvariant/bin/deeptrio/run_deeptrio"
        elif [ -f "/opt/deepvariant/bin/deeptrio/run_deeptrio" ]; then
          RUNNER="python3 /opt/deepvariant/bin/deeptrio/run_deeptrio"
        elif [ -x "/opt/deepvariant/bin/deeptrio/run_deeptrio.py" ]; then
          RUNNER="/opt/deepvariant/bin/deeptrio/run_deeptrio.py"
        elif [ -f "/opt/deepvariant/bin/deeptrio/run_deeptrio.py" ]; then
          RUNNER="python3 /opt/deepvariant/bin/deeptrio/run_deeptrio.py"
        else
          echo "Could not locate DeepTrio runner under /opt/deepvariant/bin/deeptrio" >&2
          ls -lah /opt/deepvariant/bin/deeptrio || true
          exit 1
        fi

        eval "$RUNNER" \
          --model_type WGS \
          --ref reference.fa \
          --reads_child input_bam_file.child.bam \
          --reads_parent1 input_bam_file.paternal.bam \
          --reads_parent2 input_bam_file.maternal.bam \
          --output_vcf_child child_deeptrio.vcf.gz \
          --output_vcf_parent1 parent1_deeptrio.vcf.gz \
          --output_vcf_parent2 parent2_deeptrio.vcf.gz \
          --sample_name_child ~{in_child_name} \
          --sample_name_parent1 ~{in_paternal_name} \
          --sample_name_parent2 ~{in_maternal_name} \
          --num_shards ~{in_call_cores} \
          --regions ~{in_region} \
          --intermediate_results_dir intermediate_results_dir \
          --output_gvcf_child child_deeptrio.g.vcf.gz \
          --output_gvcf_parent1 parent1_deeptrio.g.vcf.gz \
          --output_gvcf_parent2 parent2_deeptrio.g.vcf.gz
    >>>
    output {
        File output_child_vcf_file = "child_deeptrio.vcf.gz"
        File output_parent1_vcf_file = "parent1_deeptrio.vcf.gz"
        File output_parent2_vcf_file = "parent2_deeptrio.vcf.gz"
        File output_child_gvcf_file = "child_deeptrio.g.vcf.gz"
        File output_parent1_gvcf_file = "parent1_deeptrio.g.vcf.gz"
        File output_parent2_gvcf_file = "parent2_deeptrio.g.vcf.gz"
    }
    runtime {
        preemptible: 5
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        disks: "local-disk " + in_call_disk + " SSD"
        docker: in_deeptrio_docker
    }
}




