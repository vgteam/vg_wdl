version 1.0

import "giraffe_and_deepvariant.wdl" as main

### giraffe_and_deepvariant_lite.wdl ###
## Author: Jean Monlong
## Description: Core VG Giraffe mapping and DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMap {
    input {
        File? INPUT_READ_FILE_1                         # Input sample 1st read pair fastq.gz
        File? INPUT_READ_FILE_2                         # Input sample 2nd read pair fastq.gz
        File? INPUT_CRAM_FILE                           # Input CRAM file
        File? CRAM_REF                                  # Genome fasta file associated with the CRAM file
        File? CRAM_REF_INDEX                            # Index of the fasta file associated with the CRAM file
        String SAMPLE_NAME                              # The sample name
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.37.0" # VG Container used in the pipeline
        Int READS_PER_CHUNK = 30000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        File XG_FILE                                    # Path to .xg index file
        File GBWT_FILE                                  # Path to .gbwt index file
        File GGBWT_FILE                                 # Path to .gg index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
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
        Boolean SORT_GAM = false                        # Should the GAM be sorted
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 60
        Int MAP_CORES = 20
        Int MAP_MEM = 120
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        File? REFERENCE_FILE                            # (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.
        File? REFERENCE_INDEX_FILE                      # (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
        File? REFERENCE_DICT_FILE                       # (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. 
    }

    if(defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	call main.convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_cores=SPLIT_READ_CORES
	}
    }

    File read_1_file = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file])
    File read_2_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file])
    
    # Split input reads into chunks for parallelized mapping
    call main.splitReads as firstReadPair {
        input:
            in_read_file=read_1_file,
            in_pair_id="1",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }
    call main.splitReads as secondReadPair {
        input:
            in_read_file=read_2_file,
            in_pair_id="2",
            in_vg_container=VG_CONTAINER,
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call extractSubsetPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
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
        call extractReference {
            input:
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_vg_container=VG_CONTAINER,
                in_extract_mem=MAP_MEM
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call indexReference {
            input:
                in_reference_file=reference_file
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])

    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
    scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
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
            in_sample_name=SAMPLE_NAME,
            in_map_cores=MAP_CORES,
            in_map_mem=MAP_MEM
        }
        call surjectGAMtoSortedBAM {
            input:
            in_gam_file=runVGGIRAFFE.chunk_gam_file,
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_vg_container=VG_CONTAINER,
            in_map_cores=MAP_CORES,
            in_map_mem=MAP_MEM
        }
    }

    call mergeAlignmentBAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=surjectGAMtoSortedBAM.output_bam_file,
        in_map_cores=16
    }
    
    # Split merged alignment by contigs list
    call splitBAMbyPath { 
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
            in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
            in_path_list_file=pipeline_path_list_file,
            in_map_cores=MAP_CORES
    }

    ##
    ## Call variants with DeepVariant in each contig
    ##
    scatter (deepvariant_caller_input_files in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        ## Eventually shift and realign reads
        if (LEFTALIGN_BAM){
            call leftShiftBAMFile {
                input:
                in_bam_file=deepvariant_caller_input_files.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }
        File current_bam = select_first([leftShiftBAMFile.output_bam_file, deepvariant_caller_input_files.left])
        File current_bam_index = select_first([leftShiftBAMFile.output_bam_index_file, deepvariant_caller_input_files.right])
        if (REALIGN_INDELS) {
            call prepareRealignTargets {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=current_bam,
                in_bam_index_file=current_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_expansion_bases=REALIGNMENT_EXPANSION_BASES
            }
            call runAbraRealigner {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=current_bam,
                in_bam_index_file=current_bam_index,
                in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }
        File calling_bam = select_first([runAbraRealigner.indel_realigned_bam, current_bam])
        File calling_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, current_bam_index])
        ## DeepVariant calling
        call runDeepVariantMakeExamples {
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
                in_call_mem=CALL_MEM
        }
        call runDeepVariantCallVariants {
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
                in_call_mem=CALL_MEM
        }
    }

    Int vcf_disk_size =  * round(size(runDeepVariantCallVariants.output_vcf_file, 'G')) + 50
    # Merge distributed variant called VCFs
    call main.concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file,
            in_call_disk=vcf_disk_size,
            in_call_mem=10
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call main.bgzipMergedVCF {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_call_disk=vcf_disk_size,
            in_call_mem=10
    }

    Int gvcf_disk_size = 5 * round(size(runDeepVariantCallVariants.output_gvcf_file, 'G')) + 50
    # Merge distributed variant called GVCFs
    call main.concatClippedVCFChunks as concatClippedGVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_gvcf_file,
            in_call_disk=gvcf_disk_size,
            in_call_mem=10
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call main.bgzipMergedVCF as bgzipMergedGVCF{
        input:
            in_sample_name=SAMPLE_NAME + '.g',
            in_merged_vcf_file=concatClippedGVCFChunks.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_call_disk=gvcf_disk_size,
            in_call_mem=10
    }
    
    if (SORT_GAM){
        call mergeGAMandSort {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gam_chunk_files=runVGGIRAFFE.chunk_gam_file,
            in_vg_container=VG_CONTAINER,
            in_mem=MAP_MEM,
            in_cores=MAP_CORES
        }
    }
    if (!SORT_GAM) {
        call mergeGAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gam_chunk_files=runVGGIRAFFE.chunk_gam_file,
            in_vg_container=VG_CONTAINER
        }
    }

    File gam_file = select_first([mergeGAM.output_merged_gam, mergeGAMandSort.output_merged_gam])

    output {
        File output_vcf = bgzipMergedVCF.output_merged_vcf
        File output_vcf_index = bgzipMergedVCF.output_merged_vcf_index
        File output_gvcf = bgzipMergedGVCF.output_merged_vcf
        File output_gvcf_index = bgzipMergedGVCF.output_merged_vcf_index
        File output_gam = gam_file
        File? output_gam_index = mergeGAMandSort.output_merged_gam_index
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
        Int in_nb_chunks
        Int in_max_chunks
    }
    Int disk_size = round(5 * size(in_cram_file, 'G')) + 20
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
        disks: "local-disk " + disk_size + " SSD"
        docker: "jmonlong/samtools-jm:release-1.19jm0.2.2"
        preemptible: 3
    }
}

task extractSubsetPathNames {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_mem
    }
    Int disk_size = round(3 * size(in_xg_file, 'G')) + 20
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
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }
}

task extractReference {
    input {
        File in_xg_file
        File in_path_list_file
        String in_vg_container
        Int in_extract_mem
    }
    Int disk_size = round(3 * size(in_xg_file, 'G')) + 20
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
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }
}

task indexReference {
    input {
        File in_reference_file
    }
    Int disk_size = round(3 * size(in_reference_file, 'G')) + 20
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
        memory: "20 GB"
        disks: "local-disk " + disk_size + " SSD"
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
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        String in_map_mem
    }
    Int disk_size = 3 * round(size(in_xg_file, 'G') + size(in_gbwt_file, 'G') + size(in_ggbwt_file, 'G') + size(in_dist_file, 'G') + size(in_min_file, 'G') + size(in_left_read_pair_chunk_file, 'G') + size(in_right_read_pair_chunk_file, 'G')) + 20
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
          --output-format gam \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -H ~{in_gbwt_file} \
          -g ~{in_ggbwt_file} \
          -d ~{in_dist_file} \
          -m ~{in_min_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*gam")[0]
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }
}

task surjectGAMtoSortedBAM {
    input {
        File in_gam_file
        File in_xg_file
        File in_path_list_file
        String in_sample_name
        Int in_max_fragment_length
        String in_vg_container
        Int in_map_cores
        String in_map_mem
    }
    String out_prefix = basename(in_gam_file, ".gam")
    Int half_cores = in_map_cores / 2
    Int disk_size = 3 * round(size(in_xg_file, 'G') + size(in_gam_file, 'G')) + 20
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

        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_xg_file} \
          -t ~{half_cores} \
          --bam-output \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx \
          --interleaved --max-frag-len ~{in_max_fragment_length} \
          ~{in_gam_file} | samtools sort --threads ~{half_cores} \
                                    -O BAM > ~{out_prefix}.bam
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }

}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Int in_map_cores
    }
    Int disk_size = round(3 * size(in_alignment_bam_chunk_files, 'G')) + 20
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
        memory: "5 GB"
        cpu: in_map_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task mergeGAM {
    input {
        String in_sample_name
        Array[File] in_gam_chunk_files
        String in_vg_container
    }
    Int disk_size = round(4 * size(in_gam_chunk_files, 'G')) + 20
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

        cat ~{sep=" " in_gam_chunk_files} > ~{in_sample_name}.gam
    >>>
    output {
        File output_merged_gam = "~{in_sample_name}.gam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: "10 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }
}

task mergeGAMandSort {
    input {
        String in_sample_name
        Array[File] in_gam_chunk_files
        String in_vg_container
        Int in_cores
        Int in_mem
    }
    Int disk_size = round(4 * size(in_gam_chunk_files, 'G')) + 20
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

        cat ~{sep=" " in_gam_chunk_files} | vg gamsort -p -i ~{in_sample_name}.sorted.gam.gai \
                                               -t ~{in_cores} - > ~{in_sample_name}.sorted.gam
    >>>
    output {
        File output_merged_gam = "~{in_sample_name}.sorted.gam"
        File output_merged_gam_index = "~{in_sample_name}.sorted.gam.gai"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: in_vg_container
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
    }
    Int disk_size = round(3 * size(in_merged_bam_file, 'G')) + 20
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
        memory: "20 GB"
        cpu: in_map_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task leftShiftBAMFile {
    input {
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
    }
    String out_prefix = basename(in_bam_file, ".bam")
    Int disk_size = round(3 * size(in_bam_file, 'G')) + 50
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

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        bamleftalign \
            < ~{in_bam_file} \
            > ~{out_prefix}.left_shifted.bam \
            --fasta-reference reference.fa \
            --compressed
        samtools index -b ~{out_prefix}.left_shifted.bam ~{out_prefix}.left_shifted.bam.bai
    >>>
    output {
        File output_bam_file = "~{out_prefix}.left_shifted.bam"
        File output_bam_index_file = "~{out_prefix}.left_shifted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 180
        memory: "20 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/freebayes-samtools:1.2.0_1.10"
    }
}

task prepareRealignTargets {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_expansion_bases
    }
    Int disk_size = round(2 * size(in_bam_file, 'G')) + 20
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
          -nt 16 \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed

        if [ ~{in_expansion_bases} -gt 0 ]; then
            bedtools slop -i ~{in_sample_name}.${CONTIG_ID}.intervals.bed -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > ~{in_sample_name}.${CONTIG_ID}.intervals.widened.bed
            mv ~{in_sample_name}.${CONTIG_ID}.intervals.widened.bed ~{in_sample_name}.${CONTIG_ID}.intervals.bed
        fi
    >>>
    output {
        File output_target_bed_file = glob("*.bed")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0"
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
    }
    Int disk_size = round(3 * size(in_bam_file, 'G')) + 20
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
          --threads 16
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
        disks: "local-disk " + disk_size + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task runDeepVariantMakeExamples {
    input {
        String in_dv_container
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        Int in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        Int in_call_cores
        Int in_call_mem
    }
    Int disk_size = round(2 * size(in_bam_file, 'G')) + 20
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
        disks: "local-disk " + disk_size + " SSD"
        docker: in_dv_container
    }
}

task runDeepVariantCallVariants {
    input {
        String in_dv_gpu_container
        String in_sample_name
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
        File? in_model_meta_file
        File? in_model_index_file
        File? in_model_data_file
        Int in_call_cores
        Int in_call_mem
    }
    Int disk_size = 3 * round(size(in_examples_file, 'G') + size(in_nonvariant_site_tf_file, 'G') + size(in_reference_file, 'G')) + 20    
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
        disks: "local-disk " + disk_size + " SSD"
        docker: in_dv_gpu_container
    }
}

