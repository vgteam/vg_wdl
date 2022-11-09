version 1.0

import "giraffe_and_deepvariant.wdl" as main
import "giraffe_and_deepvariant_lite.wdl" as lite
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/deepvariant.wdl" as deepvariant

### giraffe_and_deepvariant_fromGAM.wdl ###
## Author: Jean Monlong
## Description: Surject a GAM and run DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow GiraffeDeepVariantFromGAM {
    input {
        File? INPUT_GAM                                 # Input GAM
        File? INPUT_GAF                                 # Input GAF
        String SAMPLE_NAME                              # The sample name
        Boolean PAIRED_END = true                    # Whether the reads are paired-end.
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.37.0" # VG Container used in the pipeline
        Int READS_PER_CHUNK = 100000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        String REFERENCE_PREFIX = ""                    # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File XG_FILE                                    # Path to .xg index file
        File? DV_MODEL_META                             # .meta file for a custom DeepVariant calling model
        File? DV_MODEL_INDEX                            # .index file for a custom DeepVariant calling model
        File? DV_MODEL_DATA                             # .data-00000-of-00001 file for a custom DeepVariant calling model
        Boolean LEFTALIGN_BAM = true                    # Whether or not to left-align reads in the BAM before DV
        Boolean REALIGN_INDELS = true                   # Whether or not to realign reads near indels before DV
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        Boolean DV_KEEP_LEGACY_AC = true                # Should DV use the legacy allele counter behavior?
        Boolean DV_NORM_READS = false                   # Should DV normalize reads itself?
        Boolean SPLIT_AND_SURJECT = false
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int VG_CORES = 20                               # cores used by vg commands
        Int VG_MEM = 120                                # memory used by vg commands to load indexes etc
        File? REFERENCE_FILE                            # (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.
        File? REFERENCE_INDEX_FILE                      # (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
        File? REFERENCE_DICT_FILE                       # (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. 
    }

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call lite.extractSubsetPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_extract_mem=VG_MEM
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
        call lite.extractReference {
            input:
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=VG_MEM
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call lite.indexReference {
            input:
                in_reference_file=reference_file
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])


    if (SPLIT_AND_SURJECT){
        ##
        ## Split GAM into chunks
        ##

        if(defined(INPUT_GAM)){
            call gautils.splitGAM {
                input:
                in_gam_file=select_first([INPUT_GAM]),
	            in_read_per_chunk=READS_PER_CHUNK
            }
        }
        if(defined(INPUT_GAF)){
            call gautils.splitGAF {
                input:
                in_gaf_file=select_first([INPUT_GAF]),
	            in_read_per_chunk=READS_PER_CHUNK
            }
        }

        Array[File] chunk_files = select_first([splitGAM.gam_chunk_files, splitGAF.gaf_chunk_files])
        
        scatter (gam_chunk_file in chunk_files) {
            call gautils.surjectGAMtoSortedBAM as surjectGAMtoSortedBAM_chunk {
                input:
                in_gam_file=gam_chunk_file,
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_sample_name=SAMPLE_NAME,
                input_is_gaf=defined(INPUT_GAF),
                in_is_paired_end=PAIRED_END,
                in_max_fragment_length=MAX_FRAGMENT_LENGTH,
                in_map_cores=VG_CORES,
                in_map_mem=VG_MEM
            }
        }

        call lite.mergeAlignmentBAMChunks {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=surjectGAMtoSortedBAM_chunk.output_bam_file,
            in_map_cores=16
        }
    }
    if (!SPLIT_AND_SURJECT){
        File ga_file = select_first([INPUT_GAM, INPUT_GAF])
        call gautils.surjectGAMtoSortedBAM as surjectGAMtoSortedBAM_whole {
            input:
            in_gam_file=ga_file,
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_is_paired_end=PAIRED_END,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            make_bam_index=true,
            input_is_gaf=defined(INPUT_GAF),
            in_map_cores=VG_CORES,
            in_map_mem=VG_MEM
        }
    }

    File merged_bam = select_first([mergeAlignmentBAMChunks.merged_bam_file, surjectGAMtoSortedBAM_whole.output_bam_file])
    File merged_bam_index = select_first([mergeAlignmentBAMChunks.merged_bam_file_index, surjectGAMtoSortedBAM_whole.output_bam_index_file])
    
    # Split merged alignment by contigs list
    call lite.splitBAMbyPath { 
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=merged_bam,
            in_merged_bam_file_index=merged_bam_index,
            in_path_list_file=pipeline_path_list_file,
            in_map_cores=VG_CORES
    }

    ##
    ## Call variants with DeepVariant in each contig
    ##
    scatter (deepvariant_caller_input_files in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        if (REFERENCE_PREFIX != "") {
            # use samtools to replace the header contigs with those from our dict.
            # this allows the header to contain contigs that are not in the graph,
            # which is more general and lets CHM13-based graphs be used to call on GRCh38
            # also, strip out contig prefixes in the BAM body
            call main.fixBAMContigNaming {
                input:
                in_bam_file=deepvariant_caller_input_files.left,
                in_ref_dict=reference_dict_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
                in_map_cores=VG_CORES,
                in_map_mem=VG_MEM
            }
        }
        File properly_named_bam_file = select_first([fixBAMContigNaming.fixed_bam_file, deepvariant_caller_input_files.left]) 
        File properly_named_bam_index_file = select_first([fixBAMContigNaming.fixed_bam_index_file, deepvariant_caller_input_files.right]) 

        ## Eventually shift and realign reads
        if (LEFTALIGN_BAM){
            call lite.leftShiftBAMFile {
                input:
                in_bam_file=properly_named_bam_file,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }
        
        File current_bam = select_first([leftShiftBAMFile.output_bam_file, properly_named_bam_file])
        File current_bam_index = select_first([leftShiftBAMFile.output_bam_index_file, properly_named_bam_index_file])
        
        if (REALIGN_INDELS) {
            call lite.prepareRealignTargets {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=current_bam,
                in_bam_index_file=current_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_expansion_bases=REALIGNMENT_EXPANSION_BASES
            }
            call lite.runAbraRealigner {
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
        call deepvariant.runDeepVariantMakeExamples {
            input:
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
        call deepvariant.runDeepVariantCallVariants {
            input:
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

    Int vcf_disk_size = 30 * round(size(runDeepVariantCallVariants.output_vcf_file, 'G')) + 50
    # Merge distributed variant called VCFs
    call main.concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file,
            in_call_disk=vcf_disk_size,
            in_call_mem=10
    }

    Int gvcf_disk_size = 30 * round(size(runDeepVariantCallVariants.output_gvcf_file, 'G')) + 50
    # Merge distributed variant called GVCFs
    call main.concatClippedVCFChunks as concatClippedGVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_gvcf_file,
            in_call_disk=gvcf_disk_size,
            in_call_mem=10
    }
    
    output {
        File output_vcf = concatClippedVCFChunks.output_merged_vcf
        File output_vcf_index = concatClippedVCFChunks.output_merged_vcf_index
        File output_gvcf = concatClippedGVCFChunks.output_merged_vcf
        File output_gvcf_index = concatClippedGVCFChunks.output_merged_vcf_index
    }   
}

########################
### TASK DEFINITIONS ###
########################

task surjectGAMtoSortedIndexedBAM {
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
                                    -O BAM > ~{out_prefix}.bam \
            && samtools index ~{out_prefix}.bam
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
        File output_bam_index_file = "~{out_prefix}.bam.bai"
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
