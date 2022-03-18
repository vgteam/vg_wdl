version 1.0

import "giraffe_and_deepvariant.wdl" as main
import "giraffe_and_deepvariant_lite.wdl" as lite

### giraffe_and_deepvariant_fromGAM.wdl ###
## Author: Jean Monlong
## Description: Surject a GAM and run DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgMultiMap {
    input {
        File INPUT_GAM                                 # Input GAM
        String SAMPLE_NAME                              # The sample name
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.37.0" # VG Container used in the pipeline
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        File XG_FILE                                    # Path to .xg index file
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
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int VG_CORES = 16                               # cores used by vg commands
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
                    in_vg_container=VG_CONTAINER,
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
                in_vg_container=VG_CONTAINER,
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

    ##
    ## Split GAM into one GAM per chromosome
    ##

    call splitGAMbyChrom {
        input:
        in_gam_file=INPUT_GAM,
        in_xg_file=XG_FILE,
        in_path_list_file=pipeline_path_list_file,
        in_sample_name=SAMPLE_NAME,
        in_mem=VG_MEM,
        in_cores=VG_CORES,
        in_vg_container=VG_CONTAINER,
    }

    ##
    ## Call variants with DeepVariant in each contig
    ##
    scatter (gam_contig_file in splitGAMbyChrom.gam_contig_files) {
        call surjectGAMtoSortedIndexedBAM {
            input:
            in_gam_file=gam_contig_file,
            in_xg_file=XG_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_vg_container=VG_CONTAINER,
            in_cores=VG_CORES,
            in_mem=VG_MEM
        }

        ## Eventually shift and realign reads
        if (LEFTALIGN_BAM){
            call lite.leftShiftBAMFile {
                input:
                in_bam_file=surjectGAMtoSortedIndexedBAM.output_bam_file,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }
        File current_bam = select_first([leftShiftBAMFile.output_bam_file, surjectGAMtoSortedIndexedBAM.output_bam_file])
        File current_bam_index = select_first([leftShiftBAMFile.output_bam_index_file, surjectGAMtoSortedIndexedBAM.output_bam_index_file])
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
        call lite.runDeepVariantMakeExamples {
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
        call lite.runDeepVariantCallVariants {
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

    # Merge distributed variant called VCFs
    call main.concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file,
            in_call_disk=50,
            in_call_mem=20
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call main.bgzipMergedVCF {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_call_disk=20,
            in_call_mem=20
    }
        
    # Merge distributed variant called GVCFs
    call main.concatClippedVCFChunks as concatClippedGVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_gvcf_file,
            in_call_disk=200,
            in_call_mem=20
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call main.bgzipMergedVCF as bgzipMergedGVCF{
        input:
            in_sample_name=SAMPLE_NAME + '.g',
            in_merged_vcf_file=concatClippedGVCFChunks.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_call_disk=100,
            in_call_mem=20
    }
    
    output {
        File output_vcf = bgzipMergedVCF.output_merged_vcf
        File output_vcf_index = bgzipMergedVCF.output_merged_vcf_index
        File output_gvcf = bgzipMergedGVCF.output_merged_vcf
        File output_gvcf_index = bgzipMergedGVCF.output_merged_vcf_index
    }   
}

########################
### TASK DEFINITIONS ###
########################

task splitGAMbyChrom {
    input {
        File in_gam_file
        File in_xg_file
        File in_path_list_file
        String in_sample_name
        Int in_mem
        Int in_cores
        String in_vg_container
    }
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

        vg chunk -C -x ~{in_xg_file} -a ~{in_gam_file} -b ~{in_sample_name} -P ~{in_path_list_file}
    >>>
    output {
        Array[File] gam_contig_files = glob("~{in_sample_name}*.gam")
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

task surjectGAMtoSortedIndexedBAM {
    input {
        File in_gam_file
        File in_xg_file
        File in_path_list_file
        String in_sample_name
        Int in_max_fragment_length
        String in_vg_container
        Int in_cores
        String in_mem
    }
    String out_prefix = basename(in_gam_file, ".gam")
    Int half_cores = in_cores / 2
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

        samtools index ~{out_prefix}.bam
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
        File output_bam_index_file = "~{out_prefix}.bam.bai"
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
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}_//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

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
