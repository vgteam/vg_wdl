version 1.0

### vg_deeptrio_calling_workflow.wdl ###
## Author: Charles Markello
## Description: Core VG DEEPTRIO workflow for calling variants on trio datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow vgDeepTrioCall {
    input {
        Boolean ABRA_REALIGN = false                                # Set to 'true' to use ABRA2 IndelRealigner instead of GATK for indel realignment
        Boolean SMALL_RESOURCES = false                             # Set to 'true' to use small resources for tiny test dataset
        File? MATERNAL_BAM_FILE                                     # Input maternal surjected .bam file
        File? MATERNAL_BAM_FILE_INDEX                               # Input maternal .bai index of surjected .bam file.
        Array[File]? MATERNAL_BAM_CONTIG_LIST                       # Input maternal bam per contig in a list of files.
        Array[File]? MATERNAL_BAM_INDEX_CONTIG_LIST                 # Input maternal bam index per contig in a list of files following the same contig order as MATERNAL_BAM_CONTIG_LIST.
        File? PATERNAL_BAM_FILE                                     # Input paternal surjected .bam file
        File? PATERNAL_BAM_FILE_INDEX                               # Input paternal .bai index of surjected .bam file.
        Array[File]? PATERNAL_BAM_CONTIG_LIST                       # Input paternal bam per contig in a list of files.
        Array[File]? PATERNAL_BAM_INDEX_CONTIG_LIST                 # Input paternal bam index per contig in a list of files following the same contig order as PATERNAL_BAM_CONTIG_LIST.
        File? CHILD_BAM_FILE                                        # Input child surjected .bam file
        File? CHILD_BAM_FILE_INDEX                                  # Input child .bai index of surjected .bam file.
        String SAMPLE_NAME_CHILD                                    # The child sample name
        String SAMPLE_NAME_MATERNAL                                 # The maternal sample name
        String SAMPLE_NAME_PATERNAL                                 # The paternal sample name
        Array[String]+ CONTIGS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        File REF_FILE                                               # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                         # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                          # Path to .dict file of the REF_FILE fasta reference
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.34.0"           # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        File? DEEPTRIO_CHILD_MODEL
        File? DEEPTRIO_PARENT_MODEL
        File? DEEPVAR_MODEL
        Int VGCALL_CORES = 32
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
    }
    # Split merged alignment by contigs list
    call splitBAMbyPath as splitChildBAMbyPath {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_merged_bam_file=CHILD_BAM_FILE,
            in_merged_bam_file_index=CHILD_BAM_FILE_INDEX,
            contigs=CONTIGS,
            in_small_resources=SMALL_RESOURCES
    }
    if (!defined(MATERNAL_BAM_CONTIG_LIST)) {
        call splitBAMbyPath as splitMaternalBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_merged_bam_file=MATERNAL_BAM_FILE,
                in_merged_bam_file_index=MATERNAL_BAM_FILE_INDEX,
                contigs=CONTIGS,
                in_small_resources=SMALL_RESOURCES
        }
    }
    if (!defined(PATERNAL_BAM_CONTIG_LIST)) {
        call splitBAMbyPath as splitPaternalBAMbyPath {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_merged_bam_file=PATERNAL_BAM_FILE,
                in_merged_bam_file_index=PATERNAL_BAM_FILE_INDEX,
                contigs=CONTIGS,
                in_small_resources=SMALL_RESOURCES
        }
    }
    # Run distributed DeepTRIO linear variant calling for each chromosomal contig
    Array[File] maternal_bams_by_contig = select_first([MATERNAL_BAM_CONTIG_LIST, splitMaternalBAMbyPath.bam_contig_files])
    Array[File] maternal_bam_indexes_by_contig = select_first([MATERNAL_BAM_INDEX_CONTIG_LIST, splitMaternalBAMbyPath.bam_contig_files_index])
    Array[Pair[File, File]] maternal_bams_and_indexes_by_contig = zip(maternal_bams_by_contig, maternal_bam_indexes_by_contig)
    Array[File] paternal_bams_by_contig = select_first([PATERNAL_BAM_CONTIG_LIST, splitPaternalBAMbyPath.bam_contig_files])
    Array[File] paternal_bam_indexes_by_contig = select_first([PATERNAL_BAM_INDEX_CONTIG_LIST, splitPaternalBAMbyPath.bam_contig_files_index])
    Array[Pair[File, File]] paternal_bams_and_indexes_by_contig = zip(paternal_bams_by_contig, paternal_bam_indexes_by_contig)
    Array[Pair[File, File]] child_bams_and_indexes_by_contig = zip(splitChildBAMbyPath.bam_contig_files, splitChildBAMbyPath.bam_contig_files_index)
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
        if (ABRA_REALIGN) {
            call runGATKRealignerTargetCreator as realignTargetCreateChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_bam_file,
                    in_bam_index_file=child_bam_file_index,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE
            }
            call runAbraRealigner as abraRealignChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_bam_file,
                    in_bam_index_file=child_bam_file_index,
                    in_target_bed_file=realignTargetCreateChild.realigner_target_bed,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE
            }
            # Only do indel realignment of parents if there isn't already a supplied indel realigned bam set
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKRealignerTargetCreator as realignTargetCreateMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_bam_file,
                        in_bam_index_file=maternal_bam_file_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
                call runAbraRealigner as abraRealignMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_bam_file,
                        in_bam_index_file=maternal_bam_file_index,
                        in_target_bed_file=realignTargetCreateMaternal.realigner_target_bed,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKRealignerTargetCreator as realignTargetCreatePaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_bam_file,
                        in_bam_index_file=paternal_bam_file_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
                call runAbraRealigner as abraRealignPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_bam_file,
                        in_bam_index_file=paternal_bam_file_index,
                        in_target_bed_file=realignTargetCreatePaternal.realigner_target_bed,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE
                }
            }
        }
        if (!ABRA_REALIGN) {
            call runGATKIndelRealigner as gatkRealignChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_bam_file,
                    in_bam_index_file=child_bam_file_index,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE,
                    in_small_resources=SMALL_RESOURCES
            }
            # Only do indel realignment of parents if there isn't already a supplied indel realigned bam set
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKIndelRealigner as gatkRealignMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_bam_file,
                        in_bam_index_file=maternal_bam_file_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE,
                        in_small_resources=SMALL_RESOURCES
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runGATKIndelRealigner as gatkRealignPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_bam_file,
                        in_bam_index_file=paternal_bam_file_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE,
                        in_small_resources=SMALL_RESOURCES
                }
            }
        }
        File child_indel_realigned_bam = select_first([abraRealignChild.indel_realigned_bam, gatkRealignChild.indel_realigned_bam])
        File child_indel_realigned_bam_index = select_first([abraRealignChild.indel_realigned_bam_index, gatkRealignChild.indel_realigned_bam_index])
        File maternal_indel_realigned_bam = select_first([abraRealignMaternal.indel_realigned_bam, gatkRealignMaternal.indel_realigned_bam, maternal_bam_file])
        File maternal_indel_realigned_bam_index = select_first([abraRealignMaternal.indel_realigned_bam_index, gatkRealignMaternal.indel_realigned_bam_index, maternal_bam_file_index])
        File paternal_indel_realigned_bam = select_first([abraRealignPaternal.indel_realigned_bam, gatkRealignPaternal.indel_realigned_bam, paternal_bam_file])
        File paternal_indel_realigned_bam_index = select_first([abraRealignPaternal.indel_realigned_bam_index, gatkRealignPaternal.indel_realigned_bam_index, paternal_bam_file_index])
        
        #TODO
        String contig_name = sub(sub(sub(child_indel_realigned_bam, ".indel_realigned.bam", ""), SAMPLE_NAME_CHILD, ""), ".", "")
        if ((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT")) {
            call runDeepVariant as callDeepVariantChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_bam_file=child_indel_realigned_bam,
                    in_bam_file_index=child_indel_realigned_bam_index,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_model=DEEPVAR_MODEL,
                    in_small_resources=SMALL_RESOURCES
            } 
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runDeepVariant as callDeepVariantMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_bam_file=maternal_indel_realigned_bam,
                        in_bam_file_index=maternal_indel_realigned_bam_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_model=DEEPVAR_MODEL,
                        in_small_resources=SMALL_RESOURCES
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runDeepVariant as callDeepVariantPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_bam_file=paternal_indel_realigned_bam,
                        in_bam_file_index=paternal_indel_realigned_bam_index,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_model=DEEPVAR_MODEL,
                        in_small_resources=SMALL_RESOURCES
                }
            }
        }
        if (!((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT"))) {
            call runDeepTrioMakeExamples {
                input:
                    in_child_name=SAMPLE_NAME_CHILD,
                    in_maternal_name=SAMPLE_NAME_MATERNAL,
                    in_paternal_name=SAMPLE_NAME_PATERNAL,
                    in_child_bam_file=child_indel_realigned_bam,
                    in_child_bam_file_index=child_indel_realigned_bam_index,
                    in_maternal_bam_file=maternal_indel_realigned_bam,
                    in_maternal_bam_file_index=maternal_indel_realigned_bam_index,
                    in_paternal_bam_file=paternal_indel_realigned_bam,
                    in_paternal_bam_file_index=paternal_indel_realigned_bam_index,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_small_resources=SMALL_RESOURCES
            }
            call runDeepTrioCallVariants as callVariantsChild {
                input:
                    in_sample_name=SAMPLE_NAME_CHILD,
                    in_sample_type="child",
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_examples_file=runDeepTrioMakeExamples.child_examples_file,
                    in_nonvariant_site_tf_file=runDeepTrioMakeExamples.child_nonvariant_site_tf_file,
                    in_model=DEEPTRIO_CHILD_MODEL,
                    in_small_resources=SMALL_RESOURCES
            }
            if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runDeepTrioCallVariants as callVariantsMaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_MATERNAL,
                        in_sample_type="parent2",
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_examples_file=runDeepTrioMakeExamples.maternal_examples_file,
                        in_nonvariant_site_tf_file=runDeepTrioMakeExamples.maternal_nonvariant_site_tf_file,
                        in_model=DEEPTRIO_PARENT_MODEL,
                        in_small_resources=SMALL_RESOURCES
                }
            }
            if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
                call runDeepTrioCallVariants as callVariantsPaternal {
                    input:
                        in_sample_name=SAMPLE_NAME_PATERNAL,
                        in_sample_type="parent1",
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_examples_file=runDeepTrioMakeExamples.paternal_examples_file,
                        in_nonvariant_site_tf_file=runDeepTrioMakeExamples.paternal_nonvariant_site_tf_file,
                        in_model=DEEPTRIO_PARENT_MODEL,
                        in_small_resources=SMALL_RESOURCES
                }
            }
        }
    }
    Array[File] childDeepVarGVCF = select_all(callDeepVariantChild.output_gvcf_file)
    Array[File] childDeepTrioGVCF = select_all(callVariantsChild.output_gvcf_file)
    Array[File] child_contig_gvcf_output_list = select_all(flatten([childDeepTrioGVCF, childDeepVarGVCF]))
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks as concatVCFChunksChild {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_clipped_vcf_chunk_files=child_contig_gvcf_output_list,
            in_small_resources=SMALL_RESOURCES
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF as bgzipVGCalledChildVCF {
        input:
            in_sample_name=SAMPLE_NAME_CHILD,
            in_merged_vcf_file=concatVCFChunksChild.output_merged_vcf,
            in_vg_container=VG_CONTAINER,
            in_small_resources=SMALL_RESOURCES
    }
    if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] maDeepVarGVCF = select_all(callDeepVariantMaternal.output_gvcf_file)
        Array[File] maDeepTrioGVCF = select_all(callVariantsMaternal.output_gvcf_file)
        Array[File] maternal_contig_gvcf_output_list = select_all(flatten([maDeepTrioGVCF, maDeepVarGVCF]))
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksMaternal {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_clipped_vcf_chunk_files=maternal_contig_gvcf_output_list,
                in_small_resources=SMALL_RESOURCES
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledMaternalVCF {
            input:
                in_sample_name=SAMPLE_NAME_MATERNAL,
                in_merged_vcf_file=concatVCFChunksMaternal.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_small_resources=SMALL_RESOURCES
        }
    }
    if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] paDeepVarGVCF = select_all(callDeepVariantPaternal.output_gvcf_file)
        Array[File] paDeepTrioGVCF = select_all(callVariantsPaternal.output_gvcf_file)
        Array[File] paternal_contig_gvcf_output_list = select_all(flatten([paDeepTrioGVCF, paDeepVarGVCF]))
        # Merge distributed variant called VCFs
        call concatClippedVCFChunks as concatVCFChunksPaternal {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_clipped_vcf_chunk_files=paternal_contig_gvcf_output_list,
                in_small_resources=SMALL_RESOURCES
        }
        # Extract either the normal or structural variant based VCFs and compress them
        call bgzipMergedVCF as bgzipVGCalledPaternalVCF {
            input:
                in_sample_name=SAMPLE_NAME_PATERNAL,
                in_merged_vcf_file=concatVCFChunksPaternal.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_small_resources=SMALL_RESOURCES
        }
    }
    
    # Collect parental indel-realigned BAM contig lists
    if (!defined(MATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] maternal_indelrealigned_bams_by_contig = maternal_indel_realigned_bam
        Array[File] maternal_indelrealigned_bam_indexes_by_contig = maternal_indel_realigned_bam_index
    }
    if (!defined(PATERNAL_BAM_INDEX_CONTIG_LIST)) {
        Array[File] paternal_indelrealigned_bams_by_contig = paternal_indel_realigned_bam
        Array[File] paternal_indelrealigned_bam_indexes_by_contig = paternal_indel_realigned_bam_index
    }
    
    output {
        File output_child_gvcf = bgzipVGCalledChildVCF.output_merged_vcf
        File output_child_gvcf_index = bgzipVGCalledChildVCF.output_merged_vcf_index
        File? output_maternal_gvcf = bgzipVGCalledMaternalVCF.output_merged_vcf
        File? output_maternal_gvcf_index = bgzipVGCalledMaternalVCF.output_merged_vcf_index
        File? output_paternal_gvcf = bgzipVGCalledPaternalVCF.output_merged_vcf
        File? output_paternal_gvcf_index = bgzipVGCalledPaternalVCF.output_merged_vcf_index
        Array[File] output_child_indelrealigned_bams = child_indel_realigned_bam
        Array[File] output_child_indelrealigned_bam_indexes = child_indel_realigned_bam_index
        Array[File]? output_maternal_indelrealigned_bams = maternal_indelrealigned_bams_by_contig
        Array[File]? output_maternal_indelrealigned_bam_indexes = maternal_indelrealigned_bam_indexes_by_contig
        Array[File]? output_paternal_indelrealigned_bams = paternal_indelrealigned_bams_by_contig
        Array[File]? output_paternal_indelrealigned_bam_indexes = paternal_indelrealigned_bam_indexes_by_contig
    }
}
########################
### TASK DEFINITIONS ###
########################
task splitBAMbyPath {
    input {
        String in_sample_name
        File? in_merged_bam_file
        File? in_merged_bam_file_index
        Array[String]+ contigs
        Boolean in_small_resources
    }
    
    Int in_cores = if in_small_resources then 2 else 32
    Int in_disk = if in_small_resources then 1 else 10
    String in_mem = if in_small_resources then "1" else "40"
    
    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai
        
        while read -r contig; do
            samtools view \
              -@ ~{in_cores} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < "~{write_lines(contigs)}"
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
    }
    runtime {
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        ln -f -s ~{in_reference_dict_file} ref.dict
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g))
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \ 
          --remove_program_records \ 
          -drf DuplicateRead \ 
          --disable_bam_indexing \ 
          -nt "32" \ 
          -R ref.fna \ 
          -L ${CONTIG_ID} \ 
          -I input_bam_file.bam \ 
          --out forIndelRealigner.intervals 
         
        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed 
    >>> 
    output { 
        File realigner_target_bed = glob("*.bed")[0] 
    } 
    runtime { 
        memory: 50 + " GB" 
        cpu: 32 
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b" 
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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g))
        java -Xmx20G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref ref.fna \
          --index \
          --threads 32
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bai")[0]
    }
    runtime {
        memory: 50 + " GB"
        cpu: 32
        docker: "dceoy/abra2:latest"
    }
}

task runGATKIndelRealigner {
    input {
        String in_sample_name
        File in_bam_file 
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Boolean in_small_resources
    }
    
    Int in_cores = if in_small_resources then 4 else 32
    Int in_disk = if in_small_resources then 1 else 50
    String in_mem = if in_small_resources then "1" else "50"
    
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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        ln -f -s ~{in_reference_dict_file} ref.dict
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g))
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "~{in_cores}" \
          -R ref.fna \
          -I input_bam_file.bam \
          -L ${CONTIG_ID} \
          --out forIndelRealigner.intervals \
        && java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
          --remove_program_records \
          -R ref.fna \
          --targetIntervals forIndelRealigner.intervals \
          -I input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bai")[0]
    }
    runtime {
        memory: in_mem + " GB"
        cpu: in_cores
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}

task runDeepVariant {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File? in_model
        Boolean in_small_resources
    }
    
    Int in_call_cores = if in_small_resources then 4 else 8
    Int in_call_disk = if in_small_resources then 2 else 40
    String in_call_mem = if in_small_resources then "5" else "64"
    
    Boolean custom_model = defined(in_model)
    
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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.indel_realigned.bam$//g))
        DEEPVAR_MODEL="/opt/models/wgs/model.ckpt"
        if [ ~{custom_model} == true ]; then
            tar -xzf ~{in_model}
            model_basename=$(basename ~{in_model})
            stripped_basename=${model_basename%%.*}
            core_filename=$(ls $stripped_basename | head -n 1)
            common_filename=${core_filename%.*}
            DEEPVAR_MODEL="${PWD}/${stripped_basename}/${common_filename}"
        fi

        /opt/deepvariant/bin/run_deepvariant \
        --make_examples_extra_args 'min_mapping_quality=1' \
        --model_type WGS \
        --customized_model ${DEEPVAR_MODEL} \
        --regions ${CONTIG_ID} \
        --ref ref.fna \
        --reads input_bam_file.bam \
        --output_vcf "~{in_sample_name}_deeptrio.vcf.gz" \
        --output_gvcf "~{in_sample_name}_deeptrio.g.vcf.gz" \
        --intermediate_results_dir tmp_deepvariant \
        --num_shards=~{in_call_cores}
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deeptrio.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deeptrio.g.vcf.gz"
    }
    runtime {
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        preemptible: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "google/deepvariant:1.1.0-gpu"
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
        Boolean in_small_resources
    }
    
    Int in_vgcall_cores = if in_small_resources then 4 else 32
    Int in_vgcall_disk = if in_small_resources then 2 else 40
    String in_vgcall_mem = if in_small_resources then "5" else "64"

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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        CONTIG_ID=($(ls ~{in_child_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_child_name}.//g | sed s/.indel_realigned.bam$//g))
        
        seq 0 $((~{in_vgcall_cores}-1)) | \
        parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/deeptrio/make_examples \
        --mode calling \
        --ref ref.fna \
        --reads_parent1 input_bam_file.paternal.bam \
        --reads_parent2 input_bam_file.maternal.bam \
        --reads input_bam_file.child.bam \
        --examples ./make_examples.tfrecord@~{in_vgcall_cores}.gz \
        --sample_name ~{in_child_name} \
        --sample_name_parent1 ~{in_paternal_name} \
        --sample_name_parent2 ~{in_maternal_name} \
        --gvcf ./gvcf.tfrecord@~{in_vgcall_cores}.gz \
        --min_mapping_quality 1 \
        --pileup_image_height_child 60 \
        --pileup_image_height_parent 40 \
        --regions ${CONTIG_ID} \
        --task {}
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
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:deeptrio-1.1.0"
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
        File? in_model
        Boolean in_small_resources
    }
    
    Int in_vgcall_cores = if in_small_resources then 4 else 8
    Int in_vgcall_disk = if in_small_resources then 2 else 40
    String in_vgcall_mem = if in_small_resources then "5" else "64"
    
    Boolean custom_model = defined(in_model)
    
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
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
        if [ ~{in_sample_type} == "child" ]; then
            EXAMPLES_FILE="make_examples_child.tfrecord@~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_child.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_child.tfrecord@~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/child/model.ckpt"
        elif [ ~{in_sample_type} == "parent1" ]; then
            EXAMPLES_FILE="make_examples_parent1.tfrecord@~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_parent1.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_parent1.tfrecord@~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/parent/model.ckpt"
        elif [ ~{in_sample_type} == "parent2" ]; then
            EXAMPLES_FILE="make_examples_parent2.tfrecord@~{in_vgcall_cores}.gz"
            OUTPUT_FILE="call_variants_output_parent2.tfrecord.gz"
            NONVARIANT_SITE_FILE="gvcf_parent2.tfrecord@~{in_vgcall_cores}.gz"
            DEEPTRIO_MODEL="/opt/models/deeptrio/wgs/parent/model.ckpt"
        fi
        if [ ~{custom_model} == true ]; then
            tar -xzf ~{in_model}
            model_basename=$(basename ~{in_model})
            stripped_basename=${model_basename%%.*}
            core_filename=$(ls $stripped_basename | head -n 1)
            common_filename=${core_filename%.*}
            DEEPTRIO_MODEL="${PWD}/${stripped_basename}/${common_filename}"
        fi
        /opt/deepvariant/bin/call_variants \
        --outfile ${OUTPUT_FILE} \
        --examples ${EXAMPLES_FILE} \
        --checkpoint ${DEEPTRIO_MODEL} && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref ref.fna \
        --infile ${OUTPUT_FILE} \
        --nonvariant_site_tfrecord_path ${NONVARIANT_SITE_FILE} \
        --outfile "~{in_sample_name}_deeptrio.vcf.gz" \
        --gvcf_outfile "~{in_sample_name}_deeptrio.g.vcf.gz"
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deeptrio.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deeptrio.g.vcf.gz"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "google/deepvariant:deeptrio-1.1.0-gpu"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Boolean in_small_resources
    }
    
    Int in_disk = if in_small_resources then 1 else 10
    String in_mem = if in_small_resources then "1" else "1"

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
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
        Boolean in_small_resources
    }
    
    Int in_cores = if in_small_resources then 2 else 2
    Int in_disk = if in_small_resources then 1 else 10
    String in_mem = if in_small_resources then "1" else "10"

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
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: in_vg_container
    }
}


