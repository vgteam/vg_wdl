version 1.0

#version draft-2
# Apparently cromwell 36.1 doesn't parse the version line in wdl files if the wdl version is a draft version.
# Cromwell 36.1 wdl parser defaults to draft-2 if a `version` line is not given.
# https://gatkforums.broadinstitute.org/wdl/discussion/15519/escape-characters-in-draft-2-v-1-0

# Working example pipeline which uses VG to map and HaplotypeCaller to call variants.   
# Tested on Googles Cloud Platorm "GCP" and works for VG container "quay.io/vgteam/vg:v1.11.0-215-ge5edc43e-t246-run".  
# WDL pipeline works with JSON input file "vg_pipeline.workingexample.inputs_tiny.json" 
# Steps in pipeline:    
# 1) Split reads into chunks.   
#    500 READS_PER_CHUNK value recommended for HG002 tiny chr 21 dataset    
#      "gs://cmarkell-vg-wdl-dev/HG002_chr21_1.tiny.fastq.gz" and "gs://cmarkell-vg-wdl-dev/HG002_chr21_2.tiny.fastq.gz".   
#    1000000 READS_PER_CHUNK value recommended for HG002 whole chr 21 dataset   
#      "gs://cmarkell-vg-wdl-dev/HG002_chr21_1.fastq.gz" and "gs://cmarkell-vg-wdl-dev/HG002_chr21_2.fastq.gz". 
# 2) Align input paired-end reads to a VG graph using either VG MAP or VG MPMAP algorithms.
#    Set vgPipeline.VGMPMAP_MODE to "true" to run the VG MPMAP algorithm on the read set.
#    Set vgPipeline.VGMPMAP_MODE to "false" to run the VG MAP algorithm on the read set.
# Either run GATK's HaplotypeCaller or use the graph-backed caller "VG Call".
#    Set vgPipeline.SURJECT_MODE to "true" to run HaplotypeCaller or the Dragen modules caller on alignments.
#       3) Surject VG alignments from GAM format to BAM format.
#       4) Merge chunked BAM alignments and preprocess them for GATK compatibility.
#       Either run HaplotypeCaller or use the Dragen module caller on the NIH Biowulf system.
#           Set vgPipeline.RUN_DRAGEN_CALLER to "false" to run GATKs HaplotypeCaller on alignments.
#               5) Run GATK HaplotypeCaller on processed BAM data.
#           Set vgPipeline.RUN_DRAGEN_CALLER to "true" to run the Dragen modules caller on alignments.
#           (Currently only works on NIH Biowulf system)
#               5) Run Dragen variant caller on processed BAM data.
#    Set vgPipeline.SURJECT_MODE to "false" to run VG Call on alignments.   
#       3) Merge graph alignments into a single GAM file
#       4) Chunk GAM file by path names as specified by the "vgPipeline.PATH_LIST_FILE"
#          vgPipeline.PATH_LIST_FILE is a text file that lists path names one per row.
#       5) Run VG Call on GAM files chunked by path name
#          Requires vgPipeline.PATH_LENGTH_FILE which is a text file that lists a tab-delimited set of path names and path lengths
#             one per row.
#       6) Merge VCFs then compress and index the merged VCF.

workflow vgMultiMapCall {
    input {
        Boolean? RUN_LINEAR_CALLER
        Boolean SURJECT_MODE = select_first([RUN_LINEAR_CALLER, true])
        Boolean? RUN_DRAGEN_CALLER
        Boolean DRAGEN_MODE = select_first([RUN_DRAGEN_CALLER, false])
        Boolean? RUN_GVCF_OUTPUT
        Boolean GVCF_MODE = select_first([RUN_GVCF_OUTPUT, false])
        Boolean? RUN_SNPEFF_ANNOTATION
        Boolean SNPEFF_ANNOTATION = select_first([RUN_SNPEFF_ANNOTATION, true])
        File? INPUT_BAM_FILE
        File? INPUT_BAM_FILE_INDEX
        File? INPUT_GAM_FILE
        File? INPUT_GAM_FILE_INDEX
        String SAMPLE_NAME
        String VG_CONTAINER
        Int CHUNK_BASES
        Int OVERLAP
        File PATH_LIST_FILE
        File PATH_LENGTH_FILE
        File XG_FILE
        File REF_FILE
        File REF_INDEX_FILE
        File REF_DICT_FILE
        File SNPEFF_DATABASE
        Int CHUNK_GAM_CORES
        Int CHUNK_GAM_DISK
        Int CHUNK_GAM_MEM
        Int VGCALL_CORES
        Int VGCALL_DISK
        Int VGCALL_MEM
        String DRAGEN_REF_INDEX_NAME
        String UDPBINFO_PATH
        String HELIX_USERNAME
    }
    
    if (SURJECT_MODE) {
        if (!DRAGEN_MODE) {
            if (!GVCF_MODE) {
                call runGATKHaplotypeCaller {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=INPUT_BAM_FILE,
                        in_bam_file_index=INPUT_BAM_FILE_INDEX,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
                call bgzipMergedVCF as bgzipGATKCalledVCF { 
                    input: 
                        in_sample_name=SAMPLE_NAME, 
                        in_merged_vcf_file=runGATKHaplotypeCaller.genotyped_vcf,
                        in_vg_container=VG_CONTAINER 
                }
            }
            if (GVCF_MODE) {
                call runGATKHaplotypeCallerGVCF {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=INPUT_BAM_FILE,
                        in_bam_file_index=INPUT_BAM_FILE_INDEX,
                        in_reference_file=REF_FILE,
                        in_reference_index_file=REF_INDEX_FILE,
                        in_reference_dict_file=REF_DICT_FILE
                }
            }
        }
        if (DRAGEN_MODE) {
            if (!GVCF_MODE) {
                call runDragenCaller {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=INPUT_BAM_FILE,
                        in_bam_file_index=INPUT_BAM_FILE_INDEX,
                        in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                        in_udpbinfo_path=UDPBINFO_PATH,
                        in_helix_username=HELIX_USERNAME
                }
            }
            if (GVCF_MODE) {
                call runDragenCallerGVCF {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_bam_file=INPUT_BAM_FILE,
                        in_bam_file_index=INPUT_BAM_FILE_INDEX,
                        in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                        in_udpbinfo_path=UDPBINFO_PATH,
                        in_helix_username=HELIX_USERNAME
                }
            }
        }
    }
    if (!SURJECT_MODE) {
        call chunkAlignmentsByPathNames {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_sorted_gam=INPUT_GAM_FILE,
                in_merged_sorted_gam_gai=INPUT_GAM_FILE_INDEX,
                in_xg_file=XG_FILE,
                in_path_list_file=PATH_LIST_FILE,
                in_chunk_bases=CHUNK_BASES,
                in_overlap=OVERLAP,
                in_vg_container=VG_CONTAINER,
                in_chunk_gam_cores=CHUNK_GAM_CORES,
                in_chunk_gam_disk=CHUNK_GAM_DISK,
                in_chunk_gam_mem=CHUNK_GAM_MEM
        }
        Array[Pair[File,File]] vg_caller_input_files_list = zip(chunkAlignmentsByPathNames.output_vg_chunks, chunkAlignmentsByPathNames.output_gam_chunks) 
        scatter (vg_caller_input_files in vg_caller_input_files_list) { 
            call runVGCaller { 
                input: 
                    in_sample_name=SAMPLE_NAME, 
                    in_chunk_bed_file=chunkAlignmentsByPathNames.output_bed_chunk_file, 
                    in_vg_file=vg_caller_input_files.left, 
                    in_gam_file=vg_caller_input_files.right, 
                    in_path_length_file=PATH_LENGTH_FILE, 
                    in_chunk_bases=CHUNK_BASES, 
                    in_overlap=OVERLAP, 
                    in_vg_container=VG_CONTAINER,
                    in_vgcall_cores=VGCALL_CORES,
                    in_vgcall_disk=VGCALL_DISK,
                    in_vgcall_mem=VGCALL_MEM
            } 
            call runVCFClipper { 
                input: 
                    in_chunk_vcf=runVGCaller.output_vcf, 
                    in_chunk_vcf_index=runVGCaller.output_vcf_index, 
                    in_chunk_clip_string=runVGCaller.clip_string 
            } 
        } 
        call concatClippedVCFChunks { 
            input: 
                in_sample_name=SAMPLE_NAME, 
                in_clipped_vcf_chunk_files=runVCFClipper.output_clipped_vcf 
        } 
        call bgzipMergedVCF as bgzipVGCalledVCF { 
            input: 
                in_sample_name=SAMPLE_NAME, 
                in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf, 
                in_vg_container=VG_CONTAINER 
        }
    }
    File variantcaller_vcf_output = select_first([bgzipVGCalledVCF.output_merged_vcf, bgzipGATKCalledVCF.output_merged_vcf, runGATKHaplotypeCallerGVCF.rawLikelihoods_gvcf, runDragenCaller.dragen_genotyped_vcf, runDragenCallerGVCF.dragen_genotyped_gvcf])
    if (SNPEFF_ANNOTATION) {
        call normalizeVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bgzip_vcf_file=variantcaller_vcf_output,
        }
        call snpEffAnnotateVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_normalized_vcf_file=normalizeVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
        }
    }
    if (!SNPEFF_ANNOTATION) {
        File final_vcf_output = variantcaller_vcf_output
    }
    output {
        File output_vcf = select_first([snpEffAnnotateVCF.output_snpeff_annotated_vcf, final_vcf_output])
    }
}


task runGATKHaplotypeCaller {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
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

        gatk HaplotypeCaller \
          --native-pair-hmm-threads 32 \
          --reference ${in_reference_file} \
          --input ${in_bam_file} \
          --output ${in_sample_name}.vcf
    }
    output {
        File genotyped_vcf = "${in_sample_name}.vcf"
    }
    runtime {
        memory: 100
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runGATKHaplotypeCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
    }
    
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        gatk HaplotypeCaller \
          --native-pair-hmm-threads 32 \
          -ERC GVCF \
          --reference ${in_reference_file} \
          --input ${in_bam_file} \
          --output ${in_sample_name}.rawLikelihoods.gvcf
    }
    output {
        File rawLikelihoods_gvcf = "${in_sample_name}.rawLikelihoods.gvcf"
    }
    runtime {
        memory: 100
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runDragenCaller {
    input {
        String in_sample_name
        File in_bam_file
        String in_dragen_ref_index_name
        String in_udpbinfo_path
        String in_helix_username
    }
     
    String bam_file_name = basename(in_bam_file)
    
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        mkdir -p /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/ && \
        cp ~{in_bam_file} /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/ && \
        DRAGEN_WORK_DIR_PATH="/staging/~{in_helix_username}/~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${DRAGEN_WORK_DIR_PATH}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} -b /staging/helix/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/~{bam_file_name} --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-sample-name ~{in_sample_name} --intermediate-results-dir ${TMP_DIR} --output-directory ${DRAGEN_WORK_DIR_PATH} --output-file-prefix ~{in_sample_name}_dragen_genotyped" && \
        mkdir /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper && chmod ug+rw -R /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${DRAGEN_WORK_DIR_PATH}/. /staging/helix/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${DRAGEN_WORK_DIR_PATH}/" && \
        mv /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper ~{in_sample_name}_dragen_genotyper && \
        rm -fr /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_genotyped_vcf = "~{in_sample_name}_dragen_genotyper/~{in_sample_name}_dragen_genotyped.vcf.gz"
    }
    runtime {
        memory: 100
    }
}

task runDragenCallerGVCF {
    input {
        String in_sample_name
        File in_bam_file
        String in_dragen_ref_index_name
        String in_udpbinfo_path
        String in_helix_username
    }
    
    String bam_file_name = basename(in_bam_file)
    
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        mkdir -p /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/ && \
        cp ~{in_bam_file} /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/ && \
        DRAGEN_WORK_DIR_PATH="/staging/~{in_helix_username}/~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${DRAGEN_WORK_DIR_PATH}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} -b /staging/helix/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/~{bam_file_name} --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true --pair-by-name=true --vc-emit-ref-confidence GVCF --vc-sample-name ~{in_sample_name} --intermediate-results-dir ${TMP_DIR} --output-directory ${DRAGEN_WORK_DIR_PATH} --output-file-prefix ~{in_sample_name}_dragen_genotyped" && \
        mkdir /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper && chmod ug+rw -R /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${DRAGEN_WORK_DIR_PATH}/. /staging/helix/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${DRAGEN_WORK_DIR_PATH}/" && \
        mv /data/~{in_udpbinfo_path}/~{in_sample_name}_dragen_genotyper ~{in_sample_name}_dragen_genotyper && \
        rm -fr /data/~{in_udpbinfo_path}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_genotyped_gvcf = "~{in_sample_name}_dragen_genotyper/~{in_sample_name}_dragen_genotyped.gvcf.gz"
    }
    runtime {
        memory: 100
    }
}


task chunkAlignmentsByPathNames {
    input {
        String in_sample_name
        File in_merged_sorted_gam
        File in_merged_sorted_gam_gai
        File in_xg_file
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
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
    
        vg chunk \
            -x ${in_xg_file} \
            -a ${in_merged_sorted_gam} \
            -c 50 \
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
        memory: in_chunk_gam_mem
        cpu: in_chunk_gam_cores
        disks: in_chunk_gam_disk
        docker: in_vg_container
    }   
}

task runVGCaller {
    input {
        String in_sample_name
        File in_chunk_bed_file
        File in_vg_file
        File in_gam_file
        File in_path_length_file
        Int in_chunk_bases
        Int in_overlap
        String in_vg_container
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
    }
    
    String chunk_tag = basename(in_vg_file, ".vg")
    command <<< 
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        PATH_NAME=$(echo ~{chunk_tag} | cut -f 4 -d '_')
        OFFSET=""
        BED_CHUNK_LINE_RECORD=""
        FIRST_Q=""
        LAST_Q=""
        while IFS=$'\t' read -ra bed_chunk_line; do
            bed_line_chunk_tag=$(echo ${bed_chunk_line[3]} | cut -f 1 -d '.')
            if [ "~{chunk_tag}" = "${bed_line_chunk_tag}" ]; then
                OFFSET=$(( bed_chunk_line[1] + 1 ))
                BED_CHUNK_LINE_RECORD="${bed_chunk_line[@]//$'\t'/ }"
            fi
            bed_line_path_name=$(echo ${bed_line_chunk_tag} | cut -f 4 -d '_')
            if [ "${PATH_NAME}" = "${bed_line_path_name}" ]; then
                if [ "${FIRST_Q}" = "" ]; then
                    FIRST_Q="${bed_chunk_line[@]//$'\t'/ }"
                fi
                LAST_Q="${bed_chunk_line[@]//$'\t'/ }"
            fi
        done < ~{in_chunk_bed_file}

        CHR_LENGTH=""
        while IFS=$'\t' read -ra path_length_line; do
            path_length_line_name=${path_length_line[0]}
            if [ "${PATH_NAME}" = "${path_length_line_name}" ]; then
                CHR_LENGTH=${path_length_line[1]}
            fi
        done < ~{in_path_length_file}

        vg index \
            ~{in_vg_file} \
            -x ~{chunk_tag}.xg \
            -t ~{in_vgcall_cores} && \
        vg filter \
            ~{in_gam_file} \
            -t 1 -r 0.9 -fu -s 1000 -m 1 -q 15 -D 999 \
            -x ~{chunk_tag}.xg > ~{chunk_tag}.filtered.gam && \
        vg augment \
            ~{in_vg_file} \
            ~{chunk_tag}.filtered.gam \
            -t ~{in_vgcall_cores} -q 10 -a pileup \
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
            -o ${OFFSET} > ~{chunk_tag}.vcf && \
        head -10000 ~{chunk_tag}.vcf | grep "^#" >> ~{chunk_tag}.sorted.vcf && \
        if [ "~(cat ~{chunk_tag}.vcf | grep -v '^#')" ]; then
            cat ~{chunk_tag}.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ~{chunk_tag}.sorted.vcf
        fi && \
        bgzip ~{chunk_tag}.sorted.vcf && \
        tabix -f -p vcf ~{chunk_tag}.sorted.vcf.gz && \

        # Compile clipping string
        chunk_tag_index=$(echo ~{chunk_tag} | cut -f 3 -d '_') && \
        clipped_chunk_offset=$(( ${chunk_tag_index} * ~{in_chunk_bases} - ${chunk_tag_index} * ~{in_overlap} )) && \
        if [ "${BED_CHUNK_LINE_RECORD}" = "${FIRST_Q}" ]; then
          START=$(( clipped_chunk_offset + 1 ))
        else
          START=$(( clipped_chunk_offset + ~{in_overlap}/2 ))
        fi
        if [ "${BED_CHUNK_LINE_RECORD}" = "${LAST_Q}" ]; then
          END=$(( clipped_chunk_offset + ~{in_chunk_bases} - 1 ))
        else
          END=$(( clipped_chunk_offset + ~{in_chunk_bases} - (~{in_overlap}/2) ))
        fi
        echo "${PATH_NAME}:${START}-${END}" > ~{chunk_tag}_clip_string.txt
    >>>
    output {
        File output_vcf = "~{chunk_tag}.sorted.vcf.gz"
        File output_vcf_index = "~{chunk_tag}.sorted.vcf.gz.tbi"
        String clip_string = read_lines("~{chunk_tag}_clip_string.txt")[0]
    }
    runtime {
        memory: in_vgcall_mem
        cpu: in_vgcall_cores
        disks: in_vgcall_disk
        docker: in_vg_container
    }
}

task runVCFClipper {
    input {
        File in_chunk_vcf
        File in_chunk_vcf_index
        String in_chunk_clip_string
    }
    
    String chunk_tag = basename(in_chunk_vcf, ".sorted.vcf.gz")
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        bcftools view ${in_chunk_vcf} -t ${in_chunk_clip_string} > ${chunk_tag}.clipped.vcf
    }
    output {
        File output_clipped_vcf = "${chunk_tag}.clipped.vcf"
    }
    runtime {
        memory: 50
        disks: 100
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
    }
    
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        bcftools concat ${sep=" " in_clipped_vcf_chunk_files} > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        memory: 50
        disks: 100
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
    }
    
    # TODO:
    #   If GVCF in in_merged_vcf_file then output_vcf_extension="gvcf" else output_vcf_extension="vcf"
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
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
        memory: 50
        disks: 100
        docker: in_vg_container
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
    }
    
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        bcftools norm -m-both --threads 16 -o ${in_sample_name}.unrolled.vcf ${in_bgzip_vcf_file}
    }
    output {
        File output_normalized_vcf = "${in_sample_name}.unrolled.vcf"
    }
    runtime {
        cpu: 16
        memory: 50
        disks: 100
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File in_snpeff_database
    }
    
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        # cause a bash script to exit immediately when a command fails
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -e -u -o pipefail
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
        cpu: 16
        memory: 50
        disks: 100
        docker: "quay.io/biocontainers/snpeff:4.3.1t--2"
    }
}

