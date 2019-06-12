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
        Boolean? RUN_VGMPMAP_ALGORITHM
        Boolean VGMPMAP_MODE = select_first([RUN_VGMPMAP_ALGORITHM, true])
        Boolean? RUN_LINEAR_CALLER
        Boolean SURJECT_MODE = select_first([RUN_LINEAR_CALLER, true])
        Boolean? RUN_DRAGEN_CALLER
        Boolean DRAGEN_MODE = select_first([RUN_DRAGEN_CALLER, false])
        Boolean? RUN_GVCF_OUTPUT
        Boolean GVCF_MODE = select_first([RUN_GVCF_OUTPUT, false])
        Boolean? RUN_SV_CALLER
        Boolean SV_CALLER_MODE = select_first([RUN_SV_CALLER, false])
        File MERGED_GAM_FILE
        File MERGED_GAM_GAI_FILE
        String SAMPLE_NAME
        String VG_CONTAINER
        Int READS_PER_CHUNK
        Int CHUNK_BASES
        Int OVERLAP
        File? PATH_LIST_FILE
        File XG_FILE
        File GCSA_FILE
        File GCSA_LCP_FILE
        File? GBWT_FILE
        File? SNARLS_FILE
        Int SPLIT_READ_CORES
        Int SPLIT_READ_DISK
        Int MAP_CORES
        Int MAP_DISK
        Int MAP_MEM
        Int MERGE_GAM_CORES
        Int MERGE_GAM_DISK
        Int MERGE_GAM_MEM
        Int MERGE_GAM_TIME
        Int CHUNK_GAM_CORES
        Int CHUNK_GAM_DISK
        Int CHUNK_GAM_MEM
        Int VGCALL_CORES
        Int VGCALL_DISK
        Int VGCALL_MEM
    }
    
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    call extractPathNames {
        input:
            in_xg_file=XG_FILE,
            in_vg_container=VG_CONTAINER
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractPathNames.output_path_list])
    
    if (!SURJECT_MODE) {
        if (!SV_CALLER_MODE) {
            call chunkAlignmentsByPathNames as chunkAlignmentsDefault {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_merged_sorted_gam=MERGED_GAM_FILE,
                    in_merged_sorted_gam_gai=MERGED_GAM_GAI_FILE,
                    in_xg_file=XG_FILE,
                    in_path_list_file=pipeline_path_list_file,
                    in_chunk_context=50,
                    in_chunk_bases=CHUNK_BASES,
                    in_overlap=OVERLAP,
                    in_vg_container=VG_CONTAINER,
                    in_chunk_gam_cores=CHUNK_GAM_CORES,
                    in_chunk_gam_disk=CHUNK_GAM_DISK,
                    in_chunk_gam_mem=CHUNK_GAM_MEM
            }
            Array[Pair[File,File]] vg_caller_input_files_list_default = zip(chunkAlignmentsDefault.output_vg_chunks, chunkAlignmentsDefault.output_gam_chunks) 
            scatter (vg_caller_input_files in vg_caller_input_files_list_default) { 
                call runVGCaller as VGCallerDefault { 
                    input: 
                        in_sample_name=SAMPLE_NAME, 
                        in_chunk_bed_file=chunkAlignmentsDefault.output_bed_chunk_file, 
                        in_vg_file=vg_caller_input_files.left, 
                        in_gam_file=vg_caller_input_files.right, 
                        in_chunk_bases=CHUNK_BASES, 
                        in_overlap=OVERLAP, 
                        in_vg_container=VG_CONTAINER,
                        in_vgcall_cores=VGCALL_CORES,
                        in_vgcall_disk=VGCALL_DISK,
                        in_vgcall_mem=VGCALL_MEM,
                        in_sv_mode=false
                } 
                call runVCFClipper as VCFClipperDefault { 
                    input: 
                        in_chunk_vcf=VGCallerDefault.output_vcf, 
                        in_chunk_vcf_index=VGCallerDefault.output_vcf_index, 
                        in_chunk_clip_string=VGCallerDefault.clip_string 
                } 
            } 
            call concatClippedVCFChunks as concatVCFChunksDefault { 
                input: 
                    in_sample_name=SAMPLE_NAME, 
                    in_clipped_vcf_chunk_files=VCFClipperDefault.output_clipped_vcf 
            } 
        }
        # Run structural variant graph calling procedure
        # Run recall options and larger vg chunk context for calling structural variants
        if (SV_CALLER_MODE) {
            # Chunk GAM alignments by contigs list
            call chunkAlignmentsByPathNames as chunkAlignmentsSV {
                input:
                in_sample_name=SAMPLE_NAME,
                in_merged_sorted_gam=MERGED_GAM_FILE,
                in_merged_sorted_gam_gai=MERGED_GAM_GAI_FILE,
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_chunk_context=200,
                in_chunk_bases=CHUNK_BASES,
                in_overlap=OVERLAP,
                in_vg_container=VG_CONTAINER,
                in_chunk_gam_cores=CHUNK_GAM_CORES,
                in_chunk_gam_disk=CHUNK_GAM_DISK,
                in_chunk_gam_mem=CHUNK_GAM_MEM
            }
            Array[Pair[File,File]] vg_caller_input_files_list_SV = zip(chunkAlignmentsSV.output_vg_chunks, chunkAlignmentsSV.output_gam_chunks)
            # Run distributed graph-based variant calling
            scatter (vg_caller_input_files in vg_caller_input_files_list_SV) {
                call runVGCaller as VGCallerSV {
                    input:
                        in_sample_name=SAMPLE_NAME,
                        in_chunk_bed_file=chunkAlignmentsSV.output_bed_chunk_file,
                        in_vg_file=vg_caller_input_files.left,
                        in_gam_file=vg_caller_input_files.right,
                        in_chunk_bases=CHUNK_BASES,
                        in_overlap=OVERLAP,
                        in_vg_container=VG_CONTAINER,
                        in_vgcall_cores=VGCALL_CORES,
                        in_vgcall_disk=VGCALL_DISK,
                        in_vgcall_mem=VGCALL_MEM,
                        in_sv_mode=true
                }
                call runVCFClipper as VCFClipperSV {
                    input:
                        in_chunk_vcf=VGCallerSV.output_vcf,
                        in_chunk_vcf_index=VGCallerSV.output_vcf_index,
                        in_chunk_clip_string=VGCallerSV.clip_string
                }
            }
            # Merge distributed variant called VCFs
            call concatClippedVCFChunks as concatVCFChunksSV {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_clipped_vcf_chunk_files=VCFClipperSV.output_clipped_vcf
            }
        }
        # Extract either the normal or structural variant based VCFs and compress them
        File concatted_vcf = select_first([concatVCFChunksDefault.output_merged_vcf, concatVCFChunksSV.output_merged_vcf])
        call bgzipMergedVCF as bgzipVGCalledVCF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_vcf_file=concatted_vcf,
                in_vg_container=VG_CONTAINER
        }
    }
    output {
        File output_vcf = select_first([bgzipVGCalledVCF.output_merged_vcf])
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
        memory: "50 GB"
        disks: "local-disk 50 SSD"
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
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
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
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
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
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
    }
}



