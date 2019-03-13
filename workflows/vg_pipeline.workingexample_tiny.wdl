#version draft-2
# Apparently cromwell 36.1 doesn't parse the version line in wdl files if the wdl version is a draft version.
# Cromwell 36.1 wdl parser defaults to draft-2 if a `version` line is not given.

# Working example pipeline which uses VG to map and HaplotypeCaller to call variants.	
# Tested on Googles Cloud Platorm "GCP" and works for VG container "quay.io/vgteam/vg:v1.11.0-215-ge5edc43e-t246-run".	
# WDL pipeline works with JSON input file "vg_pipeline.workingexample_tiny.inputs.json"	
# Steps in pipeline: 	
# 1) Split reads into chunks.	
#    500 READS_PER_CHUNK value recommended for HG002 tiny chr 21 dataset	
#      "gs://cmarkell-vg-wdl-dev/HG002_chr21_1.tiny.fastq.gz" and "gs://cmarkell-vg-wdl-dev/HG002_chr21_2.tiny.fastq.gz".	
#    1000000 READS_PER_CHUNK value recommended for HG002 whole chr 21 dataset	
#      "gs://cmarkell-vg-wdl-dev/HG002_chr21_1.fastq.gz" and "gs://cmarkell-vg-wdl-dev/HG002_chr21_2.fastq.gz".	
# 2) Align input paired-end reads to a VG graph using VG MPMAP algorithm.	
# 3) Surject VG alignments from GAM format to BAM format.	
# 4) Merge chunked BAM alignments and preprocess them for GATK compatibility.	
# 5) Run GATK HaplotypeCaller on processed BAM data.

workflow vgPipeline {
    File INPUT_READ_FILE_1
    File INPUT_READ_FILE_2
    String SAMPLE_NAME
    String VG_CONTAINER
    Int READS_PER_CHUNK
    Int CHUNK_BASES
    Int OVERLAP
    File PATH_LIST_FILE
    File PATH_LENGTH_FILE
    File VG_FILE
    File XG_FILE
    File GCSA_FILE
    File GCSA_LCP_FILE
    File GBWT_FILE
    File REF_FILE
    File REF_INDEX_FILE
    File REF_DICT_FILE
    Boolean VGMPMAP_MODE
    Boolean SURJECT_MODE

    call splitReads as firstReadPair {
    input:
        in_read_file=INPUT_READ_FILE_1,
        in_pair_id="1",
        in_vg_container=VG_CONTAINER,
        in_reads_per_chunk=READS_PER_CHUNK
    }
    call splitReads as secondReadPair {
    input:
        in_read_file=INPUT_READ_FILE_2,
        in_pair_id="2",
        in_vg_container=VG_CONTAINER,
        in_reads_per_chunk=READS_PER_CHUNK
    }
    Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
    scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
        if (VGMPMAP_MODE) {
            call runVGMPMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_sample_name=SAMPLE_NAME
            }
        } 
        if (!VGMPMAP_MODE) {
            call runVGMAP {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_xg_file=XG_FILE,
                    in_gcsa_file=GCSA_FILE,
                    in_gcsa_lcp_file=GCSA_LCP_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_sample_name=SAMPLE_NAME
            }
        }
        File vg_map_algorithm_chunk_gam_output = select_first([runVGMAP.chunk_gam_file, runVGMPMAP.chunk_gam_file])
        if (SURJECT_MODE) {
            call runSurject {
                input:
                    in_gam_chunk_file=vg_map_algorithm_chunk_gam_output,
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_sample_name=SAMPLE_NAME
            }
        }
    }
    if (SURJECT_MODE) {
        Array[File?] alignment_chunk_bam_files_maybes = runSurject.chunk_bam_file
        Array[File] alignment_chunk_bam_files_valid = select_all(alignment_chunk_bam_files_maybes)
        call mergeAlignmentBAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=alignment_chunk_bam_files_valid
        }
        call runGATKHaplotypeCaller {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
                in_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
    }
    if (!SURJECT_MODE) {
        call mergeAlignmentGAMChunks {
            input:
                in_sample_name=SAMPLE_NAME,
                in_vg_container=VG_CONTAINER,
                in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output
        }
        call chunkAlignmentsByPathNames {
            input:
                in_sample_name=SAMPLE_NAME,
                in_merged_sorted_gam=mergeAlignmentGAMChunks.merged_sorted_gam_file,
                in_merged_sorted_gam_gai=mergeAlignmentGAMChunks.merged_sorted_gam_gai_file,
                in_vg_file=VG_FILE,
                in_xg_file=XG_FILE,
                in_path_list_file=PATH_LIST_FILE,
                in_chunk_bases=CHUNK_BASES,
                in_overlap=OVERLAP,
                in_vg_container=VG_CONTAINER
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
                    in_vg_container=VG_CONTAINER 
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
        call postprocessMergedVCF { 
            input: 
                in_sample_name=SAMPLE_NAME, 
                in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf, 
                in_vg_container=VG_CONTAINER 
        }
    }
}

task splitReads {
    File in_read_file
    String in_pair_id
    String in_vg_container
    Int in_reads_per_chunk

    String dollar = "$" 
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

        let CHUNK_LINES=${in_reads_per_chunk}*4
        gzip -cd ${in_read_file} | split -l ${dollar}{CHUNK_LINES} --filter='pigz -p 32 > ${dollar}{FILE}.fq.gz' - fq_chunk_${in_pair_id}.part.
    >>> 
    output {
        Array[File] output_read_chunks = glob("fq_chunk_${in_pair_id}.part.*")
    }
    runtime {
        cpu: "32"
        docker: in_vg_container
    }
}

task runVGMAP {
    File in_left_read_pair_chunk_file
    File in_right_read_pair_chunk_file
    File in_xg_file
    File in_gcsa_file
    File in_gcsa_lcp_file
    File in_gbwt_file
    String in_vg_container
    String in_sample_name

    String dollar = "$" 
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

        READ_CHUNK_ID=($(ls ${in_left_read_pair_chunk_file} | awk -F'.' '{print ${dollar}3}'))
        vg map \
          -x ${in_xg_file} \
          -g ${in_gcsa_file} \
          -f ${in_left_read_pair_chunk_file} -f ${in_right_read_pair_chunk_file} \
          -t 32 > ${in_sample_name}.${dollar}{READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}

task runVGMPMAP {
    File in_left_read_pair_chunk_file
    File in_right_read_pair_chunk_file
    File in_xg_file
    File in_gcsa_file
    File in_gcsa_lcp_file
    File in_gbwt_file
    String in_vg_container
    String in_sample_name

    String dollar = "$" 
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

        READ_CHUNK_ID=($(ls ${in_left_read_pair_chunk_file} | awk -F'.' '{print ${dollar}3}'))
        vg mpmap \
          -S \
          -f ${in_left_read_pair_chunk_file} -f ${in_right_read_pair_chunk_file} \
          -x ${in_xg_file} \
          -g ${in_gcsa_file} \
          --gbwt-name ${in_gbwt_file} \
          --recombination-penalty 5.0 -t 32 > ${in_sample_name}.${dollar}{READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}

task runSurject {
    File in_gam_chunk_file
    File in_xg_file
    String in_vg_container
    String in_sample_name

    String dollar = "$"
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

        READ_CHUNK_ID=($(ls ${in_gam_chunk_file} | awk -F'.' '{print ${dollar}2}'))
        vg surject \
          -i \
          -x ${in_xg_file} \
          -t 32 \
          -b ${in_gam_chunk_file} > ${in_sample_name}.${dollar}{READ_CHUNK_ID}.bam
    >>>
    output {
        File chunk_bam_file = glob("*.bam")[0]
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}

task mergeAlignmentBAMChunks {
    String in_sample_name
    Array[File] in_alignment_bam_chunk_files

    String dollar = "$"
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
          -f --threads 32 \
          ${in_sample_name}_merged.bam \
          ${sep=" " in_alignment_bam_chunk_files} \
        && samtools sort \
          --threads 32 \
          ${in_sample_name}_merged.bam \
          -n \
          -O BAM \
          -o ${in_sample_name}_merged.namesorted.bam \
        && samtools fixmate \
          -O BAM \
          ${in_sample_name}_merged.namesorted.bam \
          ${in_sample_name}_merged.namesorted.fixmate.bam \
        && samtools sort \
          --threads 32 \
          ${in_sample_name}_merged.namesorted.fixmate.bam \
          -O BAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.bam \
        && samtools addreplacerg \
          -O BAM \
          -r ID:1 -r LB:lib1 -r SM:${in_sample_name} -r PL:illumina -r PU:unit1 \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
          ${in_sample_name}_merged.fixmate.positionsorted.bam \
        && samtools view \
          -@ 32 \
          -h -O SAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.sam \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
        && samtools view \
          -@ 32 \
          -h -O BAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.sam \
        && samtools index \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.bam
    >>>
    output {
        File merged_bam_file = "${in_sample_name}_merged.fixmate.positionsorted.rg.bam"
        File merged_bam_file_index = "${in_sample_name}_merged.fixmate.positionsorted.rg.bam.bai"
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task runGATKHaplotypeCaller {
    String in_sample_name
    File in_bam_file
    File in_bam_file_index
    File in_reference_file
    File in_reference_index_file
    File in_reference_dict_file

    String dollar = "$"
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

        gatk HaplotypeCaller \
          --reference ${in_reference_file} \
          --input ${in_bam_file} \
          --output ${in_sample_name}.vcf
    >>>
    output {
        File genotyped_vcf = "${in_sample_name}.vcf"
    }
    runtime {
        memory: "100G"
        cpu: "32"
        docker: "broadinstitute/gatk:latest"
    }
}

task mergeAlignmentGAMChunks {
    String in_sample_name
    String in_vg_container
    Array[File] in_alignment_gam_chunk_files

    String dollar = "$"
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
        
        cat ${sep=" " in_alignment_gam_chunk_files} > ${in_sample_name}_merged.gam \
        && vg gamsort \
            ${in_sample_name}_merged.gam \
            -i ${in_sample_name}_merged.sorted.gam.gai \
            -t 32 > ${in_sample_name}_merged.sorted.gam
    >>>
    output {
        File merged_sorted_gam_file = "${in_sample_name}_merged.sorted.gam"
        File merged_sorted_gam_gai_file = "${in_sample_name}_merged.sorted.gam.gai"
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}

task chunkAlignmentsByPathNames {
    String in_sample_name
    File in_merged_sorted_gam
    File in_merged_sorted_gam_gai
    File in_vg_file
    File in_xg_file
    File in_path_list_file
    Int in_chunk_bases
    Int in_overlap
    String in_vg_container

    String dollar = "$" 
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
    
        vg chunk \
            -x ${in_xg_file} \
            -a ${in_merged_sorted_gam} \
            -c 50 \
            -P ${in_path_list_file} \
            -g \
            -s ${in_chunk_bases} -o ${in_overlap} \
            -b call_chunk \
            -t 32 \
            -E output_bed_chunks.bed -f
    >>> 
    output {
        Array[File] output_vg_chunks = glob("call_chunk*.vg")
        Array[File] output_gam_chunks = glob("call_chunk*.gam")
        File output_bed_chunk_file = "output_bed_chunks.bed"
    }   
    runtime {
        memory: "400G"
        cpu: "32"
        disks: "local-disk 200 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }   
}

task runVGCaller {
    String in_sample_name
    File in_chunk_bed_file
    File in_vg_file
    File in_gam_file
    File in_path_length_file
    Int in_chunk_bases
    Int in_overlap
    String in_vg_container

    String chunk_tag = basename(in_vg_file, ".vg")
    String dollar = "$" 
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

        PATH_NAME=${dollar}(echo ${chunk_tag} | cut -f 4 -d '_')
        OFFSET=""
        BED_CHUNK_LINE_RECORD=""
        FIRST_Q=""
        LAST_Q=""
        while IFS=${dollar}'\t' read -ra bed_chunk_line; do
            bed_line_chunk_tag=${dollar}(echo ${dollar}{bed_chunk_line[3]} | cut -f 1 -d '.')
            if [ "${chunk_tag}" = "${dollar}{bed_line_chunk_tag}" ]; then
                OFFSET=${dollar}(( bed_chunk_line[1] + 1 ))
                BED_CHUNK_LINE_RECORD="${dollar}{bed_chunk_line[@]//${dollar}'\t'/ }"
            fi
            bed_line_path_name=${dollar}(echo ${dollar}{bed_line_chunk_tag} | cut -f 4 -d '_')
            if [ "${dollar}{PATH_NAME}" = "${dollar}{bed_line_path_name}" ]; then
                if [ "${dollar}{FIRST_Q}" = "" ]; then
                    FIRST_Q="${dollar}{bed_chunk_line[@]//${dollar}'\t'/ }"
                fi
                LAST_Q="${dollar}{bed_chunk_line[@]//${dollar}'\t'/ }"
            fi
        done < ${in_chunk_bed_file}

        CHR_LENGTH=""
        while IFS=${dollar}'\t' read -ra path_length_line; do
            path_length_line_name=${dollar}{path_length_line[0]}
            if [ "${dollar}{PATH_NAME}" = "${dollar}{path_length_line_name}" ]; then
                CHR_LENGTH=${dollar}{path_length_line[1]}
            fi
        done < ${in_path_length_file}

        vg index \
            ${in_vg_file} \
            -x ${chunk_tag}.xg \
            -t 32 && \
        vg filter \
            ${in_gam_file} \
            -t 1 -r 0.9 -fu -s 1000 -m 1 -q 15 -D 999 \
            -x ${chunk_tag}.xg > ${chunk_tag}.filtered.gam && \
        vg augment \
            ${in_vg_file} \
            ${chunk_tag}.filtered.gam \
            -t 32 -q 10 -a pileup \
            -Z ${chunk_tag}.trans \
            -S ${chunk_tag}.support > ${chunk_tag}.aug.vg && \
        vg call \
            ${chunk_tag}.aug.vg \
            -t 32 \
            -S ${in_sample_name} \
            -z ${chunk_tag}.trans \
            -s ${chunk_tag}.support \
            -b ${in_vg_file} \
            -r ${dollar}{PATH_NAME} \
            -c ${dollar}{PATH_NAME} \
            -l ${dollar}{CHR_LENGTH} \
            -o ${dollar}{OFFSET} > ${chunk_tag}.vcf && \
        head -10000 ${chunk_tag}.vcf | grep "^#" >> ${chunk_tag}.sorted.vcf && \
        if [ "$(cat ${chunk_tag}.vcf | grep -v '^#')" ]; then
            cat ${chunk_tag}.vcf | grep -v "^#" | sort -k1,1d -k2,2n >> ${chunk_tag}.sorted.vcf
        fi && \
        bgzip ${chunk_tag}.sorted.vcf && \
        tabix -f -p vcf ${chunk_tag}.sorted.vcf.gz && \

        # Compile clipping string
        chunk_tag_index=${dollar}(echo ${chunk_tag} | cut -f 3 -d '_') && \
        clipped_chunk_offset=${dollar}(( ${dollar}{chunk_tag_index} * ${in_chunk_bases} - ${dollar}{chunk_tag_index} * ${in_overlap} )) && \
        if [ "${dollar}{BED_CHUNK_LINE_RECORD}" = "${dollar}{FIRST_Q}" ]; then
          START=${dollar}(( clipped_chunk_offset + 1 ))
        else
          START=${dollar}(( clipped_chunk_offset + ${in_overlap}/2 ))
        fi
        if [ "${dollar}{BED_CHUNK_LINE_RECORD}" = "${dollar}{LAST_Q}" ]; then
          END=${dollar}(( clipped_chunk_offset + ${in_chunk_bases} - 1 ))
        else
          END=${dollar}(( clipped_chunk_offset + ${in_chunk_bases} - (${in_overlap}/2) ))
        fi
        echo "${dollar}{PATH_NAME}:${dollar}{START}-${dollar}{END}" > ${chunk_tag}_clip_string.txt
    >>>
    output {
        File output_vcf = "${chunk_tag}.sorted.vcf.gz"
        File output_vcf_index = "${chunk_tag}.sorted.vcf.gz.tbi"
        String clip_string = read_lines("${chunk_tag}_clip_string.txt")[0]
    }
    runtime {
        memory: "100G"
        cpu: "32"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}

task runVCFClipper {
    File in_chunk_vcf
    File in_chunk_vcf_index
    String in_chunk_clip_string

    String chunk_tag = basename(in_chunk_vcf, ".sorted.vcf.gz")
    String dollar = "$"
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

        bcftools view ${in_chunk_vcf} -t ${in_chunk_clip_string} > ${chunk_tag}.clipped.vcf
    >>>
    output {
        File output_clipped_vcf = "${chunk_tag}.clipped.vcf"
    }
    runtime {
        memory: "50G"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
        zones: "us-west1-b"
    }
}

task concatClippedVCFChunks {
    String in_sample_name
    Array[File] in_clipped_vcf_chunk_files

    String dollar = "$"
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

        bcftools concat ${sep=" " in_clipped_vcf_chunk_files} > ${in_sample_name}_merged.vcf
    >>>
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        memory: "50G"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
        zones: "us-west1-b"
    }
}

task postprocessMergedVCF {
    String in_sample_name
    File in_merged_vcf_file
    String in_vg_container

    String dollar = "$"
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

        bgzip -c ${in_merged_vcf_file} > ${in_sample_name}_merged.vcf.gz && \
        tabix -f -p vcf ${in_sample_name}_merged.vcf.gz
    >>>
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}_merged.vcf.gz.tbi"
    }
    runtime {
        memory: "50G"
        disks: "local-disk 100 SSD"
        docker: in_vg_container
        zones: "us-west1-b"
    }
}


