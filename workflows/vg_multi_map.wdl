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
        File INPUT_READ_FILE_1
        File INPUT_READ_FILE_2
        String SAMPLE_NAME
        String VG_CONTAINER
        Int READS_PER_CHUNK
        File XG_FILE
        File GCSA_FILE
        File GCSA_LCP_FILE
        File? GBWT_FILE
        File? SNARLS_FILE
        File REF_FILE
        File REF_INDEX_FILE
        File REF_DICT_FILE
        Int SPLIT_READ_CORES
        Int SPLIT_READ_DISK
        Int MAP_CORES
        Int MAP_DISK
        Int MAP_MEM
        Int MERGE_GAM_CORES
        Int MERGE_GAM_DISK
        Int MERGE_GAM_MEM
        Int MERGE_GAM_TIME
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

    # Distribute vg mapping opperation over each chunked read pair
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
                    in_snarls_file=SNARLS_FILE,
                    in_sample_name=SAMPLE_NAME,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
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
                    in_sample_name=SAMPLE_NAME,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        
        # Surject GAM alignment files to BAM if SURJECT_MODE set to true
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
        Array[File] alignment_chunk_bam_files = select_all(runSurject.chunk_bam_file)
        if (!VGMPMAP_MODE) {
            call mergeAlignmentBAMChunksVGMAP {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_alignment_bam_chunk_files=alignment_chunk_bam_files,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE
            }
        }
        if (VGMPMAP_MODE) {
            call mergeAlignmentBAMChunksVGMPMAP {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_alignment_bam_chunk_files=alignment_chunk_bam_files,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE
            }
        }
        File merged_bam_file_output = select_first([mergeAlignmentBAMChunksVGMAP.merged_bam_file, mergeAlignmentBAMChunksVGMPMAP.merged_bam_file])
        File merged_bam_file_index_output = select_first([mergeAlignmentBAMChunksVGMAP.merged_bam_file_index, mergeAlignmentBAMChunksVGMPMAP.merged_bam_file_index])
        call runPICARD {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=merged_bam_file_output,
                in_bam_file_index=merged_bam_file_index_output,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        call runGATKIndelRealigner {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=runPICARD.mark_dupped_reordered_bam,
                in_bam_file_index=runPICARD.mark_dupped_reordered_bam_index,
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
                in_alignment_gam_chunk_files=vg_map_algorithm_chunk_gam_output,
                in_merge_gam_cores=MERGE_GAM_CORES,
                in_merge_gam_disk=MERGE_GAM_DISK,
                in_merge_gam_mem=MERGE_GAM_MEM,
                in_merge_gam_time=MERGE_GAM_TIME
        }
    }
    output {
        File? output_bam = runGATKIndelRealigner.indel_realigned_bam
        File? output_bam_index = runGATKIndelRealigner.indel_realigned_bam_index
        File? output_gam = mergeAlignmentGAMChunks.merged_sorted_gam_file
        File? output_gam_index = mergeAlignmentGAMChunks.merged_sorted_gam_gai_file
    }
}

# Now works for WDL parser version 1.0
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
        cpu: in_split_read_cores
        disks: in_split_read_disk
        docker: in_vg_container
    }
}

task runVGMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gcsa_file
        File in_gcsa_lcp_file
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        Int in_map_mem
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
        vg map \
          -x ~{in_xg_file} \
          -g ~{in_gcsa_file} \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: in_map_mem
        cpu: in_map_cores
        disks: in_map_disk
        docker: in_vg_container
    }
}

task runVGMPMAP {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gcsa_file
        File in_gcsa_lcp_file
        File? in_gbwt_file
        File? in_snarls_file
        String in_vg_container
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        Int in_map_mem
    }
    
    Boolean gbwt_options = defined(in_gbwt_file) && defined(in_snarls_file)
    
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
        if [ ~{gbwt_options} == true ]; then
          GBWT_OPTION_STRING="--gbwt-name ~{in_gbwt_file} -s ~{in_snarls_file}"
        else
          GBWT_OPTION_STRING=""
        fi
        vg mpmap \
          -S \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -g ~{in_gcsa_file} \
          ${GBWT_OPTION_STRING} \
          --read-group "ID:1\tLB:lib1\tSM:~{in_sample_name}\tPL:illumina\tPU:unit1" \
          --sample "~{in_sample_name}" \
          --recombination-penalty 5.0 -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*.gam")[0]
    }
    runtime {
        memory: in_map_mem
        cpu: in_map_cores
        disks: in_map_disk
        docker: in_vg_container
    }
}

task runSurject {
    input {
        File in_gam_chunk_file
        File in_xg_file
        String in_vg_container
        String in_sample_name
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

        READ_CHUNK_ID=($(ls ~{in_gam_chunk_file} | awk -F'.' '{print $(NF-1)}'))
        vg surject \
          -i \
          -x ~{in_xg_file} \
          -t 32 \
          -b ~{in_gam_chunk_file} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
    >>>
    output {
        File chunk_bam_file = glob("*.bam")[0]
    }
    runtime {
        memory: 100
        cpu: 32
        disks: 100
        docker: in_vg_container
    }
}

task mergeAlignmentBAMChunksVGMAP {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
        && rm -f ${in_sample_name}_merged.bam \
        && samtools fixmate \
          -O BAM \
          ${in_sample_name}_merged.namesorted.bam \
          ${in_sample_name}_merged.namesorted.fixmate.bam \
        && rm -f ${in_sample_name}_merged.namesorted.bam \
        && samtools sort \
          --threads 32 \
          ${in_sample_name}_merged.namesorted.fixmate.bam \
          -O BAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.bam \
        && rm -f ${in_sample_name}_merged.namesorted.fixmate.bam \
        && samtools addreplacerg \
          -O BAM \
          -r ID:1 -r LB:lib1 -r SM:${in_sample_name} -r PL:illumina -r PU:unit1 \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
          ${in_sample_name}_merged.fixmate.positionsorted.bam \
        && rm -f ${in_sample_name}_merged.fixmate.positionsorted.bam \
        && samtools view \
          -@ 32 \
          -h -O SAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.sam \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
        && rm -f ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
        && samtools view \
          -@ 32 \
          -h -O BAM \
          -o ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.sam \
        && rm -f ${in_sample_name}_merged.fixmate.positionsorted.rg.sam \
        && samtools calmd \
          -b \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
          ${in_reference_file} \
          > ${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.bam \
        && rm -f ${in_sample_name}_merged.fixmate.positionsorted.rg.bam \
        && samtools index \
          ${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.bam
    }
    output {
        File merged_bam_file = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.bam"
        File merged_bam_file_index = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.bam.bai"
    }
    runtime {
        memory: 100
        cpu: 32
        disks: 100
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task mergeAlignmentBAMChunksVGMPMAP {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
        samtools merge \
          -f --threads 32 \
          ${in_sample_name}_merged.bam \
          ${sep=" " in_alignment_bam_chunk_files} \
        && samtools sort \
          --threads 32 \
          ${in_sample_name}_merged.bam \
          -O BAM \
          -o ${in_sample_name}_merged.positionsorted.bam \
        && samtools calmd \
          -b \
          ${in_sample_name}_merged.positionsorted.bam \
          ${in_reference_file} \
          > ${in_sample_name}_merged.positionsorted.mdtag.bam \
        && rm -f ${in_sample_name}_merged.positionsorted.bam \
        && samtools index \
          ${in_sample_name}_merged.positionsorted.mdtag.bam
    }
    output {
        File merged_bam_file = "${in_sample_name}_merged.positionsorted.mdtag.bam"
        File merged_bam_file_index = "${in_sample_name}_merged.positionsorted.mdtag.bam.bai"
    }
    runtime {
        memory: 100
        cpu: 32
        disks: 100
        docker: "biocontainers/samtools:v1.3_cv3"
    }
}

task runPICARD {
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

        java -Xmx4g -XX:ParallelGCThreads=32 -jar /usr/picard/picard.jar MarkDuplicates \
          VALIDATION_STRINGENCY=LENIENT \
          I=${in_bam_file} \
          O=${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.bam \
          M=marked_dup_metrics.txt 2> mark_dup_stderr.txt \
        && java -Xmx20g -XX:ParallelGCThreads=32 -jar /usr/picard/picard.jar ReorderSam \
            VALIDATION_STRINGENCY=LENIENT \
            REFERENCE=${in_reference_file} \
            INPUT=${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.bam \
            OUTPUT=${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.bam \
        && java -Xmx20g -XX:ParallelGCThreads=32 -jar /usr/picard/picard.jar BuildBamIndex \
            VALIDATION_STRINGENCY=LENIENT \
            I=${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.bam \
            O=${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.bam.bai
    }
    output {
        File mark_dupped_reordered_bam = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.bam"
        File mark_dupped_reordered_bam_index = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.bam.bai"
    }
    runtime {
        memory: 100
        cpu: 32
        docker: "broadinstitute/picard:latest"
    }
}


task runGATKIndelRealigner {
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

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          -R ${in_reference_file} \
          -I ${in_bam_file} \
          --out forIndelRealigner.intervals \
        && java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
          -R ${in_reference_file} \
          --targetIntervals forIndelRealigner.intervals \
          -I ${in_bam_file} \
          --out ${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.indel_realigned.bam
    }
    output {
        File indel_realigned_bam = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.indel_realigned.bam"
        File indel_realigned_bam_index = "${in_sample_name}_merged.fixmate.positionsorted.rg.mdtag.dupmarked.reordered.indel_realigned.bai"
    }
    runtime {
        memory: 100
        cpu: 32
        docker: "broadinstitute/gatk3:3.8-1"
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
        memory: in_merge_gam_mem
        cpu: in_merge_gam_cores
        disks: in_merge_gam_mem
        time: in_merge_gam_time
        docker: in_vg_container
    }
}



