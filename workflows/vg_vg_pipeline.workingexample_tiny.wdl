workflow vgPipeline {
  File INPUT_READ_FILE_1
  File INPUT_READ_FILE_2
  String SAMPLE_NAME
  String VG_CONTAINER
  Int READS_PER_CHUNK
  File XG_FILE
  File GCSA_FILE
  File GCSA_LCP_FILE
  File GBWT_FILE
  File REF_FILE
  File REF_INDEX_FILE
  File REF_DICT_FILE
  
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
    call runMapper {
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
    call runSurject {
      input:
        in_gam_chunk_file=runMapper.chunk_gam_file,
        in_xg_file=XG_FILE,
        in_vg_container=VG_CONTAINER,
        in_sample_name=SAMPLE_NAME
    }
  }
  call mergeAlignmentChunks {
    input:
      in_sample_name=SAMPLE_NAME,
      in_alignment_chunk_files=runSurject.chunk_bam_file,
  }
  call runGenotyper {
    input:
      in_sample_name=SAMPLE_NAME,
      in_bam_file=mergeAlignmentChunks.merged_bam_file,
      in_bam_file_index=mergeAlignmentChunks.merged_bam_file_index,
      in_reference_file=REF_FILE,
      in_reference_index_file=REF_INDEX_FILE,
      in_reference_dict_file=REF_DICT_FILE
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
    docker: in_vg_container
  }
}

task runMapper {
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
      --recombination-penalty 5.0 -t 32 > ${in_sample_name}.${dollar}{READ_CHUNK_ID}.gam \
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

task mergeAlignmentChunks {
  String in_sample_name
  Array[File] in_alignment_chunk_files
  
  String dollar = "$"
  command <<<
    samtools merge \
      -f --threads 32 \
      ${in_sample_name}_merged.bam \
      ${sep=" " in_alignment_chunk_files} \
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

task runGenotyper {
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

