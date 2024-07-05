version 1.0

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_reference_file, "G")) + 10
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
        cpu: 1
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task indexVcf {
    input {
        File in_vcf
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_vcf, "G")) + 10
    }

    command <<<
        set -eux -o pipefail
        bcftools index -t -o variants.vcf.gz.tbi ~{in_vcf}
    >>>
    output {
        File vcf_index_file = "variants.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        cpu: 1
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task fixVCFContigNaming {
    input {
        File in_vcf
        String in_prefix_to_strip = "GRCh38."
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_vcf, "G")) + 10
    }

    command <<<
    set -eux -o pipefail

    bcftools view -h ~{in_vcf} | grep contig | awk -v PREF="~{in_prefix_to_strip}" '{split($0, a, ","); gsub("##contig=<ID=", "", a[1]); chr=a[1]; gsub(PREF, "", a[1]); print chr"\t"a[1]}' > chrsrename.txt
    
    bcftools annotate --rename-chrs chrsrename.txt -o variants.vcf.gz -Oz ~{in_vcf}
    >>>
    output {
        File vcf = "variants.vcf.gz"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task removeHomRefs {
    input {
        File in_vcf
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_vcf, "G")) + 10
    }

    command <<<
    set -eux -o pipefail

    bcftools view -c 1 -o variants.vcf.gz -O z ~{in_vcf}
    >>>
    output {
        File vcf = "variants.vcf.gz"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task splitReads {
    input {
        File in_read_file
        String in_pair_id
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk = 5 * round(size(in_read_file, "G")) + 20
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

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int mem_gb = 8
        Int disk_size = 30 * round(size(in_clipped_vcf_chunk_files, "G")) + 50
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

        mkdir bcftools.tmp
        bcftools concat -n ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort -T bcftools.tmp -O z -o ${in_sample_name}.vcf.gz - && bcftools index -t -o ${in_sample_name}.vcf.gz.tbi ${in_sample_name}.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        time: 60
        memory: mem_gb + " GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
    }
}

task sortBAM {
    input {
        File in_bam_file
        File? in_ref_dict
        String in_prefix_to_strip = ""
        Int nb_cores = 16
        Int disk_size = 5 * round(size(in_bam_file, "G")) + 20
        String mem_gb = 16
    }

    String out_prefix = basename(in_bam_file, ".bam")
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

        if [ ~{in_prefix_to_strip} != "" ]
        then
            # patch the SQ fields from the dict into a new header
            samtools view -H ~{in_bam_file} | grep ^@HD > new_header.sam
            grep ^@SQ ~{in_ref_dict} | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
            samtools view -H ~{in_bam_file}  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam
            
            cat <(cat new_header.sam) <(samtools view ~{in_bam_file}) | \
                sed -e "s/~{in_prefix_to_strip}//g" | \
                samtools sort --threads ~{nb_cores} -O BAM > ~{out_prefix}.positionsorted.bam
        else
            samtools sort --threads ~{nb_cores} ~{in_bam_file} \
                     -O BAM > ~{out_prefix}.positionsorted.bam
            
        fi

        samtools index -b ~{out_prefix}.positionsorted.bam ~{out_prefix}.positionsorted.bam.bai
    >>>
    output {
        File sorted_bam = "~{out_prefix}.positionsorted.bam"
        File sorted_bam_index = "~{out_prefix}.positionsorted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: mem_gb + " GB"
        cpu: nb_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task leftShiftBAMFile {
    input {
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
        Int disk_size = round(3 * size(in_bam_file, 'G')) + 50
    }
    String out_prefix = basename(in_bam_file, ".bam")
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

task runAbraRealigner {
    input {
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk = round(3 * (size(in_bam_file, 'G') + size(in_reference_file, 'G'))) + 50
        Int memoryGb = 40
        Int threadCount = 16
    }
    String out_prefix = basename(in_bam_file, ".bam")
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

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java -Xmx~{memoryGb}G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{out_prefix}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads ~{threadCount}
    >>>
    output {
        File indel_realigned_bam = "~{out_prefix}.indel_realigned.bam"
        File indel_realigned_bam_index = "~{out_prefix}.indel_realigned.bai"
    }
    runtime {
        preemptible: 2
        time: 180
        memory: memoryGb + " GB"
        cpu: threadCount
        disks: "local-disk " + in_call_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task prepareRealignTargets {
    input {
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_expansion_bases
        Int disk_size = round(2 * size(in_bam_file, 'G')) + 20
        Int thread_count = 16
        Int mem_gb = 20
    }
    String out_prefix = basename(in_bam_file, ".bam")
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

        CONTIG_ID=`head -1 < <(samtools view input_bam_file.bam) | cut -f3`
        
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
          -nt ~{thread_count} \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{out_prefix}.intervals.bed

        if [ ~{in_expansion_bases} -gt 0 ]; then
            bedtools slop -i ~{out_prefix}.intervals.bed -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > ~{out_prefix}.intervals.widened.bed
            mv ~{out_prefix}.intervals.widened.bed ~{out_prefix}.intervals.bed
        fi
    >>>
    output {
        File output_target_bed_file = "~{out_prefix}.intervals.bed"
    }
    runtime {
        preemptible: 2
        time: 180
        memory: mem_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0"
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        String in_prefix_to_strip = ""
        Int thread_count = 16
        Int disk_size = round(3 * size(in_merged_bam_file, 'G')) + 20
        Int mem_gb = 20
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        if [ ~{in_prefix_to_strip} != "" ]
        then
            sed -e "s/~{in_prefix_to_strip}//g" ~{in_path_list_file} > paths.txt
        else
            cp ~{in_path_list_file} paths.txt
        fi
        
        while read -r contig; do
            samtools view \
              -@ ~{thread_count} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < paths.txt

        ## get unmapped reads
        mkdir unmapped
        samtools view \
                 -@ ~{thread_count} \
                 -h -O BAM \
                 -f 4 \
                 input_bam_file.bam \
                 -o unmapped/~{in_sample_name}.unmapped.bam \
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
        File bam_unmapped_file = glob("unmapped/~{in_sample_name}.*.bam")[0]
    }
    runtime {
        preemptible: 2
        memory: mem_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Int in_cores = 16
        Int disk_size = round(5 * size(in_alignment_bam_chunk_files, 'G')) + 20
        Int mem_gb = 8
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
          -f -p -c --threads ~{in_cores} \
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
        memory: mem_gb + " GB"
        cpu: in_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task convertCRAMtoFASTQ {
    input {
	    File? in_cram_file
        File? in_ref_file
        File? in_ref_index_file
        Boolean in_paired_reads = true
	    Int in_cores
        Int disk_size = round(5 * size(in_cram_file, 'G')) + 50
    }
    Int half_cores = in_cores / 2
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

    if [ ~{in_paired_reads} == true ]
    then
        samtools collate -@ ~{half_cores} --reference ~{in_ref_file} -Ouf ~{in_cram_file} | samtools fastq -@ ~{half_cores} -1 reads.R1.fastq.gz -2 reads.R2.fastq.gz -0 reads.o.fq.gz -s reads.s.fq.gz -c 1 -N -
    else
        samtools fastq -@ ~{in_cores} -o reads.R1.fastq.gz -c 1 --reference ~{in_ref_file} ~{in_cram_file}
    fi
    >>>
    output {
        File output_fastq_1_file = "reads.R1.fastq.gz"
        File? output_fastq_2_file = "reads.R2.fastq.gz"
    }
    runtime {
        preemptible: 2
        cpu: in_cores
        memory: "50 GB"
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }    
}


task kmerCountingKMC {
    input {
        File input_read_file_1
        File? input_read_file_2
        String output_file_name

        Int kmer_length
        Int max_ram = 64
	    Int nb_cores = 16
        Int disk_size = round(size(input_read_file_1, "G") + size(input_read_file_2, "G")) + 10
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

    echo ~{input_read_file_1} > scratch_file.txt
    ~{if defined(input_read_file_2) then "echo ~{input_read_file_2} >> scratch_file.txt" else ""}

    kmc -k~{kmer_length} -m~{max_ram} -okff -t~{nb_cores} @scratch_file.txt ~{output_file_name} .

    rm scratch_file.txt
    >>>
    output {
        File kff_file = output_file_name + ".kff"
    }
    runtime {
        preemptible: 2
        cpu: nb_cores
        memory: max_ram + " GB"
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/kmc:3.2.1--hf1761c0_2"
    }
}
