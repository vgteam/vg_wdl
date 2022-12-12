version 1.0

task runDeepVariantMakeExamples {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        Int in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        String in_other_makeexamples_arg = ""
        Int in_call_cores
        Int in_call_mem
        String in_dv_container = "google/deepvariant:1.3.0"
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
        ${KEEP_LEGACY_AC_ARG} ${NORM_READS_ARG} ~{in_other_makeexamples_arg} \
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
        String in_dv_gpu_container = "google/deepvariant:1.3.0-gpu"
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

