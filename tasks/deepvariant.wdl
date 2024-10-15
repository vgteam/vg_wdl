version 1.0

# Run DeepVariant example generation on a single-contig BAM named <sample>.<contig>(.indel_realigned)?(.left_shifted)?.bam
task runDeepVariantMakeExamples {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        String in_model_type = "WGS"
        Array[File] in_model_files = []
        Array[File] in_model_variables_files = []
        Int? in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        String in_other_makeexamples_arg = ""
        Boolean in_dv_is_1_7_or_newer = false # Needs to be True for DV 1.7+
        Int in_call_cores
        Int in_call_mem
        String in_dv_container = "google/deepvariant:1.5.0"
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

        MODEL_TYPE=~{in_model_type}
        MODEL_TYPE_ARGS=()
        # Determine extra make_example arguments for the given model type just like in DeepVariant's main wrapper script.
        # See <https://github.com/google/deepvariant/blob/ab068c4588a02e2167051bd9e74c0c9579462b51/scripts/run_deepvariant.py#L243-L276>
        # Except instead of building a channel list we load it form the model.
        case ${MODEL_TYPE} in
            WGS)
                if [[ ~{in_dv_is_1_7_or_newer} == false ]] ; then
                    MODEL_TYPE_ARGS+=(--channels insert_size)
                fi
                if [ ~{defined(in_min_mapq)} == true ]; then
                    # Add our min MAPQ override
                    MODEL_TYPE_ARGS+=(--min_mapping_quality ~{in_min_mapq})
                fi
                ;;

            WES)
                if [[ ~{in_dv_is_1_7_or_newer} == false ]] ; then
                    MODEL_TYPE_ARGS+=(--channels insert_size)
                fi
                if [ ~{defined(in_min_mapq)} == true ]; then
                    # Add our min MAPQ override
                    MODEL_TYPE_ARGS+=(--min_mapping_quality ~{in_min_mapq})
                fi
                ;;

            PACBIO)
                if [[ ~{in_dv_is_1_7_or_newer} == false ]] ; then
                    MODEL_TYPE_ARGS+=(--add_hp_channel)
                fi
                MODEL_TYPE_ARGS+=(--alt_aligned_pileup 'diff_channels')
                MODEL_TYPE_ARGS+=(--max_reads_per_partition 600)
                MODEL_TYPE_ARGS+=(--min_mapping_quality ~{select_first([in_min_mapq, 1])})
                MODEL_TYPE_ARGS+=(--parse_sam_aux_fields)
                MODEL_TYPE_ARGS+=(--partition_size 25000)
                MODEL_TYPE_ARGS+=(--phase_reads)
                MODEL_TYPE_ARGS+=(--pileup_image_width 199)
                MODEL_TYPE_ARGS+=(--norealign_reads)
                MODEL_TYPE_ARGS+=(--sort_by_haplotypes)
                MODEL_TYPE_ARGS+=(--track_ref_reads)
                MODEL_TYPE_ARGS+=(--vsc_min_fraction_indels 0.12)
                if [[ ~{in_dv_is_1_7_or_newer} == true ]] ; then
                    MODEL_TYPE_ARGS+=(--trim_reads_for_pileup)
                fi
                ;;

            ONT_R104)
                if [[ ~{in_dv_is_1_7_or_newer} == false ]] ; then
                    MODEL_TYPE_ARGS+=(--add_hp_channel)
                fi
                MODEL_TYPE_ARGS+=(--alt_aligned_pileup 'diff_channels')
                MODEL_TYPE_ARGS+=(--max_reads_per_partition 600)
                MODEL_TYPE_ARGS+=(--min_mapping_quality ~{select_first([in_min_mapq, 5])})
                MODEL_TYPE_ARGS+=(--parse_sam_aux_fields)
                MODEL_TYPE_ARGS+=(--partition_size 25000)
                MODEL_TYPE_ARGS+=(--phase_reads)
                MODEL_TYPE_ARGS+=(--pileup_image_width 199)
                MODEL_TYPE_ARGS+=(--norealign_reads)
                MODEL_TYPE_ARGS+=(--sort_by_haplotypes)
                MODEL_TYPE_ARGS+=(--track_ref_reads)
                MODEL_TYPE_ARGS+=(--vsc_min_fraction_snps 0.08)
                MODEL_TYPE_ARGS+=(--vsc_min_fraction_indels 0.12)
                if [[ ~{in_dv_is_1_7_or_newer} == true ]] ; then
                    MODEL_TYPE_ARGS+=(--trim_reads_for_pileup)
                fi
                ;;
            
            HYBRID_PACBIO_ILLUMINA)
                if [ ~{defined(in_min_mapq)} == true ]; then
                    # Add our min MAPQ override
                    MODEL_TYPE_ARGS+=(--min_mapping_quality ~{in_min_mapq})
                fi
                if [[ ~{in_dv_is_1_7_or_newer} == true ]] ; then
                    MODEL_TYPE_ARGS+=(--trim_reads_for_pileup)
                fi
                ;;
        esac

        # Set up the model so DV can read the channels out of it
        if [[ ~{length(in_model_files)} -gt 0 ]] ; then
            # Need to use a custom model
            mkdir model_dir
            ln -s ~{sep=" " in_model_files} model_dir/
            if [[ ~{length(in_model_variables_files)} -gt 0 ]] ; then
                # Some models (like the DV release default models) also have a "variables" subdirectory. Handle it specially.
                # TODO: Is it possible to iterate over a WDL Map in Bash so we can just send the whole structure?
                mkdir model_dir/variables
                ln -s ~{sep=" " in_model_variables_files} model_dir/variables/
            fi
        else
            # Use default models for type
            ln -s /opt/models/${MODEL_TYPE,,} model_dir
        fi
        ls -lah model_dir >&2
        CHECKPOINT_INDEX_FILES=(model_dir/*.ckpt.index)
        if [[ -e "${CHECKPOINT_INDEX_FILES[0]}" ]] ; then
            # This is a checkpoint-format model and we need to name it by passing this path without the .index
            CHECKPOINT_NAME="${CHECKPOINT_INDEX_FILES[0]%.index}"
        else
            # This is a savedmodel-format model and is named just by the directory
            CHECKPOINT_NAME="model_dir"
        fi
        CHECKPOINT_ARGS=()
        if [[ ~{in_dv_is_1_7_or_newer} == true ]] ; then
            # We only actually show the model to DV if we want to use the channels from it. Older DV can't take it here.
            CHECKPOINT_ARGS+=(--checkpoint "${CHECKPOINT_NAME}")
        fi

        seq 0 $((~{in_call_cores}-1)) | \
        parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
        --mode calling \
        "${CHECKPOINT_ARGS[@]}" \
        --ref reference.fa \
        --reads input_bam_file.bam \
        --examples ./make_examples.tfrecord@~{in_call_cores}.gz \
        --sample_name ~{in_sample_name} \
        --gvcf ./gvcf.tfrecord@~{in_call_cores}.gz \
        ${KEEP_LEGACY_AC_ARG} ${NORM_READS_ARG} "${MODEL_TYPE_ARGS[@]}" ~{in_other_makeexamples_arg} \
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

# Run DeepVariant calling AND postprocessing on a file of examples, on GPUs.
task runDeepVariantCallVariants {
    input {
        String in_sample_name
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
        String in_model_type = "WGS"
        Array[File] in_model_files = []
        Array[File] in_model_variables_files = []
        # Not available on DV1.5
        Array[String] in_haploid_contigs = []
        # Not available on DV1.5
        File? in_par_regions_bed_file
        Int in_call_cores
        Int in_call_mem
        Boolean in_use_gpus = true
        String in_dv_gpu_container = "google/deepvariant:1.5.0-gpu"
    }
    Int disk_size = 5 * round(size(in_examples_file, 'G') + size(in_nonvariant_site_tf_file, 'G') + size(in_reference_file, 'G')) + 50
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

        MODEL_TYPE=~{in_model_type}

        # Set up the model
        if [[ ~{length(in_model_files)} -gt 0 ]] ; then
            # Need to use a custom model
            mkdir model_dir
            ln -s ~{sep=" " in_model_files} model_dir/
            if [[ ~{length(in_model_variables_files)} -gt 0 ]] ; then
                # Some models (like the DV release default models) also have a "variables" subdirectory. Handle it specially.
                # TODO: Is it possible to iterate over a WDL Map in Bash so we can just send the whole structure?
                mkdir model_dir/variables
                ln -s ~{sep=" " in_model_variables_files} model_dir/variables/
            fi
        else
            # Use default models for type
            ln -s /opt/models/${MODEL_TYPE,,} model_dir
        fi
        ls -lah model_dir >&2
        CHECKPOINT_INDEX_FILES=(model_dir/*.ckpt.index)
        if [[ -e "${CHECKPOINT_INDEX_FILES[0]}" ]] ; then
            # This is a checkpoint-format model and we need to name it by passing this path without the .index
            CHECKPOINT_NAME="${CHECKPOINT_INDEX_FILES[0]%.index}"
        else
            # This is a savedmodel-format model and is named just by the directory
            CHECKPOINT_NAME="model_dir"
        fi
        
        /opt/deepvariant/bin/call_variants \
        --outfile call_variants_output.tfrecord.gz \
        --examples "make_examples.tfrecord@~{in_call_cores}.gz" \
        --checkpoint "${CHECKPOINT_NAME}" && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref reference.fa \
        --infile call_variants_output.tfrecord.gz \
        --nonvariant_site_tfrecord_path "gvcf.tfrecord@~{in_call_cores}.gz" \
        --outfile "~{in_sample_name}_deepvariant.vcf.gz" \
        ~{if length(in_haploid_contigs) > 0 then "--haploid_contigs" else ""} ~{sep=',' in_haploid_contigs} \
        ~{"--par_regions_bed " + in_par_regions_bed_file} \
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
        gpu: in_use_gpus
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + disk_size + " SSD"
        docker: in_dv_gpu_container
    }
}


