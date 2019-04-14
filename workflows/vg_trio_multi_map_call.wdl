version 1.0
#version draft-2
# Apparently cromwell 36.1 doesn't parse the version line in wdl files if the wdl version is a draft version.
# Cromwell 36.1 wdl parser defaults to draft-2 if a `version` line is not given.

# TODO:
# - Use subworkflow functionality in new vg_trio_sub.wdl workflow
#   - DOCS: https://cromwell.readthedocs.io/en/stable/SubWorkflows/
# - define DRAGEN joint genotypeing task
# - possibly define genotype phasing task
#     - Start with just the regular mpmap gcsa and xg construction

import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow

workflow vgTrioPipeline {
    input {
        File MATERNAL_INPUT_READ_FILE_1
        File MATERNAL_INPUT_READ_FILE_2
        File PATERNAL_INPUT_READ_FILE_1
        File PATERNAL_INPUT_READ_FILE_2
        File PROBAND_INPUT_READ_FILE_1
        File PROBAND_INPUT_READ_FILE_2
        String SAMPLE_NAME_MATERNAL
        String SAMPLE_NAME_PATERNAL
        String SAMPLE_NAME_PROBAND
        String VG_CONTAINER
        Int READS_PER_CHUNK
        Int CHUNK_BASES
        Int OVERLAP
        File PATH_LIST_FILE
        File PATH_LENGTH_FILE
        File XG_FILE
        File GCSA_FILE
        File GCSA_LCP_FILE
        File GBWT_FILE
        File SNARLS_FILE
        File REF_FILE
        File REF_INDEX_FILE
        File REF_DICT_FILE
        File SNPEFF_DATABASE
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
        String DRAGEN_REF_INDEX_NAME
        String UDPBINFO_PATH
        String HELIX_USERNAME
        Boolean RUN_VGMPMAP_ALGORITHM
        Boolean RUN_DRAGEN_CALLER
        Array[String]+ CONTIGS = ["20", "21"]
        File REF_FASTA_GZ
        Boolean USE_HAPLOTYPES = false
    }
    
    call vgMultiMapCallWorkflow.vgMultiMapCall as maternalMapCallWorkflow {
        input:
            INPUT_READ_FILE_1=MATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=MATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            PATH_LENGTH_FILE=PATH_LENGTH_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            RUN_VGMPMAP_ALGORITHM=RUN_VGMPMAP_ALGORITHM,
            RUN_DRAGEN_CALLER=false,
            RUN_LINEAR_CALLER=true,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false
    }
    call vgMultiMapCallWorkflow.vgMultiMapCall as paternalMapCallWorkflow {
        input:
            INPUT_READ_FILE_1=PATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            PATH_LENGTH_FILE=PATH_LENGTH_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            RUN_VGMPMAP_ALGORITHM=RUN_VGMPMAP_ALGORITHM,
            RUN_DRAGEN_CALLER=false,
            RUN_LINEAR_CALLER=true,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false
    }
    call vgMultiMapCallWorkflow.vgMultiMapCall as probandMapCallWorkflow {
        input:
            INPUT_READ_FILE_1=PROBAND_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PROBAND_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            PATH_LENGTH_FILE=PATH_LENGTH_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            RUN_VGMPMAP_ALGORITHM=RUN_VGMPMAP_ALGORITHM,
            RUN_DRAGEN_CALLER=false,
            RUN_LINEAR_CALLER=true,
            RUN_GVCF_OUTPUT=true,
            RUN_SNPEFF_ANNOTATION=false
    }
    call runGATKCombineGenotypeGVCFs {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_gvcf_file_maternal=maternalMapCallWorkflow.output_vcf,
            in_gvcf_file_paternal=paternalMapCallWorkflow.output_vcf,
            in_gvcf_file_proband=probandMapCallWorkflow.output_vcf,
            in_reference_file=REF_FILE,
            in_reference_index_file=REF_INDEX_FILE,
            in_reference_dict_file=REF_DICT_FILE
    
    }   
    call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipGATKGVCF {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_merged_vcf_file=runGATKCombineGenotypeGVCFs.joint_genotyped_vcf,
            in_vg_container=VG_CONTAINER
    }
    
    call runSplitJointGenotypedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=bgzipGATKGVCF.output_merged_vcf,
            joint_genotyped_vcf_index=bgzipGATKGVCF.output_merged_vcf_index,
            contigs=CONTIGS
    }
    call vgConstructWorkflow.vg_construct_and_index as constructGraphIndexWorkflow {
        input:
            graph_name="parental_graph",
            ref_fasta_gz=REF_FASTA_GZ,
            contigs=CONTIGS,
            contigs_vcf_gz=runSplitJointGenotypedVCF.contig_vcfs,
            use_haplotypes=USE_HAPLOTYPES,
            vg_docker=VG_CONTAINER
    }

    call vgMultiMapCallWorkflow.vgMultiMapCall as probandMapCallWorkflow2ndIteration {
        input:
            INPUT_READ_FILE_1=PROBAND_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PROBAND_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            PATH_LENGTH_FILE=PATH_LENGTH_FILE,
            XG_FILE=constructGraphIndexWorkflow.xg,
            GCSA_FILE=constructGraphIndexWorkflow.gcsa,
            GCSA_LCP_FILE=constructGraphIndexWorkflow.gcsa_lcp,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            RUN_VGMPMAP_ALGORITHM=RUN_VGMPMAP_ALGORITHM,
            RUN_DRAGEN_CALLER=RUN_DRAGEN_CALLER,
            RUN_LINEAR_CALLER=true,
            RUN_GVCF_OUTPUT=false,
            RUN_SNPEFF_ANNOTATION=true
    }
    
    output {
        File output_cohort_vcf = bgzipGATKGVCF.output_merged_vcf
        File? output_maternal_bam = maternalMapCallWorkflow.output_bam
        File? output_paternal_bam = paternalMapCallWorkflow.output_bam
        File? output_proband_bam = probandMapCallWorkflow.output_bam
        File output_final_proband_vcf = probandMapCallWorkflow2ndIteration.output_vcf
        File? output_final_proband_bam = probandMapCallWorkflow2ndIteration.output_bam
    }
}

task runGATKCombineGenotypeGVCFs {
    input {
        String in_sample_name
        File in_gvcf_file_maternal
        File in_gvcf_file_paternal
        File in_gvcf_file_proband
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

        gatk CombineGVCFs \
          --reference ${in_reference_file} \
          -V ${in_gvcf_file_maternal} -V ${in_gvcf_file_paternal} -V ${in_gvcf_file_proband} \
          --output ${in_sample_name}_cohort.combined.gvcf \
        && gatk GenotypeGVCFs \
          --reference ${in_reference_file} \
          --variant ${in_sample_name}_cohort.combined.gvcf \
          --output ${in_sample_name}_cohort.jointgenotyped.vcf
    }
    output {
        File joint_genotyped_vcf = "${in_sample_name}_cohort.jointgenotyped.vcf"
    }
    runtime {
        memory: 100
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runSplitJointGenotypedVCF {
    input {
        String in_proband_sample_name
        String in_maternal_sample_name
        String in_paternal_sample_name
        File joint_genotyped_vcf
        File joint_genotyped_vcf_index
        Array[String]+ contigs
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

        while read -r contig; do
            bcftools view -O z -r "${contig}" -s ~{in_maternal_sample_name},~{in_paternal_sample_name} ~{joint_genotyped_vcf} > "${contig}.vcf.gz"
            echo "${contig}.vcf.gz" >> contig_vcf_list.txt
        done < "~{write_lines(contigs)}"
    >>>
    output {
        Array[File]+ contig_vcfs = read_lines("contig_vcf_list.txt")
    }
    runtime {
        memory: 50
        disks: 100
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}

