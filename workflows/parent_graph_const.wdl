version 1.0

import "./vg_construct_and_index.wdl" as vgConstructWorkflow
import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow

workflow vgParentGraphConstruction {
    input {
        String SAMPLE_NAME_MATERNAL
        String SAMPLE_NAME_PATERNAL
        String SAMPLE_NAME_PROBAND
        File JOINT_GENOTYPED_VCF
        File JOINT_GENOTYPED_VCF_INDEX
        File MATERNAL_BAM
        File MATERNAL_BAM_INDEX
        File PATERNAL_BAM
        File PATERNAL_BAM_INDEX
        File PROBAND_BAM
        File PROBAND_BAM_INDEX
        File PED_FILE
        File GEN_MAP_FILES              # zipped list of genetic recombination rate files
        String GRAPH_NAME
        File REF_FASTA_GZ
        Array[String]+ contigs = ["21"]
        Boolean use_haplotypes = true
        Boolean make_snarls = true 
        Boolean use_decoys = true 
        String decoy_regex = ">GL"
        String vg_docker = "quay.io/vgteam/vg:dev-v1.14.0-242-g39ac300b-t293-run"
    }
    call runSplitJointGenotypedVCF as splitJointGenotypedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=JOINT_GENOTYPED_VCF,
            joint_genotyped_vcf_index=JOINT_GENOTYPED_VCF_INDEX,
            contigs=contigs,
            filter_parents=false
    }
    scatter (contig_pair in zip(contigs, splitJointGenotypedVCF.contig_vcfs)) {
        call runWhatsHapPhasing {
            input:
                in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                joint_genotyped_vcf=contig_pair.right,
                in_maternal_bam=MATERNAL_BAM,
                in_maternal_bam_index=MATERNAL_BAM_INDEX,
                in_paternal_bam=PATERNAL_BAM,
                in_paternal_bam_index=PATERNAL_BAM_INDEX,
                in_proband_bam=PROBAND_BAM,
                in_proband_bam_index=PROBAND_BAM_INDEX,
                in_ped_file=PED_FILE,
                in_genetic_map=GEN_MAP_FILES,
                in_contig=contig_pair.left
        }
    }
    call vgMultiMapCallWorkflow.concatClippedVCFChunks as concatCohortPhasedVCFs {
        input:
            in_sample_name="${SAMPLE_NAME_PROBAND}_cohort",
            in_clipped_vcf_chunk_files=runWhatsHapPhasing.phased_cohort_vcf
    }
    call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipCohortPhasedVCF {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_merged_vcf_file=concatCohortPhasedVCFs.output_merged_vcf,
            in_vg_container=vg_docker
    }
    call runSplitJointGenotypedVCF as splitPhasedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=bgzipCohortPhasedVCF.output_merged_vcf,
            joint_genotyped_vcf_index=bgzipCohortPhasedVCF.output_merged_vcf_index,
            contigs=contigs,
            filter_parents=true
    }
    call vgConstructWorkflow.vg_construct_and_index as constructGraphIndexWorkflow {
        input:
            graph_name=GRAPH_NAME,
            ref_fasta_gz=REF_FASTA_GZ,
            contigs=contigs,
            contigs_vcf_gz=splitPhasedVCF.contig_vcfs,
            use_haplotypes=use_haplotypes,
            make_snarls=make_snarls,
            use_decoys=use_decoys,
            decoy_regex=decoy_regex,
            vg_docker=vg_docker
    }
    
    output {
        File vg = constructGraphIndexWorkflow.vg
        File xg = constructGraphIndexWorkflow.xg
        File gcsa = constructGraphIndexWorkflow.gcsa
        File gcsa_lcp = constructGraphIndexWorkflow.gcsa_lcp
        File? gbwt = constructGraphIndexWorkflow.gbwt
        File? snarls = constructGraphIndexWorkflow.snarls
    }
}

task runWhatsHapPhasing {
    input {
        String in_cohort_sample_name
        File joint_genotyped_vcf
        File in_maternal_bam
        File in_maternal_bam_index
        File in_paternal_bam
        File in_paternal_bam_index
        File in_proband_bam
        File in_proband_bam_index
        File in_ped_file
        File in_genetic_map
        String in_contig
    }
    
    command <<<
        set -exu -o pipefail
        
        tar -xvf ~{in_genetic_map}
        if [[ ~{in_contig} == "Y" || ~{in_contig} == "MT" ]]; then
            GENMAP_OPTION_STRING=""
        elif [ ~{in_contig} == "X" ]; then
            GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chrX_nonPAR_combined_b37.txt --chromosome X"
        else
            GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chr~{in_contig}_combined_b37.txt --chromosome ~{in_contig}"
        fi
        whatshap phase \
            --indels \
            --ped ~{in_ped_file} \
            ${GENMAP_OPTION_STRING} \
            -o ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf \
            ~{joint_genotyped_vcf} ~{in_proband_bam} ~{in_maternal_bam} ~{in_paternal_bam}
    >>>
    
    output {
        File phased_cohort_vcf = "~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/whatshap:0.18--py37h6bb024c_0"
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
        Boolean filter_parents
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
        
        if [ ~{filter_parents} == true ]; then
          SAMPLE_FILTER_STRING="-s ~{in_maternal_sample_name},~{in_paternal_sample_name}"
        else
          SAMPLE_FILTER_STRING=""
        fi
        
        while read -r contig; do
            bcftools view -O z -r "${contig}" ${SAMPLE_FILTER_STRING} ~{joint_genotyped_vcf} > "${contig}.vcf.gz"
            echo "${contig}.vcf.gz" >> contig_vcf_list.txt
        done < "~{write_lines(contigs)}"
    >>>
    output {
        Array[File]+ contig_vcfs = read_lines("contig_vcf_list.txt")
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    }
}



