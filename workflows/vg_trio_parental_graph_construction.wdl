version 1.0

### vg_trio_parental_graph_construction.wdl ###
# Author: Charles Markello
# Description: Workflow for constructing  mother-father graph references.
#              Designed as the 3rd step in a pedigree-backed graph alignment pipeline.

import "https://github.com/vgteam/vg_wdl/blob/master/workflows/vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "https://github.com/vgteam/vg_wdl/blob/master/workflows/vg_construct_and_index.wdl" as vgConstructWorkflow

###########################
### WORKFLOW DEFINITION ###
###########################
workflow vgTrioPipeline {
    input {
        File COHORT_JOINT_VCF                               # Input joint-genotyped .vcf for mother-father-proband trio
        File COHORT_JOINT_VCF_INDEX                         # Input .vcf index of joint-genotyped mother-father-proband trio
        File MATERNAL_INPUT_BAM_FILE                        # Input maternal surjected .bam file
        File MATERNAL_INPUT_BAM_FILE_INDEX                  # Input maternal .bai index of surjected .bam file
        File PATERNAL_INPUT_BAM_FILE                        # Input paternal surjected .bam file
        File PATERNAL_INPUT_BAM_FILE_INDEX                  # Input paternal .bai index of surjected .bam file
        File PROBAND_INPUT_BAM_FILE                         # Input proband surjected .bam file
        File PROBAND_INPUT_BAM_FILE_INDEX                   # Input proband .bai index of surjected .bam file
        String SAMPLE_NAME_MATERNAL                         # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                         # Sample name for the father
        String SAMPLE_NAME_PROBAND                          # Sample name for the proband
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.19.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        Array[String]+ CONTIGS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]                                                          # List of contig strings used to construct new graph references with. Must match contigs in COHORT_JOINT_VCF
        File REF_FASTA_GZ                                   # Cannonical reference file .fasta.gz used for graph reference construction (only grch37/hg19 currently
                                                            #   supported)
        File PED_FILE                                       # .ped file describing the familial relationship of the trio samples
                                                            #   (https://whatshap.readthedocs.io/en/latest/guide.html#ped-file-format)
        File GEN_MAP_FILES                                  # .tar file containing recombination rates for each contig. Obtained from the 1000 Genomes genetic map 
                                                            #   (ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz)
        String GRAPH_NAME                                   # Name of the parental graph reference
        Boolean USE_HAPLOTYPES = true                       # Set to 'true' to construct the GBWT index which incorporates haplotype information into the graph.
        Boolean MAKE_SNARLS = false                         # Set to 'true' to construct the SNARLS index which incorporates indexes of "bubble" structures in the graph.
        Boolean USE_DECOYS = true                           # Set to 'true' to include decoy contigs from the FASTA reference into the graph reference.
        String DECOY_REGEX = ">GL\|>NC_007605\|>hs37d5"     # grep regular expression string that is used to extract decoy contig ids. USE_DECOYS must be set to 'true' to use this option.
    }
    
    
    #####################################
    ## Run parental graph construction ##
    #####################################
    if (USE_HAPLOTYPES) {
        call runSplitJointGenotypedVCF as splitJointGenotypedVCF {
            input:
                in_proband_sample_name=SAMPLE_NAME_PROBAND,
                in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
                in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
                joint_genotyped_vcf=COHORT_JOINT_VCF,
                joint_genotyped_vcf_index=COHORT_JOINT_VCF_INDEX,
                contigs=CONTIGS,
                filter_parents=false
        }
        scatter (contig_pair in zip(CONTIGS, splitJointGenotypedVCF.contig_vcfs)) {
            call runWhatsHapPhasing {
                input:
                    in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                    joint_genotyped_vcf=contig_pair.right,
                    in_maternal_bam=MATERNAL_INPUT_BAM_FILE,
                    in_maternal_bam_index=MATERNAL_INPUT_BAM_FILE_INDEX,
                    in_paternal_bam=PATERNAL_INPUT_BAM_FILE,
                    in_paternal_bam_index=PATERNAL_INPUT_BAM_FILE_INDEX,
                    in_proband_bam=PROBAND_INPUT_BAM_FILE,
                    in_proband_bam_index=PROBAND_INPUT_BAM_FILE_INDEX,
                    in_ped_file=PED_FILE,
                    in_genetic_map=GEN_MAP_FILES,
                    in_contig=contig_pair.left,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE
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
                in_vg_container=VG_CONTAINER
        }
    }
    File cohort_vcf_output = select_first([bgzipCohortPhasedVCF.output_merged_vcf, COHORT_JOINT_VCF])
    File cohort_vcf_index_output = select_first([bgzipCohortPhasedVCF.output_merged_vcf_index, COHORT_JOINT_VCF_INDEX])
    call runSplitJointGenotypedVCF as splitPhasedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=cohort_vcf_output,
            joint_genotyped_vcf_index=cohort_vcf_index_output,
            contigs=CONTIGS,
            filter_parents=true
    }
    call vgConstructWorkflow.vg_construct_and_index as constructGraphIndexWorkflow {
        input:
            graph_name=GRAPH_NAME,
            ref_fasta_gz=REF_FASTA_GZ,
            contigs=CONTIGS,
            contigs_vcf_gz=splitPhasedVCF.contig_vcfs,
            use_haplotypes=USE_HAPLOTYPES,
            make_snarls=MAKE_SNARLS,
            use_decoys=USE_DECOYS,
            decoy_regex=DECOY_REGEX,
            vg_docker=VG_CONTAINER
    }
    
    output {
        File parental_xg_index = constructGraphIndexWorkflow.xg
        File parental_gcsa_index = constructGraphIndexWorkflow.gcsa
        File parental_gcsa_lcp_index = constructGraphIndexWorkflow.gcsa_lcp
        File? parental_gbwt_index = constructGraphIndexWorkflow.gbwt
        File? parental_snarls_index = constructGraphIndexWorkflow.snarls
    }
}

########################
### TASK DEFINITIONS ###
########################


# Split a vcf into a list of vcfs each representing a contig
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
        set -exu -o pipefail
        
        if [ ~{filter_parents} == true ]; then
          SAMPLE_FILTER_STRING="-s ~{in_maternal_sample_name},~{in_paternal_sample_name}"
        else
          SAMPLE_FILTER_STRING=""
        fi
        ln -s ~{joint_genotyped_vcf} input_vcf_file.vcf.gz
        ln -s ~{joint_genotyped_vcf_index} input_vcf_file.vcf.gz.tbi
        while read -r contig; do
            if [[ ~{filter_parents} == true && ${contig} == "MT" ]]; then
                bcftools view -O z -r "${contig}" -s ~{in_maternal_sample_name} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            elif [[ ~{filter_parents} == true && ${contig} == "Y" ]]; then
                bcftools view -O z -r "${contig}" -s ~{in_paternal_sample_name} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            else
                bcftools view -O z -r "${contig}" ${SAMPLE_FILTER_STRING} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            fi
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

task runWhatsHapPhasing {
    input {
        String in_cohort_sample_name
        File joint_genotyped_vcf
        File? in_maternal_bam
        File? in_maternal_bam_index
        File? in_paternal_bam
        File? in_paternal_bam_index
        File? in_proband_bam
        File? in_proband_bam_index
        File in_ped_file
        File in_genetic_map
        String in_contig
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
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
        ln -s ~{in_maternal_bam} input_maternal.bam
        ln -s ~{in_maternal_bam_index} input_maternal.bam.bai
        ln -s ~{in_paternal_bam} input_paternal.bam
        ln -s ~{in_paternal_bam_index} input_paternal.bam.bai
        ln -s ~{in_proband_bam} input_proband.bam
        ln -s ~{in_proband_bam_index} input_proband.bam.bai
        whatshap phase \
            --reference ~{in_reference_file} \
            --indels \
            --ped ~{in_ped_file} \
            ${GENMAP_OPTION_STRING} \
            -o ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf \
            ~{joint_genotyped_vcf} input_proband.bam input_maternal.bam input_paternal.bam \
        && bgzip ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf
    >>>

    output {
        File phased_cohort_vcf = "~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf.gz"
    }
    runtime {
        memory: 50 + " GB"
        disks: "local-disk 100 SSD"
        docker: "quay.io/biocontainers/whatshap:0.18--py37h6bb024c_0"
    }
}

