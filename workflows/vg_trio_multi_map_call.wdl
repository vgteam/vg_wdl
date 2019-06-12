version 1.0

### vg_trio_multi_map_call.wdl ###
## Author: Charles Markello
## Description: Trio-backed VG mapping and variant calling workflow for mother-father-child trio datasets.
## Reference: https://github.com/vgteam/vg/wiki

import "./vg_multi_map_call.wdl" as vgMultiMapCallWorkflow
import "./vg_multi_map.wdl" as vgMultiMapWorkflow
import "./vg_multi_call.wdl" as vgMultiCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow

workflow vgTrioPipeline {
    input {
        File MATERNAL_INPUT_READ_FILE_1                     # Input maternal 1st read pair fastq.gz
        File MATERNAL_INPUT_READ_FILE_2                     # Input maternal 2nd read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_1                     # Input paternal 1st read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_2                     # Input paternal 2nd read pair fastq.gz
        File PROBAND_INPUT_READ_FILE_1                      # Input proband 1st read pair fastq.gz
        File PROBAND_INPUT_READ_FILE_2                      # Input proband 2nd read pair fastq.gz
        String SAMPLE_NAME_MATERNAL                         # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                         # Sample name for the father
        String SAMPLE_NAME_PROBAND                          # Sample name for the proband
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.16.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 20000000                      # Number of reads contained in each mapping chunk (20000000 for wgs)
        Int CHUNK_BASES = 50000000                          # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                                  # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File GCSA_FILE                                      # Path to .gcsa index file
        File GCSA_LCP_FILE                                  # Path to .gcsa.lcp index file
        File? GBWT_FILE                                     # (OPTIONAL) Path to .gbwt index file
        File? SNARLS_FILE                                   # (OPTIONAL) Path to .snarls index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        File SNPEFF_DATABASE                                # Path to snpeff database .zip file for snpEff annotation functionality.
        Int SPLIT_READ_CORES = 32
        Int SPLIT_READ_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAM_CORES = 56
        Int MERGE_GAM_DISK = 400
        Int MERGE_GAM_MEM = 100
        Int MERGE_GAM_TIME = 1200
        Int CHUNK_GAM_CORES = 32
        Int CHUNK_GAM_DISK = 400
        Int CHUNK_GAM_MEM = 100
        Int VGCALL_CORES = 8
        Int VGCALL_DISK = 40
        Int VGCALL_MEM = 64
        String DRAGEN_REF_INDEX_NAME                        # Dragen module based reference index directory (e.g. "hs37d5_v7")
        String UDPBINFO_PATH                                # Udp data directory to use for Dragen module (e.g. "Udpbinfo", nih biowulf system only)
        String HELIX_USERNAME                               # The nih helix username which holds a user directory in UDPBINFO_PATH
        Array[String]+ CONTIGS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]
        File REF_FASTA_GZ
        File PED_FILE
        File GEN_MAP_FILES
        String GRAPH_NAME
        Boolean VGMPMAP_MODE                            # Set to 'false' to use "VG MAP" or set to 'true' to use "VG MPMAP" algorithm
        Boolean DRAGEN_MODE                             # Set to 'true' to use the Dragen modules variant caller. Set to 'false' to use GATK HaplotypeCallers genotyper.
        Boolean USE_HAPLOTYPES = true                   # Set to 'true' to construct the GBWT index which incorporates haplotype information into the graph.
        Boolean MAKE_SNARLS = true                      # Set to 'true' to construct the SNARLS index which incorporates indexes of "bubble" structures in the graph.
        Boolean USE_DECOYS = true                       # Set to 'true' to include decoy contigs from the FASTA reference into the graph reference.
        Boolean SNPEFF_ANNOTATION = true                # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        String DECOY_REGEX = ">GL"                      # grep regular expression string that is used to extract decoy contig ids. USE_DECOYS must be set to 'true'
    }
    
    ###################################
    ## Run mapping workflows on Trio ##
    ###################################
    call vgMultiMapWorkflow.vgMultiMapCall as maternalMapWorkflow {
        input:
            INPUT_READ_FILE_1=MATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=MATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            VGMPMAP_MODE=VGMPMAP_MODE,
            SURJECT_MODE=true
    }
    call vgMultiMapWorkflow.vgMultiMapCall as paternalMapWorkflow {
        input:
            INPUT_READ_FILE_1=PATERNAL_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PATERNAL_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            VGMPMAP_MODE=VGMPMAP_MODE,
            SURJECT_MODE=true
    }
    call vgMultiMapWorkflow.vgMultiMapCall as probandMapWorkflow {
        input:
            INPUT_READ_FILE_1=PROBAND_INPUT_READ_FILE_1,
            INPUT_READ_FILE_2=PROBAND_INPUT_READ_FILE_2,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            VG_CONTAINER=VG_CONTAINER,
            READS_PER_CHUNK=READS_PER_CHUNK,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            GCSA_FILE=GCSA_FILE,
            GCSA_LCP_FILE=GCSA_LCP_FILE,
            GBWT_FILE=GBWT_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MERGE_GAM_CORES=MERGE_GAM_CORES,
            MERGE_GAM_DISK=MERGE_GAM_DISK,
            MERGE_GAM_MEM=MERGE_GAM_MEM,
            MERGE_GAM_TIME=MERGE_GAM_TIME,
            VGMPMAP_MODE=VGMPMAP_MODE,
            SURJECT_MODE=true
    }
    
    ###########################################
    ## Run variant calling workflows on Trio ##
    ###########################################
    call vgMultiCallWorkflow.vgMultiMapCall as maternalCallWorkflow {
        input:
            INPUT_BAM_FILE=maternalMapWorkflow.output_bam,
            INPUT_BAM_FILE_INDEX=maternalMapWorkflow.output_bam_index,
            SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            SURJECT_MODE=true,
            DRAGEN_MODE=DRAGEN_MODE,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false
    }
    call vgMultiCallWorkflow.vgMultiMapCall as paternalCallWorkflow {
        input:
            INPUT_BAM_FILE=paternalMapWorkflow.output_bam,
            INPUT_BAM_FILE_INDEX=paternalMapWorkflow.output_bam_index,
            SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
            VG_CONTAINER=VG_CONTAINER,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            SURJECT_MODE=true,
            DRAGEN_MODE=DRAGEN_MODE,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT= if DRAGEN_MODE then maternalCallWorkflow.output_vcf else "null"
    }
    call vgMultiCallWorkflow.vgMultiMapCall as probandCallWorkflow {
        input:
            INPUT_BAM_FILE=probandMapWorkflow.output_bam,
            INPUT_BAM_FILE_INDEX=probandMapWorkflow.output_bam_index,
            SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            VG_CONTAINER=VG_CONTAINER,
            CHUNK_BASES=CHUNK_BASES,
            OVERLAP=OVERLAP,
            PATH_LIST_FILE=PATH_LIST_FILE,
            XG_FILE=XG_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SNPEFF_DATABASE=SNPEFF_DATABASE,
            CHUNK_GAM_CORES=CHUNK_GAM_CORES,
            CHUNK_GAM_DISK=CHUNK_GAM_DISK,
            CHUNK_GAM_MEM=CHUNK_GAM_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM,
            DRAGEN_REF_INDEX_NAME=DRAGEN_REF_INDEX_NAME,
            UDPBINFO_PATH=UDPBINFO_PATH,
            HELIX_USERNAME=HELIX_USERNAME,
            SURJECT_MODE=true,
            DRAGEN_MODE=DRAGEN_MODE,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false,
            PREVIOUS_WORKFLOW_OUTPUT= if DRAGEN_MODE then paternalCallWorkflow.output_vcf else "null"
    }
    
    ###############################
    ## Run trio joint genotyping ##
    ###############################
    if (!RUN_DRAGEN_CALLER) {
        call runGATKCombineGenotypeGVCFs as gatkJointGenotyper1st {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandCallWorkflow.output_vcf,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipGATKGVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_merged_vcf_file=gatkJointGenotyper1st.joint_genotyped_vcf,
                in_vg_container=VG_CONTAINER
        }
    }
    if (RUN_DRAGEN_CALLER) {
        call runDragenJointGenotyper as dragenJointGenotyper1st {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandCallWorkflow.output_vcf,
                in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                in_udp_data_dir=UDPBINFO_PATH,
                in_helix_username=HELIX_USERNAME
        }
    }
    
    #####################################
    ## Run parental graph construction ##
    #####################################
    File output_joint_genotyped_vcf = select_first([bgzipGATKGVCF.output_merged_vcf, dragenJointGenotyper1st.dragen_joint_genotyped_vcf])
    File output_joint_genotyped_vcf_index = select_first([bgzipGATKGVCF.output_merged_vcf_index, dragenJointGenotyper1st.dragen_joint_genotyped_vcf_index])
    call runSplitJointGenotypedVCF as splitJointGenotypedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=output_joint_genotyped_vcf,
            joint_genotyped_vcf_index=output_joint_genotyped_vcf_index,
            contigs=CONTIGS,
            filter_parents=false
    }
    scatter (contig_pair in zip(CONTIGS, splitJointGenotypedVCF.contig_vcfs)) {
        call runWhatsHapPhasing {
            input:
                in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                joint_genotyped_vcf=contig_pair.right,
                in_maternal_bam=maternalMapWorkflow.output_bam,
                in_maternal_bam_index=maternalMapWorkflow.output_bam_index,
                in_paternal_bam=paternalMapWorkflow.output_bam,
                in_paternal_bam_index=paternalMapWorkflow.output_bam_index,
                in_proband_bam=probandMapWorkflow.output_bam,
                in_proband_bam_index=probandMapWorkflow.output_bam_index,
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
    call runSplitJointGenotypedVCF as splitPhasedVCF {
        input:
            in_proband_sample_name=SAMPLE_NAME_PROBAND,
            in_maternal_sample_name=SAMPLE_NAME_MATERNAL,
            in_paternal_sample_name=SAMPLE_NAME_PATERNAL,
            joint_genotyped_vcf=bgzipCohortPhasedVCF.output_merged_vcf,
            joint_genotyped_vcf_index=bgzipCohortPhasedVCF.output_merged_vcf_index,
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
    
    ##################################################################################
    ## Run mapping and calling workflow of proband reads against the parental graph ##
    ##################################################################################
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
            XG_FILE=constructGraphIndexWorkflow.xg,
            GCSA_FILE=constructGraphIndexWorkflow.gcsa,
            GCSA_LCP_FILE=constructGraphIndexWorkflow.gcsa_lcp,
            GBWT_FILE=constructGraphIndexWorkflow.gbwt,
            SNARLS_FILE=constructGraphIndexWorkflow.snarls,
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
            VGMPMAP_MODE=VGMPMAP_MODE,
            DRAGEN_MODE=DRAGEN_MODE,
            SURJECT_MODE=true,
            GVCF_MODE=true,
            SNPEFF_ANNOTATION=false
    }
    
    #######################################################
    ## Run 2nd trio joint genotyping on new proband GVCF ##
    #######################################################
    if (!DRAGEN_MODE) {
        call runGATKCombineGenotypeGVCFs as gatkJointGenotyper2nd {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandMapCallWorkflow2ndIteration.output_vcf,
                in_reference_file=REF_FILE,
                in_reference_index_file=REF_INDEX_FILE,
                in_reference_dict_file=REF_DICT_FILE
        }
        call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzip2ndGATKGVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_merged_vcf_file=gatkJointGenotyper2nd.joint_genotyped_vcf,
                in_vg_container=VG_CONTAINER
        }
    }
    if (DRAGEN_MODE) {
        call runDragenJointGenotyper as dragenJointGenotyper2nd {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_gvcf_file_maternal=maternalCallWorkflow.output_vcf,
                in_gvcf_file_paternal=paternalCallWorkflow.output_vcf,
                in_gvcf_file_proband=probandMapCallWorkflow2ndIteration.output_vcf,
                in_dragen_ref_index_name=DRAGEN_REF_INDEX_NAME,
                in_udp_data_dir=UDPBINFO_PATH,
                in_helix_username=HELIX_USERNAME
        }
    }
    File output_2nd_joint_genotyped_vcf = select_first([bgzip2ndGATKGVCF.output_merged_vcf, dragenJointGenotyper2nd.dragen_joint_genotyped_vcf])
    File output_2nd_joint_genotyped_vcf_index = select_first([bgzip2ndGATKGVCF.output_merged_vcf_index, dragenJointGenotyper2nd.dragen_joint_genotyped_vcf_index])
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION) {
        call vgMultiMapCallWorkflow.normalizeVCF as normalizeCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_bgzip_vcf_file=output_2nd_joint_genotyped_vcf,
        }
        call vgMultiMapCallWorkflow.snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_normalized_vcf_file=normalizeCohortVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
        }
    }
    if (!SNPEFF_ANNOTATION) {
        File final_vcf_output = output_2nd_joint_genotyped_vcf
    }
    
    output {
        File output_cohort_vcf = select_first([snpEffAnnotateCohortVCF.output_snpeff_annotated_vcf, final_vcf_output])
        File? output_maternal_bam = maternalMapWorkflow.output_bam
        File? output_maternal_bam_index = maternalMapWorkflow.output_bam_index
        File? output_paternal_bam = paternalMapWorkflow.output_bam
        File? output_paternal_bam_index = paternalMapWorkflow.output_bam_index
        File output_final_proband_vcf = probandMapCallWorkflow2ndIteration.output_vcf
        File? output_final_proband_bam = probandMapCallWorkflow2ndIteration.output_bam
        File? output_final_proband_bam_index = probandMapCallWorkflow2ndIteration.output_bam_index
    }
}

########################
### TASK DEFINITIONS ###
########################

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
        
        gatk IndexFeatureFile \
            -F ${in_gvcf_file_maternal} \
        && gatk IndexFeatureFile \
            -F ${in_gvcf_file_paternal} \
        && gatk IndexFeatureFile \
            -F ${in_gvcf_file_proband} \
        && gatk CombineGVCFs \
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
        memory: 100 + " GB"
        cpu: 32
        docker: "broadinstitute/gatk:4.1.1.0"
    }
}

task runDragenJointGenotyper {
    input {
        String in_sample_name
        File in_gvcf_file_maternal
        File in_gvcf_file_paternal
        File in_gvcf_file_proband
        String in_dragen_ref_index_name
        String in_udp_data_dir
        String in_helix_username
    }

    String maternal_gvcf_file_name = basename(in_gvcf_file_maternal)
    String paternal_gvcf_file_name = basename(in_gvcf_file_paternal)
    String proband_gvcf_file_name = basename(in_gvcf_file_proband)

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
        
        ## Copy input GVCFs into directory that Dragen can access
        UDP_DATA_DIR_PATH="~{in_udp_data_dir}/usr/~{in_helix_username}"
        mkdir -p /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        cp ~{in_gvcf_file_maternal} ~{in_gvcf_file_paternal} ~{in_gvcf_file_proband} /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/ && \
        JOINT_GENOTYPE_DRAGEN_WORK_DIR="/staging/~{in_helix_username}/output_cohort_joint_call_~{in_sample_name}" && \
        TMP_DIR="/staging/~{in_helix_username}/tmp" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "mkdir -p ${TMP_DIR}" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "dragen -f -r /staging/~{in_dragen_ref_index_name} --enable-joint-genotyping true --intermediate-results-dir ${TMP_DIR} --output-directory ${JOINT_GENOTYPE_DRAGEN_WORK_DIR} --output-file-prefix cohort_joint_genotyped_~{in_sample_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{maternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{paternal_gvcf_file_name} --variant /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_cohort_gvcfs/~{proband_gvcf_file_name}" && \
        mkdir /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && chmod ug+rw -R /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "cp -R ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/. /staging/helix/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper" && \
        ssh -t ~{in_helix_username}@helix.nih.gov ssh 165.112.174.51 "rm -fr ${JOINT_GENOTYPE_DRAGEN_WORK_DIR}/" && \
        mv /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_dragen_joint_genotyper ~{in_sample_name}_dragen_joint_genotyper && \
        rm -fr /data/${UDP_DATA_DIR_PATH}/~{in_sample_name}_surjected_bams/
    >>>
    output {
        File dragen_joint_genotyped_vcf = "~{in_sample_name}_dragen_joint_genotyper/cohort_joint_genotyped_~{in_sample_name}.vcf.gz"
        File dragen_joint_genotyped_vcf_index = "~{in_sample_name}_dragen_joint_genotyper/cohort_joint_genotyped_~{in_sample_name}.vcf.gz.tbi"
    }
    runtime {
        memory: 50 + " GB"
    }
}

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

        while read -r contig; do
            if [[ ~{filter_parents} == true && ${contig} == "MT" ]]; then
                bcftools view -O z -r "${contig}" "-s ~{in_maternal_sample_name}" ~{joint_genotyped_vcf} > "${contig}.vcf.gz"
            else
                bcftools view -O z -r "${contig}" ${SAMPLE_FILTER_STRING} ~{joint_genotyped_vcf} > "${contig}.vcf.gz"
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
        whatshap phase \
            --reference ~{in_reference_file} \
            --indels \
            --ped ~{in_ped_file} \
            ${GENMAP_OPTION_STRING} \
            -o ~{in_cohort_sample_name}_cohort_~{in_contig}.phased.vcf \
            ~{joint_genotyped_vcf} ~{in_proband_bam} ~{in_maternal_bam} ~{in_paternal_bam} \
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

