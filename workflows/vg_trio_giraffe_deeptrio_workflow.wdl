version 1.0

### vg_trio_giraffe_deeptrio_workflow.wdl ###
## Author: Charles Markello
## Description: Trio-backed VG mapping and variant calling workflow for mother-father-child trio datasets using giraffe and deeptrio platforms.
## Reference: https://github.com/vgteam/vg/wiki

import "./vg_multi_map.wdl" as vgMultiMapWorkflow
import "./vg_deeptrio_calling_workflow.wdl" as vgDeepTrioCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow
import "./vg_indel_realign.wdl" as vgIndelRealignmentWorkflow

workflow vgTrioPipeline {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        description: "Core VG mapping and variant calling workflow for pedigree datasets using giraffe and deeptrio platforms."
    }
    input {
        File MATERNAL_INPUT_READ_FILE_1                     # Input maternal 1st read pair fastq.gz
        File MATERNAL_INPUT_READ_FILE_2                     # Input maternal 2nd read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_1                     # Input paternal 1st read pair fastq.gz
        File PATERNAL_INPUT_READ_FILE_2                     # Input paternal 2nd read pair fastq.gz
        Array[File]+ SIBLING_INPUT_READ_FILE_1_LIST         # Input 1st read pair list where the proband file is listed first followed by sibling fastq.gz
        Array[File]+ SIBLING_INPUT_READ_FILE_2_LIST         # Input 2nd read pair list where the proband file is listed first followed by sibling fastq.gz
        Array[String]+ SAMPLE_NAME_SIBLING_LIST             # Sample name for siblings where the proband file is listed first followed by the siblings
        String SAMPLE_NAME_MATERNAL                         # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                         # Sample name for the father
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.28.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 10000000                      # Number of reads contained in each mapping chunk (20000000 for wgs)
        Int CHUNK_BASES = 50000000                          # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                                  # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File GBWT_FILE                                      # Path to .gbwt index file
        File GGBWT_FILE                                     # Path to .gg index file
        File DIST_FILE                                      # Path to .dist index file
        File MIN_FILE                                       # Path to .min index file
        File SNARLS_FILE                                    # Path to .snarls index file
        File REF_FILE                                       # Path to .fa cannonical reference fasta (only grch37/hg19 currently supported)
        File REF_INDEX_FILE                                 # Path to .fai index of the REF_FILE fasta reference
        File REF_DICT_FILE                                  # Path to .dict file of the REF_FILE fasta reference
        File? SNPEFF_DATABASE                               # Path to snpeff database .zip file for snpEff annotation functionality.
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
        Array[String]+ CONTIGS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]
        File REF_FASTA_GZ
        File PED_FILE
        File? GEN_MAP_FILES
        String GRAPH_NAME
        Boolean GIRAFFE_INDEXES = true                  # Set to 'true' to construct the GBWT index which incorporates haplotype information into the graph.
        Boolean USE_DECOYS = true                       # Set to 'true' to include decoy contigs from the FASTA reference into the graph reference.
        Boolean SNPEFF_ANNOTATION = true                # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        Boolean CLEANUP_FILES = true                    # Set to 'false' to turn off intermediate file cleanup.
        Boolean ABRA_REALIGN = false                    # Set to 'true' to use GATK IndelRealigner instead of ABRA2 for indel realignment
        String DECOY_REGEX = ">GL\|>NC_007605\|>hs37d5\|>hs38d1_decoys\|>chrEBV\|>chrUn\|>chr\([1-2][1-9]\|[1-9]\|Y\)_" # grep regular expression string that is used to extract decoy contig ids. USE_DECOYS must be set to 'true'
    }
    
    File PROBAND_INPUT_READ_FILE_1 = SIBLING_INPUT_READ_FILE_1_LIST[0]  # Input proband 1st read pair fastq.gz
    File PROBAND_INPUT_READ_FILE_2 = SIBLING_INPUT_READ_FILE_2_LIST[0]  # Input proband 2nd read pair fastq.gz
    String SAMPLE_NAME_PROBAND = SAMPLE_NAME_SIBLING_LIST[0]            # Sample name for the proband
    
    # Extract path names and path lengths from xg file if PATH_LIST_FILE input not provided
    if (!defined(PATH_LIST_FILE)) {
        call extractPathNames {
            input:
                in_xg_file=XG_FILE,
                in_vg_container=VG_CONTAINER
        }
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractPathNames.output_path_list])
        
    #######################################################
    ## Run mapping and variant calling workflows on Trio ##
    #######################################################
    call vgMultiMapWorkflow.vgMultiMap as maternalMapWorkflow {
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
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
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
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true
    }
    call vgMultiMapWorkflow.vgMultiMap as paternalMapWorkflow {
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
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true
    }
    call vgMultiMapWorkflow.vgMultiMap as probandMapWorkflow {
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
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
            SNARLS_FILE=SNARLS_FILE,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            SPLIT_READ_CORES=SPLIT_READ_CORES,
            SPLIT_READ_DISK=SPLIT_READ_DISK,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true
    }
    
    #########################
    ## Run trio genotyping ##
    #########################
    call vgDeepTrioCallWorkflow.vgDeepTrioCall as firstTrioCallWorkflow {
        input:
            MATERNAL_BAM_FILE=maternalMapWorkflow.merged_bam_file_output,
            MATERNAL_BAM_FILE_INDEX=maternalMapWorkflow.merged_bam_file_index_output,
            PATERNAL_BAM_FILE=paternalMapWorkflow.merged_bam_file_output,
            PATERNAL_BAM_FILE_INDEX=paternalMapWorkflow.merged_bam_file_index_output,
            CHILD_BAM_FILE=probandMapWorkflow.merged_bam_file_output,
            CHILD_BAM_FILE_INDEX=probandMapWorkflow.merged_bam_file_index_output,
            CHILD_SAMPLE_NAME=SAMPLE_NAME_PROBAND,
            MATERNAL_SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
            PATERNAL_SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
            PATH_LIST_FILE=pipeline_path_list_file,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM
    }
    
    # Collect parental indel-realigned BAM contig lists
    Array[File] maternal_bams_by_contig = firstTrioCallWorkflow.maternal_indelrealigned_bams_by_contig
    Array[File] maternal_bam_indexes_by_contig = firstTrioCallWorkflow.maternal_indelrealigned_bam_indexes_by_contig
    Array[File] paternal_bams_by_contig = firstTrioCallWorkflow.paternal_indelrealigned_bams_by_contig
    Array[File] paternal_bam_indexes_by_contig = firstTrioCallWorkflow.paternal_indelrealigned_bam_indexes_by_contig
    
    call mergeIndelRealignedBAMs as mergeMaternalIndelRealignedBams{
        input:  
            in_sample_name=SAMPLE_NAME_MATERNAL,
            in_alignment_bam_chunk_files=maternal_bams_by_contig
    }
    call mergeIndelRealignedBAMs as mergePaternalIndelRealignedBams{
        input:  
            in_sample_name=SAMPLE_NAME_PATERNAL,
            in_alignment_bam_chunk_files=paternal_bams_by_contig
    }
    
    call runDeepVariantJointGenotyper as deepVarJointGenotyper1st {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_gvcf_file_maternal=firstTrioCallWorkflow.output_maternal_gvcf,
            in_gvcf_file_paternal=firstTrioCallWorkflow.output_paternal_gvcf,
            in_gvcf_files_siblings=[firstTrioCallWorkflow.output_child_gvcf],
            in_vgcall_cores=VGCALL_CORES,
            in_vgcall_disk=VGCALL_DISK,
            in_vgcall_mem=VGCALL_MEM
    }
    
    #####################################
    ## Run parental graph construction ##
    #####################################
    File output_joint_genotyped_vcf = deepVarJointGenotyper1st.joint_genotyped_vcf
    File output_joint_genotyped_vcf_index = deepVarJointGenotyper1st.joint_genotyped_vcf_index
    if (GIRAFFE_INDEXES) {
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
                    in_maternal_bam=maternalMapCallWorkflow.output_bam,
                    in_maternal_bam_index=maternalMapCallWorkflow.output_bam_index,
                    in_paternal_bam=paternalMapCallWorkflow.output_bam,
                    in_paternal_bam_index=paternalMapCallWorkflow.output_bam_index,
                    in_proband_bam=probandMapCallWorkflow.output_bam,
                    in_proband_bam_index=probandMapCallWorkflow.output_bam_index,
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
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_clipped_vcf_chunk_files=runWhatsHapPhasing.phased_cohort_vcf,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call vgMultiMapCallWorkflow.bgzipMergedVCF as bgzipCohortPhasedVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_merged_vcf_file=concatCohortPhasedVCFs.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    File cohort_vcf_output = select_first([bgzipCohortPhasedVCF.output_merged_vcf, output_joint_genotyped_vcf])
    File cohort_vcf_index_output = select_first([bgzipCohortPhasedVCF.output_merged_vcf_index, output_joint_genotyped_vcf_index])
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
            giraffe_indexes=GIRAFFE_INDEXES,
            use_decoys=USE_DECOYS,
            decoy_regex=DECOY_REGEX,
            vg_docker=VG_CONTAINER
    }
    
    ######################################################################################################
    ## Run mapping and variant calling workflow of proband and sibling reads against the parental graph ##
    ######################################################################################################
    Array[Pair[File,File]] read_pair_files_list = zip(SIBLING_INPUT_READ_FILE_1_LIST, SIBLING_INPUT_READ_FILE_2_LIST)
    scatter (read_pair_set in zip(read_pair_files_list, SAMPLE_NAME_SIBLING_LIST)) {
        Pair[File,File] read_pair_files = read_pair_set.left
        call vgMultiMapWorkflow.vgMultiMap as secondIterationSiblingMap {
            input:
                INPUT_READ_FILE_1=read_pair_files.left,
                INPUT_READ_FILE_2=read_pair_files.right,
                SAMPLE_NAME=read_pair_set.right,
                VG_CONTAINER=VG_CONTAINER,
                READS_PER_CHUNK=READS_PER_CHUNK,
                PATH_LIST_FILE=PATH_LIST_FILE,
                XG_FILE=constructGraphIndexWorkflow.xg,
                GCSA_FILE=constructGraphIndexWorkflow.gcsa,
                GCSA_LCP_FILE=constructGraphIndexWorkflow.gcsa_lcp,
                GBWT_FILE=constructGraphIndexWorkflow.gbwt,
                GGBWT_FILE=constructGraphIndexWorkflow.ggbwt,
                DIST_FILE=constructGraphIndexWorkflow.dist,
                MIN_FILE=constructGraphIndexWorkflow.min,
                SNARLS_FILE=constructGraphIndexWorkflow.snarls,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                SPLIT_READ_CORES=SPLIT_READ_CORES,
                SPLIT_READ_DISK=SPLIT_READ_DISK,
                MAP_CORES=MAP_CORES,
                MAP_DISK=MAP_DISK,
                MAP_MEM=MAP_MEM,
                MAPPER="GIRAFFE",
                CLEANUP_FILES=CLEANUP_FILES,
                SURJECT_MODE=true
        }
        
        call vgDeepTrioCallWorkflow.vgDeepTrioCall as secondTrioCallWorkflow {
            input:
                MATERNAL_BAM_CONTIG_LIST=maternal_bams_by_contig,
                MATERNAL_BAM_INDEX_CONTIG_LIST=maternal_bam_indexes_by_contig,
                PATERNAL_BAM_CONTIG_LIST=paternal_bams_by_contig,
                PATERNAL_BAM_INDEX_CONTIG_LIST=paternal_bam_indexes_by_contig,
                CHILD_BAM_FILE=secondIterationSiblingMap.merged_bam_file_output,
                CHILD_BAM_FILE_INDEX=secondIterationSiblingMap.merged_bam_file_index_output,
                CHILD_SAMPLE_NAME=read_pair_set.right,
                MATERNAL_SAMPLE_NAME=SAMPLE_NAME_MATERNAL,
                PATERNAL_SAMPLE_NAME=SAMPLE_NAME_PATERNAL,
                PATH_LIST_FILE=pipeline_path_list_file,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                MAP_CORES=MAP_CORES,
                MAP_DISK=MAP_DISK,
                MAP_MEM=MAP_MEM,
                VGCALL_CORES=VGCALL_CORES,
                VGCALL_DISK=VGCALL_DISK,
                VGCALL_MEM=VGCALL_MEM
        }
        call mergeIndelRealignedBAMs as mergeSiblingIndelRealignedBams{
            input:
                in_sample_name=read_pair_set.right,
                in_alignment_bam_chunk_files=secondTrioCallWorkflow.output_child_indelrealigned_bams
        }
        
    }
    #######################################################
    ## Run 2nd trio joint genotyping on new proband GVCF ##
    #######################################################
    Array[File] output_sibling_bam_list = select_all(mergeSiblingIndelRealignedBams.merged_indel_realigned_bam_file)
    Array[File] output_sibling_bam_index_list = select_all(mergeSiblingIndelRealignedBams.merged_indel_realigned_bam_file_index)
    Array[File] gvcf_files_siblings = select_all(secondTrioCallWorkflow.output_child_gvcf)
    
    call runDeepVariantJointGenotyper as deepVarJointGenotyper2nd {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_gvcf_file_maternal=firstTrioCallWorkflow.output_maternal_gvcf,
            in_gvcf_file_paternal=firstTrioCallWorkflow.output_paternal_gvcf,
            in_gvcf_files_siblings=gvcf_files_siblings,
            in_vgcall_cores=VGCALL_CORES,
            in_vgcall_disk=VGCALL_DISK,
            in_vgcall_mem=VGCALL_MEM
    }
    
    File output_2nd_joint_genotyped_vcf = deepVarJointGenotyper2nd.joint_genotyped_vcf
    File output_2nd_joint_genotyped_vcf_index = deepVarJointGenotyper2nd.joint_genotyped_vcf_index
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION && defined(SNPEFF_DATABASE)) {
        call vgMultiMapCallWorkflow.normalizeVCF as normalizeCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_bgzip_vcf_file=output_2nd_joint_genotyped_vcf,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call vgMultiMapCallWorkflow.snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_normalized_vcf_file=normalizeCohortVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
    }
    if (!SNPEFF_ANNOTATION) {
        File final_vcf_output = output_2nd_joint_genotyped_vcf
    }
    
    output {
        File output_cohort_vcf = select_first([snpEffAnnotateCohortVCF.output_snpeff_annotated_vcf, final_vcf_output])
        File output_maternal_bam = mergeMaternalIndelRealignedBams.merged_indel_realigned_bam_file
        File output_maternal_bam_index = mergeMaternalIndelRealignedBams.merged_indel_realigned_bam_file_index
        File output_paternal_bam = mergePaternalIndelRealignedBams.merged_indel_realigned_bam_file
        File output_paternal_bam_index = mergePaternalIndelRealignedBams.merged_indel_realigned_bam_file_index
        Array[File] output_gvcf_files_siblings = gvcf_files_siblings
        Array[File] output_sibling_bam_list = output_sibling_bam_list
        Array[File] output_sibling_bam_index_list = output_sibling_bam_index_list
    }
}

########################
### TASK DEFINITIONS ###
########################
task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        while IFS=$'\t' read -ra path_list_line; do
            path_name="${path_list_line[0]}"
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${path_name} \
              -o ~{in_sample_name}.${path_name}.bam \
            && samtools index \
              ~{in_sample_name}.${path_name}.bam
        done < ~{in_path_list_file}
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
        Array[Pair[File, File]] bams_and_indexes_by_contig = zip(bam_contig_files, bam_contig_files_index)
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runDeepVariantJointGenotyper {
    input {
        String in_sample_name
        File in_gvcf_file_maternal
        File in_gvcf_file_paternal
        Array[File]+ in_gvcf_files_siblings
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    command <<<
        set -exu -o pipefail
        
        tabix -f -p vcf ~{in_gvcf_file_maternal}
        tabix -f -p vcf ~{in_gvcf_file_paternal}
        for sibling_gvcf_file in ~{sep=" " in_gvcf_files_siblings} ; do
            tabix -f -p vcf "${sibling_gvcf_file}"
        done
        
        /usr/local/bin/glnexus_cli \
        --config DeepVariant_unfiltered \
        --threads ~{in_vgcall_cores} \
        ~{in_gvcf_file_maternal} \
        ~{in_gvcf_file_paternal} \
        ~{sep=" " in_gvcf_files_siblings} \
        | bcftools view - | bgzip -c > ~{in_sample_name}_cohort.jointgenotyped.vcf.gz \
        && tabix -f -p vcf ~{in_sample_name}_cohort.jointgenotyped.vcf.gz
    >>>
    output {
        File joint_genotyped_vcf = "~{in_sample_name}_cohort.jointgenotyped.vcf.gz"
        File joint_genotyped_vcf_index = "~{in_sample_name}_cohort.jointgenotyped.vcf.gz.tbi"
    }
    runtime {
        memory: in_vgcall_mem + " GB"
        cpu: in_vgcall_cores
        docker: "quay.io/mlin/glnexus:v1.2.7"
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
        ln -s ~{joint_genotyped_vcf} input_vcf_file.vcf.gz
        ln -s ~{joint_genotyped_vcf_index} input_vcf_file.vcf.gz.tbi
        while read -r contig; do
            if [[ ~{filter_parents} == true && (${contig} == "MT" || ${contig} == chr"M") ]]; then
                bcftools view -O z -r "${contig}" -s ~{in_maternal_sample_name} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            elif [[ ~{filter_parents} == true && (${contig} == "Y" || ${contig} == chr"Y") ]]; then
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
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task mergeIndelRealignedBAMs {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
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
          -f -p -c --threads "$(nproc --all)" \
          ~{in_sample_name}_merged.indel_realigned.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.indel_realigned.bam
    >>>
    output {
        File merged_indel_realigned_bam_file = "~{in_sample_name}_merged.indel_realigned.bam"
        File merged_indel_realigned_bam_file_index = "~{in_sample_name}_merged.indel_realigned.bam.bai"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        disks: "local-disk 100 SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
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
        File? in_genetic_map
        String in_contig
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
    }
    
    Boolean genetic_map_available = defined(in_genetic_map)
    
    command <<<
        set -exu -o pipefail
        
        GENMAP_OPTION_STRING=""
        if [ ~{genetic_map_available} == true ]; then
            tar -xvf ~{in_genetic_map}
            if [[ ~{in_contig} == "Y" || ~{in_contig} == "MT" || ~{in_contig} == "ABOlocus" ]]; then
                GENMAP_OPTION_STRING=""
            elif [ ~{in_contig} == "X" ]; then
                GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chrX_nonPAR_combined_b37.txt --chromosome X"
            else
                GENMAP_OPTION_STRING="--genmap genetic_map_GRCh37/genetic_map_chr~{in_contig}_combined_b37.txt --chromosome ~{in_contig}"
            fi
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
        docker: "quay.io/biocontainers/whatshap@sha256:cf82de1173a35a0cb063469a602eff2e8999b4cfc0f0ee9cef0dbaedafa5ab6c"
    }
}

