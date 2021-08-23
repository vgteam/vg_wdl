version 1.0

### vg_trio_giraffe_deeptrio_workflow.wdl ###
## Author: Charles Markello
## Description: Trio-backed VG mapping and variant calling workflow for mother-father-child trio datasets using giraffe and deeptrio platforms.
## Reference: https://github.com/vgteam/vg/wiki

import "./vg_multi_map.wdl" as vgMultiMapWorkflow
import "./vg_deeptrio_calling_workflow.wdl" as vgDeepTrioCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow

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
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.31.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 10000000                      # Number of reads contained in each mapping chunk (20000000 for wgs)
        Int CHUNK_BASES = 50000000                          # Number of bases to chunk .gam alignment files for variant calling
        Int OVERLAP = 2000                                  # Number of overlapping bases between each .gam chunk
        File? PATH_LIST_FILE                                # (OPTIONAL) Text file where each line is a path name in the XG index
        File XG_FILE                                        # Path to .xg index file
        File GBWT_FILE                                      # Path to .gbwt index file
        File GGBWT_FILE                                     # Path to .gg index file
        File DIST_FILE                                      # Path to .dist index file
        File MIN_FILE                                       # Path to .min index file
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
        Array[String]+ CONTIGS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        File REF_FASTA_GZ
        File PED_FILE
        File EAGLE_DATA
        File? GEN_MAP_FILES
        File? DEEPTRIO_CHILD_MODEL
        File? DEEPTRIO_PARENT_MODEL
        File? DEEPVAR_MODEL
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
            GBWT_FILE=GBWT_FILE,
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
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
            GBWT_FILE=GBWT_FILE,
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
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
            GBWT_FILE=GBWT_FILE,
            GGBWT_FILE=GGBWT_FILE,
            DIST_FILE=DIST_FILE,
            MIN_FILE=MIN_FILE,
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
            MATERNAL_BAM_FILE=maternalMapWorkflow.output_bam,
            MATERNAL_BAM_FILE_INDEX=maternalMapWorkflow.output_bam_index,
            PATERNAL_BAM_FILE=paternalMapWorkflow.output_bam,
            PATERNAL_BAM_FILE_INDEX=paternalMapWorkflow.output_bam_index,
            CHILD_BAM_FILE=probandMapWorkflow.output_bam,
            CHILD_BAM_FILE_INDEX=probandMapWorkflow.output_bam_index,
            SAMPLE_NAME_CHILD=SAMPLE_NAME_PROBAND,
            SAMPLE_NAME_MATERNAL=SAMPLE_NAME_MATERNAL,
            SAMPLE_NAME_PATERNAL=SAMPLE_NAME_PATERNAL,
            CONTIGS=CONTIGS,
            REF_FILE=REF_FILE,
            REF_INDEX_FILE=REF_INDEX_FILE,
            REF_DICT_FILE=REF_DICT_FILE,
            DEEPTRIO_CHILD_MODEL=DEEPTRIO_CHILD_MODEL,
            DEEPTRIO_PARENT_MODEL=DEEPTRIO_PARENT_MODEL,
            DEEPVAR_MODEL=DEEPVAR_MODEL,
            MAP_CORES=MAP_CORES,
            MAP_DISK=MAP_DISK,
            MAP_MEM=MAP_MEM,
            VGCALL_CORES=VGCALL_CORES,
            VGCALL_DISK=VGCALL_DISK,
            VGCALL_MEM=VGCALL_MEM
    }
    
    # Collect parental indel-realigned BAM contig lists
    Array[File] maternal_bams_by_contig = select_first([firstTrioCallWorkflow.output_maternal_indelrealigned_bams])
    Array[File] maternal_bam_indexes_by_contig = select_first([firstTrioCallWorkflow.output_maternal_indelrealigned_bam_indexes])
    Array[File] paternal_bams_by_contig = select_first([firstTrioCallWorkflow.output_paternal_indelrealigned_bams])
    Array[File] paternal_bam_indexes_by_contig = select_first([firstTrioCallWorkflow.output_paternal_indelrealigned_bam_indexes])
    Array[File] child_bams_by_contig = firstTrioCallWorkflow.output_child_indelrealigned_bams
    Array[File] child_bam_indexes_by_contig = firstTrioCallWorkflow.output_child_indelrealigned_bam_indexes
    
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
        call runPrepPhasing {
            input:
                in_eagle_data=EAGLE_DATA,
                in_contigs=CONTIGS
        }
        Array[Pair[File,File]] maternal_bam_index_by_contigs_pair = zip(maternal_bams_by_contig, maternal_bam_indexes_by_contig)
        Array[Pair[File,File]] paternal_bam_index_by_contigs_pair = zip(paternal_bams_by_contig, paternal_bam_indexes_by_contig)
        Array[Pair[File,File]] child_bam_index_by_contigs_pair = zip(child_bams_by_contig, child_bam_indexes_by_contig)
        Array[Pair[Pair[File,File],Pair[Pair[File,File],Pair[File,File]]]] trio_bam_index_by_contigs_pair = zip(child_bam_index_by_contigs_pair,zip(maternal_bam_index_by_contigs_pair,paternal_bam_index_by_contigs_pair))
        scatter (contig_pair in zip(CONTIGS, zip(splitJointGenotypedVCF.contig_vcfs, trio_bam_index_by_contigs_pair))) {
            if (defined(runPrepPhasing.eagle_data[contig_pair.left])) {
                call runEaglePhasing {
                    input:
                        in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                        joint_genotyped_vcf=contig_pair.right.left,
                        in_eagle_bcf=runPrepPhasing.eagle_data[contig_pair.left][0],
                        in_eagle_bcf_index=runPrepPhasing.eagle_data[contig_pair.left][1],
                        in_contig=contig_pair.left,
                        in_vgcall_cores=VGCALL_CORES,
                        in_vgcall_disk=VGCALL_DISK,
                        in_vgcall_mem=VGCALL_MEM
                }
            }
            File pre_whatshap_vcf_file = select_first([runEaglePhasing.phased_cohort_vcf, contig_pair.right.left])
            call runWhatsHapPhasing {
                input:
                    in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                    joint_genotyped_vcf=pre_whatshap_vcf_file,
                    in_maternal_bam=contig_pair.right.right.right.left.left,
                    in_maternal_bam_index=contig_pair.right.right.right.left.right,
                    in_paternal_bam=contig_pair.right.right.right.right.left,
                    in_paternal_bam_index=contig_pair.right.right.right.right.right,
                    in_proband_bam=contig_pair.right.right.left.left,
                    in_proband_bam_index=contig_pair.right.right.left.right,
                    in_ped_file=PED_FILE,
                    in_genetic_map=GEN_MAP_FILES,
                    in_contig=contig_pair.left,
                    in_reference_file=REF_FILE,
                    in_reference_index_file=REF_INDEX_FILE,
                    in_reference_dict_file=REF_DICT_FILE
            }
        }
        call concatClippedVCFChunks as concatCohortPhasedVCFs {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_clipped_vcf_chunk_files=runWhatsHapPhasing.phased_cohort_vcf,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call bgzipMergedVCF as bgzipCohortPhasedVCF {
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
                GBWT_FILE=constructGraphIndexWorkflow.gbwt,
                GGBWT_FILE=constructGraphIndexWorkflow.ggbwt,
                DIST_FILE=constructGraphIndexWorkflow.dist,
                MIN_FILE=constructGraphIndexWorkflow.min,
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
                CHILD_BAM_FILE=secondIterationSiblingMap.output_bam,
                CHILD_BAM_FILE_INDEX=secondIterationSiblingMap.output_bam_index,
                SAMPLE_NAME_CHILD=read_pair_set.right,
                SAMPLE_NAME_MATERNAL=SAMPLE_NAME_MATERNAL,
                SAMPLE_NAME_PATERNAL=SAMPLE_NAME_PATERNAL,
                CONTIGS=CONTIGS,
                REF_FILE=REF_FILE,
                REF_INDEX_FILE=REF_INDEX_FILE,
                REF_DICT_FILE=REF_DICT_FILE,
                DEEPTRIO_CHILD_MODEL=DEEPTRIO_CHILD_MODEL,
                DEEPTRIO_PARENT_MODEL=DEEPTRIO_PARENT_MODEL,
                DEEPVAR_MODEL=DEEPVAR_MODEL,
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
    Array[File] sibling_bam_list = mergeSiblingIndelRealignedBams.merged_indel_realigned_bam_file
    Array[File] sibling_bam_index_list = mergeSiblingIndelRealignedBams.merged_indel_realigned_bam_file_index
    Array[File] gvcf_files_siblings = secondTrioCallWorkflow.output_child_gvcf
    
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
        call normalizeVCF as normalizeCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_bgzip_vcf_file=output_2nd_joint_genotyped_vcf,
                in_vgcall_cores=VGCALL_CORES,
                in_vgcall_disk=VGCALL_DISK,
                in_vgcall_mem=VGCALL_MEM
        }
        call snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
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
        Array[File] output_sibling_bam_list = sibling_bam_list
        Array[File] output_sibling_bam_index_list = sibling_bam_index_list
    }
}

########################
### TASK DEFINITIONS ###
########################
task extractPathNames {
    input {
        File in_xg_file
        String in_vg_container
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt
    }
    output {
        File output_path_list = "path_list.txt"
    }
    runtime {
        memory: "50 GB"
        disks: "local-disk 50 SSD"
        docker: in_vg_container
    }
}

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
        File? in_gvcf_file_maternal
        File? in_gvcf_file_paternal
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
        Array[File]? in_alignment_bam_chunk_files
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

task runPrepPhasing {
    input {
        File in_eagle_data
        Array[String]+ in_contigs
    }
    
    command <<<
        set -exu -o pipefail
        tar -xf ~{in_eagle_data}
        json_string=""
        while read -r contig; do
            if [[ ~{in_contigs[0]} == *"chr"* ]]; then
                if [ ${contig} == "chrX" ]; then
                    json_string="${json_string}\"${contig}\":[\"eagle_data_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_${contig}.filtered.eagle2-phased.bcf\", \"eagle_data_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_${contig}.filtered.eagle2-phased.bcf.csi\"],"
                elif [[ ! ${contig} =~ ^(chrY|chrM)$ ]]; then
                    json_string="${json_string}\"${contig}\":[\"eagle_data_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_${contig}.filtered.shapeit2-duohmm-phased.bcf\", \"eagle_data_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_${contig}.filtered.shapeit2-duohmm-phased.bcf.csi\"],"
                fi
            else
                if [[ ! ${contig} =~ ^(Y|MT|ABOlocus)$ ]]; then
                    json_string="${json_string}\"${contig}\":[\"eagle_data/ALL.chr${contig}.phase3_integrated.20130502.genotypes.bcf\", \"eagle_data/ALL.chr${contig}.phase3_integrated.20130502.genotypes.bcf.csi\"],"
                fi
            fi
        done < "~{write_lines(in_contigs)}"
        json_string_output="${json_string%?}"
        json_string_output="{${json_string_output}}"
        echo "${json_string_output}"
    >>>
    
    output {
        Map[String, Array[File]] eagle_data = read_json(stdout())
    }
    
    runtime {
        memory: "10 GB"
        disks: "local-disk 70 SSD"
        docker: "ubuntu:latest"
    }
}

task runEaglePhasing {
    input {
        String in_cohort_sample_name
        File joint_genotyped_vcf
        File in_eagle_bcf
        File in_eagle_bcf_index
        String in_contig
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
    }
    
    command <<<
        set -exu -o pipefail
        ln -s ~{in_eagle_bcf} input_eagle.bcf
        ln -s ~{in_eagle_bcf_index} input_eagle.bcf.csi
        tabix -p vcf ~{joint_genotyped_vcf}
        gen_map_file="/usr/src/app/genetic_map_hg38_withX.txt.gz"
        if [[ ~{in_contig} != *"chr"* ]]; then
            gen_map_file="/usr/src/app/genetic_map_hg19_withX.txt.gz"
        fi
        /usr/src/app/eagle \
        --outputUnphased \
        --geneticMapFile ${gen_map_file} \
        --outPrefix "~{in_cohort_sample_name}_cohort_~{in_contig}.eagle_phased" \
        --numThreads ~{in_vgcall_cores} \
        --vcfRef input_eagle.bcf \
        --vcfTarget ~{joint_genotyped_vcf} \
        --chrom ~{in_contig}
    >>>
    
    output {
        File phased_cohort_vcf = "~{in_cohort_sample_name}_cohort_~{in_contig}.eagle_phased.vcf.gz"
    }
    
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/cmarkello/eagle"
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
        ln -s ~{in_maternal_bam} input_m_bam_file.bam
        ln -s ~{in_maternal_bam_index} input_m_bam_file.bam.bai
        ln -s ~{in_paternal_bam} input_p_bam_file.bam
        ln -s ~{in_paternal_bam_index} input_p_bam_file.bam.bai
        ln -s ~{in_proband_bam} input_c_bam_file.bam
        ln -s ~{in_proband_bam_index} input_c_bam_file.bam.bai
        
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
            ~{joint_genotyped_vcf} input_c_bam_file.bam input_m_bam_file.bam input_p_bam_file.bam \
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

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        for vcf_file in ${sep=" " in_clipped_vcf_chunk_files} ; do
            bcftools index "$vcf_file"
        done
        bcftools concat -a ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort - > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        time: 60
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
        Int in_vgcall_disk
        Int in_vgcall_mem
    }

    # TODO:
    #   If GVCF in in_merged_vcf_file then output_vcf_extension="gvcf" else output_vcf_extension="vcf"
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

        bgzip -c ${in_merged_vcf_file} > ${in_sample_name}_merged.vcf.gz && \
        tabix -f -p vcf ${in_sample_name}_merged.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}_merged.vcf.gz.tbi"
    }
    runtime {
        time: 30
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: in_vg_container
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        bcftools norm -m-both --threads ~{in_vgcall_cores} -o ~{in_sample_name}.unrolled.vcf ~{in_bgzip_vcf_file}
    >>>
    output {
        File output_normalized_vcf = "~{in_sample_name}.unrolled.vcf"
    }
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File? in_snpeff_database
        Int in_vgcall_cores
        Int in_vgcall_disk
        Int in_vgcall_mem
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

        unzip ~{in_snpeff_database}
        snpEff -Xmx40g -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data GRCh37.75 ~{in_normalized_vcf_file} > ~{in_sample_name}.snpeff.unrolled.vcf
    >>>
    output {
        File output_snpeff_annotated_vcf = "~{in_sample_name}.snpeff.unrolled.vcf"
    }
    runtime {
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:4.3.1t--2"
    }
}


