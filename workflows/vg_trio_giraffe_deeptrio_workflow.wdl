version 1.0

### vg_trio_giraffe_deeptrio_workflow.wdl ###
## Author: Charles Markello
## Description: Trio-backed VG mapping and variant calling workflow for mother-father-child trio datasets using giraffe and deeptrio platforms.
## Reference: https://github.com/vgteam/vg/wiki

import "./vg_multi_map.wdl" as vgMultiMapWorkflow
import "./vg_deeptrio_calling_workflow.wdl" as vgDeepTrioCallWorkflow
import "./vg_construct_and_index.wdl" as vgConstructWorkflow
#import "https://raw.githubusercontent.com/vgteam/vg_wdl/master/workflows/vg_multi_map.wdl" as vgMultiMapWorkflow
#import "https://raw.githubusercontent.com/vgteam/vg_wdl/master/workflows/vg_deeptrio_calling_workflow.wdl" as vgDeepTrioCallWorkflow
#import "https://raw.githubusercontent.com/vgteam/vg_wdl/master/workflows/vg_construct_and_index.wdl" as vgConstructWorkflow

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
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.34.0"   # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        Int READS_PER_CHUNK = 200000000                      # Number of reads contained in each mapping chunk (20000000 for wgs)
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
        Array[String]+ CONTIGS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        File PED_FILE
        File EAGLE_DATA
        File? GEN_MAP_FILES
        File? DEEPTRIO_CHILD_MODEL
        File? DEEPTRIO_PARENT_MODEL
        File? DEEPVAR_MODEL
        String GRAPH_NAME
        Boolean SMALL_RESOURCES = false                 # Set to 'true' to use small resources for tiny test dataset
        Boolean GIRAFFE_INDEXES = true                  # Set to 'true' to construct the GBWT index which incorporates haplotype information into the graph.
        Boolean USE_DECOYS = true                       # Set to 'true' to include decoy contigs from the FASTA reference into the graph reference.
        Boolean SNPEFF_ANNOTATION = true                # Set to 'true' to run snpEff annotation on the joint genotyped VCF.
        Boolean CLEANUP_FILES = false                   # Set to 'true' to turn on intermediate file cleanup.
        Boolean ABRA_REALIGN = false                    # Set to 'true' to use GATK IndelRealigner instead of ABRA2 for indel realignment.
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
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true,
            SMALL_RESOURCES=SMALL_RESOURCES
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
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true,
            SMALL_RESOURCES=SMALL_RESOURCES
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
            MAPPER="GIRAFFE",
            CLEANUP_FILES=CLEANUP_FILES,
            SURJECT_MODE=true,
            SMALL_RESOURCES=SMALL_RESOURCES
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
            SMALL_RESOURCES=SMALL_RESOURCES
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
            in_alignment_bam_chunk_files=maternal_bams_by_contig,
            in_small_resources=SMALL_RESOURCES
    }
    call mergeIndelRealignedBAMs as mergePaternalIndelRealignedBams{
        input:  
            in_sample_name=SAMPLE_NAME_PATERNAL,
            in_alignment_bam_chunk_files=paternal_bams_by_contig,
            in_small_resources=SMALL_RESOURCES
    }
    
    call runDeepVariantJointGenotyper as deepVarJointGenotyper1st {
        input:
            in_sample_name=SAMPLE_NAME_PROBAND,
            in_gvcf_file_maternal=firstTrioCallWorkflow.output_maternal_gvcf,
            in_gvcf_file_paternal=firstTrioCallWorkflow.output_paternal_gvcf,
            in_gvcf_files_siblings=[firstTrioCallWorkflow.output_child_gvcf],
            in_small_resources=SMALL_RESOURCES
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
                filter_parents=false,
                in_small_resources=SMALL_RESOURCES
        }
        call runPrepPhasing {
            input:
                in_eagle_data=EAGLE_DATA
        }
        Array[Pair[File,File]] eagle_vcf_and_index_pairs = zip(runPrepPhasing.eagle_data, runPrepPhasing.eagle_data_index)
        
        #TODO
        #  1) map eagle contigs with eagle vcf/index pair array index
        #  2) use output for eagle_vcf_contig_map index
        #  3) if( defined(eagle_vcf_contig_map[contig_pair.left] )
        #  4) in_eagle_bcf = eagle_vcf_and_index_pairs[eagle_vcf_contig_map[contig_pair.left]].left
        #  5) in_eagle_bcf_index = eagle_vcf_and_index_pairs[eagle_vcf_contig_map[contig_pair.left]].right
        call runMakeContigMAP {
            input:
                in_eagle_contigs=runPrepPhasing.eagle_contigs
        }
        Array[Pair[File,File]] maternal_bam_index_by_contigs_pair = zip(maternal_bams_by_contig, maternal_bam_indexes_by_contig)
        Array[Pair[File,File]] paternal_bam_index_by_contigs_pair = zip(paternal_bams_by_contig, paternal_bam_indexes_by_contig)
        Array[Pair[File,File]] child_bam_index_by_contigs_pair = zip(child_bams_by_contig, child_bam_indexes_by_contig)
        Array[Pair[Pair[File,File],Pair[Pair[File,File],Pair[File,File]]]] trio_bam_index_by_contigs_pair = zip(child_bam_index_by_contigs_pair,zip(maternal_bam_index_by_contigs_pair,paternal_bam_index_by_contigs_pair))
        scatter (contig_pair in zip(splitJointGenotypedVCF.contig_vcfs_contig_list, zip(splitJointGenotypedVCF.contig_vcfs, trio_bam_index_by_contigs_pair))) {
            if (contig_pair.left != "chrY" && contig_pair.left != "chrM" && contig_pair.left != "Y" && contig_pair.left != "M") { 
                Int eagle_vcf_contig_index = runMakeContigMAP.eagle_vcf_contig_map[contig_pair.left]
                call runEaglePhasing {
                    input:
                        in_cohort_sample_name=SAMPLE_NAME_PROBAND,
                        joint_genotyped_vcf=contig_pair.right.left,
                        in_eagle_bcf=eagle_vcf_and_index_pairs[eagle_vcf_contig_index].left,
                        in_eagle_bcf_index=eagle_vcf_and_index_pairs[eagle_vcf_contig_index].right,
                        in_contig=contig_pair.left,
                        in_small_resources=SMALL_RESOURCES
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
                    in_small_resources=SMALL_RESOURCES
            }
        }
        call concatClippedVCFChunks as concatCohortPhasedVCFs {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_clipped_vcf_chunk_files=runWhatsHapPhasing.phased_cohort_vcf,
                in_small_resources=SMALL_RESOURCES
        }
        call bgzipMergedVCF as bgzipCohortPhasedVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_merged_vcf_file=concatCohortPhasedVCFs.output_merged_vcf,
                in_vg_container=VG_CONTAINER,
                in_small_resources=SMALL_RESOURCES
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
            filter_parents=true,
            in_small_resources=SMALL_RESOURCES
    }
    call vgConstructWorkflow.vg_construct_and_index as constructGraphIndexWorkflow {
        input:
            graph_name=GRAPH_NAME,
            ref_fasta_gz=REF_FILE,
            contigs=splitPhasedVCF.contig_vcfs_contig_list,
            contigs_vcf_gz=splitPhasedVCF.contig_vcfs,
            giraffe_indexes=GIRAFFE_INDEXES,
            use_decoys=USE_DECOYS,
            decoy_regex=DECOY_REGEX,
            vg_docker=VG_CONTAINER,
            in_small_resources=SMALL_RESOURCES
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
                MAPPER="GIRAFFE",
                CLEANUP_FILES=CLEANUP_FILES,
                SURJECT_MODE=true,
                SMALL_RESOURCES=SMALL_RESOURCES
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
                SMALL_RESOURCES=SMALL_RESOURCES
        }
        call mergeIndelRealignedBAMs as mergeSiblingIndelRealignedBams{
            input:
                in_sample_name=read_pair_set.right,
                in_alignment_bam_chunk_files=secondTrioCallWorkflow.output_child_indelrealigned_bams,
                in_small_resources=SMALL_RESOURCES
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
            in_small_resources=SMALL_RESOURCES
    }
    
    File output_2nd_joint_genotyped_vcf = deepVarJointGenotyper2nd.joint_genotyped_vcf
    File output_2nd_joint_genotyped_vcf_index = deepVarJointGenotyper2nd.joint_genotyped_vcf_index
    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION && defined(SNPEFF_DATABASE)) {
        call normalizeVCF as normalizeCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_bgzip_vcf_file=output_2nd_joint_genotyped_vcf,
                in_small_resources=SMALL_RESOURCES
        }
        call snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
            input:
                in_sample_name=SAMPLE_NAME_PROBAND,
                in_normalized_vcf_file=normalizeCohortVCF.output_normalized_vcf,
                in_snpeff_database=SNPEFF_DATABASE,
                in_small_resources=SMALL_RESOURCES
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

task runDeepVariantJointGenotyper {
    input {
        String in_sample_name
        File? in_gvcf_file_maternal
        File? in_gvcf_file_paternal
        Array[File]+ in_gvcf_files_siblings
        Boolean in_small_resources
    }

    Int in_vgcall_cores = if in_small_resources then 6 else 6
    Int in_vgcall_disk = if in_small_resources then 1 else 50
    String in_vgcall_mem = if in_small_resources then "1" else "50"

    command <<<
        set -exu -o pipefail
        
        tabix -f -p vcf ~{in_gvcf_file_maternal}
        tabix -f -p vcf ~{in_gvcf_file_paternal}
        for sibling_gvcf_file in ~{sep=" " in_gvcf_files_siblings} ; do
            tabix -f -p vcf "${sibling_gvcf_file}"
        done
        mkdir -p tmp
        --temp-dir ./tmp
        /usr/local/bin/glnexus_cli \
        --mem-gbytes ~{in_vgcall_mem} \
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
        preemptible: 1
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
        Boolean in_small_resources
    }
    
    Int in_cores = if in_small_resources then 6 else 6
    Int in_disk = if in_small_resources then 1 else 25
    String in_mem = if in_small_resources then "1" else "50"

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
                bcftools view --threads ~{in_cores} -O z -r "${contig}" -s ~{in_maternal_sample_name} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            elif [[ ~{filter_parents} == true && (${contig} == "Y" || ${contig} == chr"Y") ]]; then
                bcftools view --threads ~{in_cores} -O z -r "${contig}" -s ~{in_paternal_sample_name} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            else
                bcftools view --threads ~{in_cores} -O z -r "${contig}" ${SAMPLE_FILTER_STRING} input_vcf_file.vcf.gz > "${contig}.vcf.gz"
            fi
        done < ~{write_lines(contigs)}
        rm input_vcf_file.vcf.gz input_vcf_file.vcf.gz.tbi
        ls *.vcf.gz | sed s/.vcf.gz//g
    >>>
    output {
        Array[File]+ contig_vcfs = glob("*.vcf.gz")
        Array[String]+ contig_vcfs_contig_list = read_lines(stdout())
    }
    runtime {
        preemptible: 3
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task mergeIndelRealignedBAMs {
    input {
        String in_sample_name
        Array[File]? in_alignment_bam_chunk_files
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 8 else 10
    Int in_disk = if in_small_resources then 1 else 100
    String in_mem = if in_small_resources then "4" else "50"

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
        preemptible: 1
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runPrepPhasing {
    input {
        File in_eagle_data
    }
    
    command <<<
        set -exu -o pipefail
        mkdir eagle_data
        tar -xzf ~{in_eagle_data} -C eagle_data --strip-components 1
        for FILENAME in eagle_data/*.bcf; do
            if [[ ~{in_eagle_data} == *"grch38"* ]]; then
                echo $FILENAME | rev | cut -f1 -d'/' | rev | sed s/CCDG_14151_B01_GRM_WGS_2020-08-05_//g | sed s/.filtered.shapeit2-duohmm-phased.bcf$//g | sed s/.filtered.eagle2-phased.bcf$//g
            else
                echo $FILENAME | rev | cut -f1 -d'/' | rev | sed s/ALL.chr//g | sed s/.phase3_integrated.20130502.genotypes.bcf$//g
            fi
        done
    >>>
    
    output {
        Array[File]+ eagle_data = glob("eagle_data/*.bcf")
        Array[File]+ eagle_data_index = glob("eagle_data/*.bcf.csi")
        Array[String]+ eagle_contigs = read_lines(stdout())
    }
    
    runtime {
        preemptible: 1
        memory: "1 GB"
        disks: "local-disk 140 SSD"
        docker: "ubuntu:latest"
    }
}

task runMakeContigMAP {
    input {
        Array[String]+ in_eagle_contigs
    }
    Int in_eagle_contigs_length = length(in_eagle_contigs)
    command <<<
        python <<CODE
        contigs = "~{sep='\t' in_eagle_contigs}"
        contig_list = contigs.split()
        for i in range(~{in_eagle_contigs_length}):
            eagle_contig_name = contig_list[i]
            print("{}\t{}".format(eagle_contig_name,i))
        CODE
    >>>
    output {
        Map[String, Int] eagle_vcf_contig_map = read_map(stdout())
    }
    runtime {
        preemptible: 1
        docker: "python:3.9-slim-bullseye"
    }
}

task runEaglePhasing {
    input {
        String in_cohort_sample_name
        File joint_genotyped_vcf
        File in_eagle_bcf
        File in_eagle_bcf_index
        String in_contig
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 8
    Int in_disk = if in_small_resources then 50 else 50
    String in_mem = if in_small_resources then "10" else "20"
    
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
        --numThreads ~{in_cores} \
        --vcfRef input_eagle.bcf \
        --vcfTarget ~{joint_genotyped_vcf} \
        --chrom ~{in_contig}
    >>>
    
    output {
        File phased_cohort_vcf = "~{in_cohort_sample_name}_cohort_~{in_contig}.eagle_phased.vcf.gz"
    }
    
    runtime {
        preemptible: 1
        time: 300
        cpu: in_cores
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
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
        Boolean in_small_resources
    }

    Int in_cores = if in_small_resources then 2 else 8
    Int in_disk = if in_small_resources then 50 else 50
    String in_mem = if in_small_resources then "20" else "20"
    
    Boolean genetic_map_available = defined(in_genetic_map)
    
    
    command <<<
        set -exu -o pipefail
        ln -s ~{in_maternal_bam} input_m_bam_file.bam
        ln -s ~{in_maternal_bam_index} input_m_bam_file.bam.bai
        ln -s ~{in_paternal_bam} input_p_bam_file.bam
        ln -s ~{in_paternal_bam_index} input_p_bam_file.bam.bai
        ln -s ~{in_proband_bam} input_c_bam_file.bam
        ln -s ~{in_proband_bam_index} input_c_bam_file.bam.bai
        gzip -dc ~{in_reference_file} > ref.fna
        ln -f -s ~{in_reference_index_file} ref.fna.fai
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
            --reference ref.fna \
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
        preemptible: 1
        time: 300
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/whatshap@sha256:cf82de1173a35a0cb063469a602eff2e8999b4cfc0f0ee9cef0dbaedafa5ab6c"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Boolean in_small_resources
    }

    Int in_vgcall_disk = if in_small_resources then 1 else 25
    String in_vgcall_mem = if in_small_resources then "1" else "50"

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
        mkdir -p tmp
        for vcf_file in ${sep=" " in_clipped_vcf_chunk_files} ; do
            bcftools index "$vcf_file"
        done
        bcftools concat -a ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort --temp-dir ./tmp - > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        preemptible: 1
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
        Boolean in_small_resources
    }

    Int in_vgcall_disk = if in_small_resources then 1 else 25
    String in_vgcall_mem = if in_small_resources then "1" else "50"

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
        preemptible: 1
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
        Boolean in_small_resources
    }

    Int in_vgcall_cores = if in_small_resources then 6 else 6
    Int in_vgcall_disk = if in_small_resources then 1 else 25
    String in_vgcall_mem = if in_small_resources then "1" else "50"

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
        preemptible: 1
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
        Boolean in_small_resources
    }

    Int in_vgcall_cores = if in_small_resources then 6 else 6
    Int in_vgcall_disk = if in_small_resources then 10 else 25
    String in_vgcall_mem = if in_small_resources then "10" else "50"

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
        database_ref="GRCh38.99"
        if [[ "~{in_snpeff_database}" != *"GRCh38"* ]]; then
            database_ref="GRCh37.75"
        fi
        snpEff -Xmx40g -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data ${database_ref} ~{in_normalized_vcf_file} > ~{in_sample_name}.snpeff.unrolled.vcf
        bgzip ~{in_sample_name}.snpeff.unrolled.vcf
    >>>
    output {
        File output_snpeff_annotated_vcf = "~{in_sample_name}.snpeff.unrolled.vcf.gz"
    }
    runtime {
        preemptible: 1
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
    }
}


