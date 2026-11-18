version 1.0

import "./giraffe_and_deepvariant.wdl" as gdv_wf
import "./aardvark_evaluation.wdl" as aardvark_wf
import "../tasks/plot_aardvarks.wdl" as plot

workflow GraphAcceptanceTest {

    meta {
        description: "### Graph acceptance testing workflow \n Workflow for comparing the calling accuracy of two pangenome graphs, a \"baseline\" (the known-good release) and a \"candidate\" (the graph under test), using short-read HG002 sequencing data and truth set. \n\n The same reads are mapped to each graph with vg Giraffe and called with pangenome-aware [DeepVariant](https://github.com/google/deepvariant), using that same graph as the DeepVariant pangenome. Both call sets are then evaluated against one truth set with [Aardvark](https://github.com/PacificBiosciences/aardvark). The vg version, DeepVariant version, and all settings are the same for both runs, so only the graph differs. \n\n Both runs use the same linear reference FASTA, so both graphs must contain the same reference paths."
    }

    parameter_meta {
        BASELINE_GBZ_FILE: "Path to .gbz file for the baseline graph, the known-good release to compare against."
        BASELINE_DIST_FILE: "(OPTIONAL) Path to .dist index for the baseline graph. Built if not provided."
        BASELINE_MIN_FILE: "(OPTIONAL) Path to .min index for the baseline graph. Built if not provided."
        BASELINE_HAPL_FILE: "(OPTIONAL) Path to .hapl file for haplotype sampling of the baseline graph. Built if not provided."

        CANDIDATE_GBZ_FILE: "Path to .gbz file for the candidate graph, the graph under test. Must contain the same reference paths as the baseline graph."
        CANDIDATE_DIST_FILE: "(OPTIONAL) Path to .dist index for the candidate graph. Built if not provided."
        CANDIDATE_MIN_FILE: "(OPTIONAL) Path to .min index for the candidate graph. Built if not provided."
        CANDIDATE_HAPL_FILE: "(OPTIONAL) Path to .hapl file for haplotype sampling of the candidate graph. Built if not provided."

        INPUT_READ_FILE: "Input sample reads, as an interleaved paired fastq.gz"

        SAMPLE_NAME: "The sample name. Default is 'HG002'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        INTERLEAVED_READS: "Are paired reads interleaved in a single FASTQ? Default is 'true'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 million."
        TRUTH_VCF: "Path to .vcf.gz of truth calls to evaluate both runs against"
        TRUTH_VCF_INDEX: "(OPTIONAL) Tabix index for TRUTH_VCF. Made if not provided."
        REGIONS_BED: "BED of confident regions for Aardvark to evaluate in"
        STRATIFICATION_ARCHIVE: "(OPTIONAL) tar.gz of a GIAB-style stratification folder to break the Aardvark results down by"
        REFERENCE_FILE: "FASTA of the linear reference, used by both runs and by Aardvark. Contig names must match the graph reference paths with REFERENCE_PREFIX removed."
        REFERENCE_INDEX_FILE: "(OPTIONAL) .fai index for REFERENCE_FILE. Made if not provided."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (e.g. 'GRCh38#0#')"
        HAPLOTYPE_SAMPLING: "Whether or not to use haplotype sampling, before running giraffe and for the DeepVariant pangenome. Default is 'true'."

        BAN_SAMPLE: "(OPTIONAL) Sample whose haplotypes are excluded from haplotype sampling, in both mapping and DeepVariant. Set this to the sample being called (e.g. 'HG002') if it is in the graphs, so the truth sample doesn't help call itself."
        DV_PANGENOME_REF_NAME: "(OPTIONAL) Reference sample name in the graphs for pangenome-aware DeepVariant (e.g. 'GRCh38')"

        DV_PANGENOME_HAPLOTYPE_NUMBER: "Number of haplotypes to sample for the DeepVariant pangenome. Default is 32."
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB: "(OPTIONAL) Shared memory size, in GB, for pangenome-aware DeepVariant"
        DV_MODEL_TYPE: "Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA."
        DV_USE_GPUS: "Should DeepVariant use GPUs for calling variants? Default is 'true'."
        DV_NO_GPU_DOCKER: "Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be a pangenome-aware DeepVariant image."
        DV_GPU_DOCKER: "Container image to use when running DeepVariant for steps that benefit from GPUs. Must be a pangenome-aware DeepVariant image."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        INDEX_MINIMIZER_MEM: "Memory, in GB, to use when making the minimizer index. Default is 320."
        HAPLOTYPE_INDEXING_MEM: "Memory, in GB, to use for haplotype sampling indexing tasks. Default is 200."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
        EVAL_CORES: "Number of cores to use when evaluating variant calls. Default is 8."
        EVAL_MEM: "Memory, in GB, to use when evaluating variant calls. Default is 60."
        VG_DOCKER: "Container image to use when running vg, for both runs"
    }

    # ? means the input is optional
    input {
        File BASELINE_GBZ_FILE
        File? BASELINE_DIST_FILE
        File? BASELINE_MIN_FILE
        File? BASELINE_HAPL_FILE

        File CANDIDATE_GBZ_FILE
        File? CANDIDATE_DIST_FILE
        File? CANDIDATE_MIN_FILE
        File? CANDIDATE_HAPL_FILE

        File INPUT_READ_FILE

        String SAMPLE_NAME = "HG002"
        Boolean PAIRED_READS = true
        Boolean INTERLEAVED_READS = true
        Int READS_PER_CHUNK = 20000000
        File TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File REGIONS_BED
        File? STRATIFICATION_ARCHIVE
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        String REFERENCE_PREFIX = ""
        Boolean HAPLOTYPE_SAMPLING = true

        String? BAN_SAMPLE
        String? DV_PANGENOME_REF_NAME

        Int DV_PANGENOME_HAPLOTYPE_NUMBER = 32
        Int? DV_PANGENOME_SHARED_MEMORY_SIZE_GB
        String DV_MODEL_TYPE = "WGS"
        Boolean DV_USE_GPUS = true
        String? DV_NO_GPU_DOCKER
        String? DV_GPU_DOCKER
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Int INDEX_MINIMIZER_MEM = 320
        Int HAPLOTYPE_INDEXING_MEM = 200
        Int CALL_CORES = 8
        Int CALL_MEM = 50
        Int EVAL_CORES = 8
        Int EVAL_MEM = 60
        String VG_DOCKER = "quay.io/vgteam/vg:v1.77.0"
    }

    ###
    # Baseline graph: map and call
    ###

    call gdv_wf.GiraffeDeepVariant as baselineRun {
        input:
        INPUT_READ_FILE_1=INPUT_READ_FILE,
        GBZ_FILE=BASELINE_GBZ_FILE,
        DIST_FILE=BASELINE_DIST_FILE,
        MIN_FILE=BASELINE_MIN_FILE,
        HAPL_FILE=BASELINE_HAPL_FILE,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_GAF=false,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        READS_PER_CHUNK=READS_PER_CHUNK,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        BAN_SAMPLE=BAN_SAMPLE,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_PANGENOME_GBZ=BASELINE_GBZ_FILE,
        DV_PANGENOME_REF_NAME=DV_PANGENOME_REF_NAME,
        DV_PANGENOME_HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        DV_PANGENOME_HAPLOTYPE_NUMBER=DV_PANGENOME_HAPLOTYPE_NUMBER,
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB=DV_PANGENOME_SHARED_MEMORY_SIZE_GB,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        MAP_CORES=MAP_CORES,
        MAP_MEM=MAP_MEM,
        INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        VG_DOCKER=VG_DOCKER
    }

    ###
    # Candidate graph: map and call, with the exact same settings
    ###

    call gdv_wf.GiraffeDeepVariant as candidateRun {
        input:
        INPUT_READ_FILE_1=INPUT_READ_FILE,
        GBZ_FILE=CANDIDATE_GBZ_FILE,
        DIST_FILE=CANDIDATE_DIST_FILE,
        MIN_FILE=CANDIDATE_MIN_FILE,
        HAPL_FILE=CANDIDATE_HAPL_FILE,
        SAMPLE_NAME=SAMPLE_NAME,
        OUTPUT_GAF=false,
        PAIRED_READS=PAIRED_READS,
        INTERLEAVED_READS=INTERLEAVED_READS,
        READS_PER_CHUNK=READS_PER_CHUNK,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REFERENCE_PREFIX=REFERENCE_PREFIX,
        HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        BAN_SAMPLE=BAN_SAMPLE,
        DV_MODEL_TYPE=DV_MODEL_TYPE,
        DV_PANGENOME_GBZ=CANDIDATE_GBZ_FILE,
        DV_PANGENOME_REF_NAME=DV_PANGENOME_REF_NAME,
        DV_PANGENOME_HAPLOTYPE_SAMPLING=HAPLOTYPE_SAMPLING,
        DV_PANGENOME_HAPLOTYPE_NUMBER=DV_PANGENOME_HAPLOTYPE_NUMBER,
        DV_PANGENOME_SHARED_MEMORY_SIZE_GB=DV_PANGENOME_SHARED_MEMORY_SIZE_GB,
        DV_USE_GPUS=DV_USE_GPUS,
        DV_NO_GPU_DOCKER=DV_NO_GPU_DOCKER,
        DV_GPU_DOCKER=DV_GPU_DOCKER,
        MAP_CORES=MAP_CORES,
        MAP_MEM=MAP_MEM,
        INDEX_MINIMIZER_MEM=INDEX_MINIMIZER_MEM,
        HAPLOTYPE_INDEXING_MEM=HAPLOTYPE_INDEXING_MEM,
        CALL_CORES=CALL_CORES,
        CALL_MEM=CALL_MEM,
        VG_DOCKER=VG_DOCKER
    }

    ###
    # Evaluate both call sets against the same truth set
    ###

    call aardvark_wf.AardvarkEvaluation as baselineEval {
        input:
        QUERY_VCF=baselineRun.output_vcf,
        QUERY_VCF_INDEX=baselineRun.output_vcf_index,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REGIONS_BED=REGIONS_BED,
        STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
        SAMPLE_NAME=SAMPLE_NAME,
        THREADS=EVAL_CORES,
        EVAL_MEM=EVAL_MEM
    }

    call aardvark_wf.AardvarkEvaluation as candidateEval {
        input:
        QUERY_VCF=candidateRun.output_vcf,
        QUERY_VCF_INDEX=candidateRun.output_vcf_index,
        TRUTH_VCF=TRUTH_VCF,
        TRUTH_VCF_INDEX=TRUTH_VCF_INDEX,
        REFERENCE_FILE=REFERENCE_FILE,
        REFERENCE_INDEX_FILE=REFERENCE_INDEX_FILE,
        REGIONS_BED=REGIONS_BED,
        STRATIFICATION_ARCHIVE=STRATIFICATION_ARCHIVE,
        SAMPLE_NAME=SAMPLE_NAME,
        THREADS=EVAL_CORES,
        EVAL_MEM=EVAL_MEM
    }
    ###
    # Plot the two Aardvark summaries side by side
    ###

    call plot.plotAardvarkComparison {
        input:
        baseline_summary=baselineEval.aardvark_summary,
        candidate_summary=candidateEval.aardvark_summary,
        baseline_name=basename(BASELINE_GBZ_FILE, ".gbz"),
        candidate_name=basename(CANDIDATE_GBZ_FILE, ".gbz")
    }

    output {
        Array[File] baseline_aardvark_all_files = baselineEval.aardvark_all_files
        Array[File] candidate_aardvark_all_files = candidateEval.aardvark_all_files
        File comparison_plot = plotAardvarkComparison.plot
    }
}