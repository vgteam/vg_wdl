version 1.0

import "../tasks/variant_evaluation.wdl" as eval
import "../tasks/bioinfo_utils.wdl" as utils

workflow HappyEvaluation {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Evaluate small variants using hap.py"
    }

    parameter_meta {
        VCF: "bgzipped VCF with variant calls"
        VCF_INDEX: "If specified, use this tabix index for the VCF instead of indexing it"
        TRUTH_VCF: "bgzipped VCF with truthset"
        TRUTH_VCF_INDEX: "Tabix index for the truth VCF"
        REFERENCE_FILE: "Use this FASTA reference."
        EVALUATION_REGIONS_BED: "BED to restrict comparison against TRUTH_VCF to"
        REFERENCE_INDEX_FILE: "If specified, use this .fai index instead of indexing the reference file."
    }
    
    input {
        File VCF
        File? VCF_INDEX
        File TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File? EVALUATION_REGIONS_BED
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
    }
    
    # To evaluate the VCF we need a template of the reference
    call eval.buildReferenceTemplate {
        input:
        in_reference_file=REFERENCE_FILE
    }

    ## index the reference FASTA if needed
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference {
            input:
            in_reference_file=REFERENCE_FILE
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])

    ## index the call VCF if needed
    if (!defined(VCF_INDEX)) {
        call utils.indexVcf {
            input:
            in_vcf=VCF
        }
    }
    File vcf_index = select_first([VCF_INDEX, indexVcf.vcf_index_file])

    ## index the truth VCF if needed
    if (!defined(TRUTH_VCF_INDEX)) {
        call utils.indexVcf as indexTruthVcf{
            input:
            in_vcf=TRUTH_VCF
        }
    }
    File truth_vcf_index = select_first([TRUTH_VCF_INDEX, indexTruthVcf.vcf_index_file])

    # Direct vcfeval comparison makes an archive with FP and FN VCFs
    call eval.compareCalls {
        input:
        in_sample_vcf_file=VCF,
        in_sample_vcf_index_file=vcf_index,
        in_truth_vcf_file=TRUTH_VCF,
        in_truth_vcf_index_file=truth_vcf_index,
        in_template_archive=buildReferenceTemplate.output_template_archive,
        in_evaluation_regions_file=EVALUATION_REGIONS_BED
    }
    
    # Hap.py comparison makes accuracy results stratified by SNPs and indels
    call eval.compareCallsHappy {
        input:
        in_sample_vcf_file=VCF,
        in_sample_vcf_index_file=vcf_index,
        in_truth_vcf_file=TRUTH_VCF,
        in_truth_vcf_index_file=truth_vcf_index,
        in_reference_file=REFERENCE_FILE,
        in_reference_index_file=reference_index_file,
        in_evaluation_regions_file=EVALUATION_REGIONS_BED
    }
    
    output {
        File output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File output_vcfeval_evaluation_summary = compareCalls.output_evaluation_summary
        File output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_happy_evaluation_summary = compareCallsHappy.output_evaluation_summary
    }   
}
