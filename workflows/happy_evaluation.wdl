version 1.0

import "../tasks/variant_evaluation.wdl" as eval
import "../tasks/bioinfo_utils.wdl" as utils

workflow HappyEvaluation {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Evaluate small variants using hap.py. More information at [https://github.com/vgteam/vg_wdl/tree/gbz#happy-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#happy-workflow)."
    }

    parameter_meta {
        VCF: "bgzipped VCF with variant calls"
        VCF_INDEX: "(Optional) If specified, use this tabix index for the VCF instead of indexing it"
        TRUTH_VCF: "bgzipped VCF with truthset"
        TRUTH_VCF_INDEX: "(Optional) If specified, use this index for the truth VCF instead of indexing it"
        REFERENCE_FILE: "Use this FASTA reference."
        REFERENCE_INDEX_FILE: "(Optional) If specified, use this .fai index instead of indexing the reference file."
        EVALUATION_REGIONS_BED: "(Optional) BED to restrict comparison against TRUTH_VCF to"
        REFERENCE_PREFIX: "(Optional) Remove this off the beginning of sequence names in the VCF"
        REMOVE_HOM_REFS: "(Optional) Should homozygous ref calls be removed? (might help if hap.py segfaults). Default 'false'."
        RUN_STANDALONE_VCFEVAL: "whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)"
    }
    
    input {
        File VCF
        File? VCF_INDEX
        File TRUTH_VCF
        File? TRUTH_VCF_INDEX
        File REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? EVALUATION_REGIONS_BED
        String REFERENCE_PREFIX = ""
        Boolean REMOVE_HOM_REFS = false
        Boolean RUN_STANDALONE_VCFEVAL = true
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
    
    if (REFERENCE_PREFIX != "") {
        # use samtools to replace the header contigs with those from our dict.
        # this allows the header to contain contigs that are not in the graph,
        # which is more general and lets CHM13-based graphs be used to call on GRCh38
        # also, strip out contig prefixes in the BAM body
        call utils.fixVCFContigNaming {
            input:
            in_vcf=VCF,
            in_prefix_to_strip=REFERENCE_PREFIX
        }
    }
    File fixed_vcf_file = select_first([fixVCFContigNaming.vcf, VCF]) 
    
    if (REMOVE_HOM_REFS) {
        # sometimes blocks of homozygous variants make hap.py segfault
        # removing them is safer
        call utils.removeHomRefs {
            input:
            in_vcf=fixed_vcf_file
        }
    }
    File vcf_file = select_first([removeHomRefs.vcf, fixed_vcf_file]) 
    
    ## index the call VCF if needed
    if (!defined(VCF_INDEX) || REFERENCE_PREFIX != "" || REMOVE_HOM_REFS) {
        call utils.indexVcf {
            input:
            in_vcf=vcf_file
        }
    }
    File vcf_index = select_first([indexVcf.vcf_index_file, VCF_INDEX])

    ## index the truth VCF if needed
    if (!defined(TRUTH_VCF_INDEX)) {
        call utils.indexVcf as indexTruthVcf{
            input:
            in_vcf=TRUTH_VCF
        }
    }
    File truth_vcf_index = select_first([TRUTH_VCF_INDEX, indexTruthVcf.vcf_index_file])
    
    if (RUN_STANDALONE_VCFEVAL) {
        # Direct vcfeval comparison makes an archive with FP and FN VCFs
        call eval.compareCalls {
            input:
            in_sample_vcf_file=vcf_file,
            in_sample_vcf_index_file=vcf_index,
            in_truth_vcf_file=TRUTH_VCF,
            in_truth_vcf_index_file=truth_vcf_index,
            in_template_archive=buildReferenceTemplate.output_template_archive,
            in_evaluation_regions_file=EVALUATION_REGIONS_BED
        }
    }
    
    # Hap.py comparison makes accuracy results stratified by SNPs and indels,
    # and also an archive with FP and FN VCFs.
    call eval.compareCallsHappy {
        input:
        in_sample_vcf_file=vcf_file,
        in_sample_vcf_index_file=vcf_index,
        in_truth_vcf_file=TRUTH_VCF,
        in_truth_vcf_index_file=truth_vcf_index,
        in_reference_file=REFERENCE_FILE,
        in_reference_index_file=reference_index_file,
        in_evaluation_regions_file=EVALUATION_REGIONS_BED
    }
    
    output {
        File? output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File? output_vcfeval_evaluation_summary = compareCalls.output_evaluation_summary
        File output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_happy_evaluation_summary = compareCallsHappy.output_evaluation_summary
    }   
}
