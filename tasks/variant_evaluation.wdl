version 1.0

task buildReferenceTemplate {
    input {
        File in_reference_file
    }
    command <<<
        set -eux -o pipefail
    
        rtg format -o template.sdf "~{in_reference_file}"
        tar -czf template.sdf.tar.gz template.sdf/
    >>>
    output {
        File output_template_archive = "template.sdf.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "realtimegenomics/rtg-tools:3.12.1"
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + 10 + " SSD" 
    }
}

task compareCalls {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_template_archive
        File? in_evaluation_regions_file
        Int in_disk = 3 * round(size(in_sample_vcf_file, "G") + size(in_truth_vcf_file, "G")) + 20
        Int in_mem = 16
        Int in_cores = 8
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Set up template; we assume it drops a "template.sdf"
        tar -xf "~{in_template_archive}"
    
        rtg vcfeval \
            --baseline truth.vcf.gz \
            --calls sample.vcf.gz \
            ~{"--evaluation-regions=" + in_evaluation_regions_file} \
            --template template.sdf \
            --threads ~{in_cores} \
            --output vcfeval_results
            
        tar -czf vcfeval_results.tar.gz vcfeval_results/
    >>>
    output {
        File output_evaluation_summary = "vcfeval_results/summary.txt"
        File output_evaluation_archive = "vcfeval_results.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "realtimegenomics/rtg-tools:3.12.1"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        memory: in_mem + " GB"
    }
}

task compareCallsHappy {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_reference_file
        File in_reference_index_file
        File? in_evaluation_regions_file
        Int in_disk = 3 * round(size(in_sample_vcf_file, "G") + size(in_truth_vcf_file, "G") + size(in_reference_file, "G")) + 20
        Int in_mem = 16
        Int in_cores = 8
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        
        mkdir happy_results
   
        /opt/hap.py/bin/hap.py \
            truth.vcf.gz \
            sample.vcf.gz \
            ~{"-f " + in_evaluation_regions_file} \
            --reference reference.fa \
            --threads ~{in_cores} \
            --engine=vcfeval \
            -o happy_results/eval
    
        tar -czf happy_results.tar.gz happy_results/
    >>>
    output {
        File output_evaluation_summary = "happy_results/eval.summary.csv"
        File output_evaluation_archive = "happy_results.tar.gz"
    }
    runtime {
        preemptible: 2
        docker: "jmcdani20/hap.py:v0.3.12"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        memory: in_mem + " GB"
    }
}
